#include "bwmapper.h"
#include <time.h>

char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

char revert[256]  = {[0 ... 255] = 0, [0] = '@', [1] = 'A', [2] = 'C', [3] = 'G', [4] = 'N', [5] = 'T'};

char rcode[256] = {[0 ... 255] = 0,
                   ['a'] = 'T', ['c'] = 'G', ['g'] = 'C', ['t'] = 'A', ['u'] = 'A', ['n'] = 'N',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N' };

int visited = 0;

void SIGSEGV_handler(int sig) {
   void *array[10];
   size_t size;

   // get void*'s for all entries on the stack
   size = backtrace(array, 10);

   // print out all the frames to stderr
   fprintf(stderr, "Error: signal %d:\n", sig);
   backtrace_symbols_fd(array, size, STDERR_FILENO);
   exit(1);
}


int main(int argc, char * argv[]) {

   // TODO: parametrize
   int opt_reverse  = 1;
   int opt_verbose  = 1;
   int tau = 1;

   // Backtrace handler
   //   signal(SIGSEGV, SIGSEGV_handler); 

   if (argc < 3) {
      fprintf(stderr, "usage: bwmapper {query <file.fastq> | index} <genome file>\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1],"query") == 0) {
      char * filename = malloc(strlen(argv[3])+5);
      strcpy(filename, argv[3]);

      // Open index file.
      strcpy(filename+strlen(argv[3]), ".fmi");
      int fd = open(filename, O_RDONLY);
      if (fd == -1) {
         fprintf(stderr, "error opening '%s': %s\n", filename, strerror(errno));
         exit(EXIT_FAILURE);
      }

      // Parse query filename.
      FILE * queryfile = fopen(argv[2], "r");
      if (queryfile == NULL) {
         fprintf(stderr, "error: could not open file.\n");
         return EXIT_FAILURE;
      }

      // Read file.
      fprintf(stderr, "reading query file...\n");
      seqstack_t * seqs = read_file(queryfile, opt_reverse, opt_verbose); // TODO: set reverse and verbose options.
      if (seqs == NULL) {
         return EXIT_FAILURE;
      }

      // Load index.
      int mflags = MAP_PRIVATE | MAP_POPULATE;
      long idxsize = lseek(fd, 0, SEEK_END);
      lseek(fd, 0, SEEK_SET);
      long * indexp = mmap(NULL, idxsize, PROT_READ, mflags, fd, 0);
      if (indexp == NULL) {
         fprintf(stderr, "error opening index file: %s.\n", strerror(errno));
         return EXIT_FAILURE;
      }
      
      // Read FM index format.
      index_t index;
      if(format_FMindex(indexp, &index))
         return EXIT_FAILURE;

      // Read chromosome index.
      strcpy(filename+strlen(argv[3]), ".index");
      chr_t * chr = read_CHRindex(filename);
      if (chr == NULL) return EXIT_FAILURE;

      // Hitmap is too slow...
      clock_t tstart = clock();
      hitmap(0, index, chr, seqs);
      // Let's try something faster, tau=0 seeds.
      //fastseed(seqs, index, chr, KMER_SIZE);
      
      fprintf(stderr, "query time: %ldus\n", ((clock()-tstart)*1000000)/CLOCKS_PER_SEC);
      
   }
   else if (strcmp(argv[1],"index") == 0) {
      write_index(argv[2]);
   }
   else 
      fprintf(stderr, "usage: bwmapper {query <seq> | index} <genome file>\n");
   
   return 0;
}


/*********************/
/** query functions **/
/*********************/

/*********************/
/** index functions **/
/*********************/


chr_t *
read_CHRindex
(
 char * filename
)
{
   // Files
   FILE * input = fopen(filename,"r");
   if (input == NULL) {
      fprintf(stderr, "error in 'read_CHRindex' (fopen): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   int chrcount = 0;
   int structsize = CHRSTR_SIZE;

   long *  start = malloc(structsize*sizeof(long));
   char ** names = malloc(structsize*sizeof(char*));

   // File read vars
   ssize_t rlen;
   size_t  sz = BUFFER_SIZE;
   char  * buffer = malloc(BUFFER_SIZE);
   int lineno = 0;

   // Read chromosome entries.
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      lineno++;
      char *name;
      int i = 0;
      while (buffer[i] != '\t' && buffer[i] != '\n') i++;
      if (buffer[i] == '\n') {
         fprintf(stderr, "illegal format in %s (line %d) - ignoring chromosome: %s\n", filename, lineno, buffer);
         continue;
      }
      
      if (buffer[rlen-1] == '\n') buffer[rlen-1] = 0;

      // Realloc stacks if necessary.
      if (chrcount >= structsize) {
         structsize *= 2;
         start = realloc(start, structsize*sizeof(long));
         names = realloc(names, structsize*sizeof(char*));
         if (start == NULL || names == NULL) {
            fprintf(stderr, "error in 'read_CHRindex' (realloc): %s\n", strerror(errno));
            exit(EXIT_FAILURE);
         }
      }

      // Save chromosome start position.
      buffer[i] = 0;
      start[chrcount] = atol(buffer);
      // Save chromosome name.
      name = buffer + i + 1;
      names[chrcount] = malloc(strlen(name)+1);
      strcpy(names[chrcount], name);
      // Inc.
      chrcount++;
   }

   // Return chromosome index structure.
   chr_t * chrindex = malloc(sizeof(chr_t));
   chrindex->nchr = chrcount;
   chrindex->start= start;
   chrindex->name = names;

   fclose(input);
   return chrindex;
}


int
format_FMindex
(
 long    * index,
 index_t * fmindex
)
{
   // index structure:
   // DATA:
   //   C[NUM_BASES]  (8-byte)
   //   gsize         (8-byte)
   //   genome[gsize] (1-byte)
   //   POS[gsize]    (8-byte)
   //   ooc(1): pos:size:val[0]:val[1]:val[2]:...:val[pos-1]          (8-byte)
   //   ...
   //   occ(NUM_BASES): pos:size:val[0]:val[1]:val[2]:...:val[pos-1]  (8-byte)
   long * start = index;
   long pointer_offset = 0;

   // C values.
   fmindex->c = start;
   pointer_offset += NUM_BASES;

   // Genome size.
   long gsize = fmindex->gsize = start[pointer_offset++];

   // Genome bases. (8-bit data).
   char * genome = (char *)(index + pointer_offset);
   fmindex->genome = genome;
   
   // Go back to (8-byte data).
   pointer_offset = 0;
   start = (long *) (genome + gsize);
   
   // Reverse BWT (pos values).
   fmindex->pos = start + pointer_offset;
   pointer_offset += gsize;

   // Occurrences.
   list_t * occ = fmindex->occ = malloc(NUM_BASES * sizeof(list_t));
   if (occ == NULL)
      return -1;

   for (int i = 0; i < NUM_BASES; i++) {
      occ[i].max = start[pointer_offset++];
      occ[i].val = start + pointer_offset;
      pointer_offset += occ[i].max;
   }

   return 0;
}

ssize_t
write_index
(
 char * filename
)
{
   long gsize;

   // Allocate structures.
   long * C, * pos;
   vstack_t ** occ = malloc(NUM_BASES * sizeof(vstack_t *));

   // Parse genome and convert to integers. (0..NBASES-1)
   fprintf(stdout, "compacting genome...\n");
   char * genome = compact_genome(filename, &gsize);
   C = compute_c(genome, gsize);
   fprintf(stdout, "done\n");

   // Output files.
   char * fmfile = malloc(strlen(filename)+5);
   strcpy(fmfile, filename);
   strcpy(fmfile + strlen(filename), ".fmi");
   int fd = open(fmfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   if (fd == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", fmfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   free(fmfile);

   // FM index
   bwt_index(genome, gsize, &pos, occ);

   // Write index.
   size_t s = 0;
   while (s < NUM_BASES*sizeof(long)) s += write(fd, C + s/sizeof(long), NUM_BASES*sizeof(long) - s);
   s = 0;
   while (s < sizeof(long)) s += write(fd, &gsize, sizeof(long));
   s = 0;
   while (s < gsize*sizeof(char)) s += write(fd, genome + s/sizeof(char), gsize*sizeof(char) - s);
   s = 0;
   while (s < gsize*sizeof(long)) s += write(fd, pos + s/sizeof(long), gsize*sizeof(long) - s);
   for (int i = 0; i < NUM_BASES; i++) {
      fprintf(stdout, "occ[%d] size = %ld (%ld bytes).\n", i, occ[i]->pos, occ[i]->pos*sizeof(long));
      s = 0;
      while (s < sizeof(long)) s += write(fd, &(occ[i]->pos), sizeof(long));
      fprintf(stdout, "occ[%d]->pos: written %ld bytes.\n", i, s);
      s = 0;
      while (s < occ[i]->pos * sizeof(long)) s += write(fd,((long*) &(occ[i]->val[0])) + s/sizeof(long), occ[i]->pos*sizeof(long) - s);
      fprintf(stdout, "occ[%d]->val: written %ld bytes.\n", i, s);
   }

   return s;
}


char *
compact_genome
(
 char * filename,
 long * genomesize
)
{
   char * fileindex = malloc(strlen(filename)+7);
   strcpy(fileindex, filename);
   strcpy(fileindex + strlen(filename), ".index");

   // Files
   FILE * input = fopen(filename,"r");
   FILE * output = fopen(fileindex, "w");

   if (input == NULL) {
      fprintf(stderr, "error in 'compact_genome' (fopen): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // File read vars
   ssize_t rlen;
   size_t  sz = BUFFER_SIZE;
   char * buffer = malloc(BUFFER_SIZE);

   // Genome storage vars
   char * genome = malloc(GENOME_SIZE);
   long gbufsize = GENOME_SIZE;
   genome[0] = '@';
   *genomesize = 1;
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      if (buffer[0] == '>') {
         buffer[rlen-1] = 0;
         int k = 0;
         while (buffer[k] != ' ' && buffer[k] != 0) k++;
         buffer[k] = 0;
         fprintf(output,"%ld\t%s\n",*genomesize,buffer+1);
      }
      else {
         if (gbufsize < *genomesize + rlen) {
            while (gbufsize < *genomesize + rlen) gbufsize *= 2;
            genome = realloc(genome, gbufsize);
            if (genome == NULL) {
               fprintf(stderr, "error in 'compact_genome' (realloc): %s\n", strerror(errno));
               exit(EXIT_FAILURE);
            }
         }
         int newline = buffer[rlen-1] == '\n';
         strncpy(genome + *genomesize, buffer, rlen-newline);
         *genomesize += rlen-newline;
      }
   }

   // Realloc buffer.
   genome = realloc(genome, *genomesize+1);

   // Reverse genome.
   long end = *genomesize - 1;
   char tmp;
   for (int i = 0 ; i < *genomesize/2; i++) {
      tmp = genome[i];
      genome[i] = genome[end-i];
      genome[end-i] = tmp;
   }
   genome[*genomesize] = 0;

   // Free memory.
   free(buffer);
   free(fileindex);
   
   return genome;
}


void
bwt_index
(
 char      * genome,
 long        gsize,
 long     ** pos,
 vstack_t ** occ
)
{
   // Sort prefixes (will reuse stacks[i] to compute C)
   fprintf(stdout,"suffix sort...\n");
   long * values = malloc((gsize+3)*sizeof(long));
   for (long i = 0; i < gsize; i++) values[i] = (long)translate[(int)genome[i]];
   values[gsize] = values[gsize+1] = values[gsize+2] = 0;
   long * sa = malloc(gsize*sizeof(long));
   suffixArray(values, sa, gsize, NUM_BASES-1);
   free(values);
   //   long * sa = dc3(genome);
   fprintf(stdout,"done.\n");

   // Fill compacted occ.
   vstack_t * stacks[NUM_BASES];
   for (int i = 0; i < NUM_BASES; i++) stacks[i] = new_stack(gsize/NUM_BASES);

   for (long i = 0; i < gsize; i++) push(stacks+translate[(int)genome[(sa[i] == 0 ? gsize-1 : sa[i] - 1)]], i);
   
   // Save stacks to output.
   *pos = sa;
   
   fprintf(stdout, "saving occ structures...\n");
   for (int i = 0; i < NUM_BASES; i++) {
      stacks[i] = realloc(stacks[i], sizeof(vstack_t) + stacks[i]->pos*sizeof(long));
      if (stacks[i] == NULL) {
         fprintf(stderr, "error in bwt_index (realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      stacks[i]->size = stacks[i]->pos;
      occ[i] = stacks[i];
   }
   fprintf(stdout, "done\n");
}


long *
compute_c
(
 char * genome,
 long gsize
)
{
   long * cnt = calloc(NUM_BASES, sizeof(long));
   for (long i = 0; i < gsize; i++) cnt[(int)translate[(int)genome[i]]]++;
   long * c = malloc(NUM_BASES * sizeof(long));
   c[0] = 0;
   for (int i = 1; i < NUM_BASES; i++) c[i] = c[i-1] + cnt[i-1];
   free(cnt);
   return c;
}




/*********************/
/** misc  functions **/
/*********************/


seqstack_t *
read_file
(
   FILE      * inputf,
   const int   reverse,
   const int   verbose
)
{
   // TODO: Set a maximum read size in MB, when this limit is reached return.
   // Once the read sequences are processed, read again and process, continuing
   // at the point of the file where the last 'read_file' returned.

   // Read first line of the file to guess format.
   // Store in global variable FORMAT.
   char c = fgetc(inputf);

   if (c == EOF) return NULL;
   if (ungetc(c, inputf) == EOF) {
      return NULL;
   }

   format_t format;

   switch(c) {
   case '>':
      if (verbose) fprintf(stderr, "FASTA format detected\n");
      format = FASTA;
      break;
   case '@':
      if (verbose) fprintf(stderr, "FASTQ format detected\n");
      format = FASTQ;
      break;
   default:
      if (verbose) fprintf(stderr, "raw format detected\n");
      format = RAW;
   }

   ssize_t nread;
   size_t linesz = BUFFER_SIZE;
   size_t sublinesz = BUFFER_SIZE;
   size_t tempsz = BUFFER_SIZE;
   char *line = malloc(BUFFER_SIZE * sizeof(char));
   char *subline = malloc(BUFFER_SIZE * sizeof(char));
   char *temp = malloc(BUFFER_SIZE * sizeof(char));
   if (line == NULL) {
      return NULL;
   }
   seqstack_t *seqstack = new_seqstack(SEQSTACK_SIZE);
   if (seqstack == NULL) {
      return NULL;
   }

   char *seq  = NULL;
   char *tag  = NULL;
   int lineno = 0;

   while ((nread = getline(&line, &linesz, inputf)) != -1) {
      lineno++;
      if (line[nread-1] == '\n') line[nread-1] = '\0';
      switch(format) {
      case RAW:
         seq = line;
         tag = line;
         break;
      case FASTA:
         if (line[0] == '>') {
            if((nread = getline(&subline, &sublinesz, inputf)) == -1) {
               fprintf(stderr, "incorrect FASTA format: line %d\n", lineno);
               continue;
            }
            lineno++;
            if (subline[nread-1] == '\n') subline[nread-1] = '\0';
            tag = line;
            seq = subline;
         } else continue;
         break;
      case FASTQ:
         if (line[0] == '@') {
            if((nread = getline(&subline, &sublinesz, inputf)) == -1) {
               fprintf(stderr, "incorrect FASTQ format: line %d\n", lineno);
               continue;
            }
            lineno++;
            if (subline[nread-1] == '\n') subline[nread-1] = '\0';
            tag = line;
            seq = subline;
            for (int i = 0; i < 2; i++) {
               if((nread = getline(&temp, &tempsz, inputf)) == -1) {
                  fprintf(stderr, "incorrect FASTQ format: line %d\n", lineno);
                  continue;
               }
               lineno++;
            }
         } else continue;
         break;
      }

      size_t seqlen = strlen(seq);
      if (seqlen > MAXSEQLEN) {
         fprintf(stderr, "max sequence length exceeded (%d)\n", MAXSEQLEN);
         fprintf(stderr, "offending sequence:\n%s\n", seq);
         continue;
      }
      
      
      if (seq_push(&seqstack, tag, seq, reverse)) return NULL;
   }

   free(line);
   free(subline);
   free(temp);
   return seqstack;
}
