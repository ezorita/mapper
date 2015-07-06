#include "bwmapper.h"

char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

char revert[256]  = {[0 ... 255] = 0, [0] = '@', [1] = 'A', [2] = 'C', [3] = 'G', [4] = 'N', [5] = 'T'};

char rcode[256] = {[0 ... 255] = 0,
                   ['a'] = 'T', ['c'] = 'G', ['g'] = 'C', ['t'] = 'A', ['u'] = 'A', ['n'] = 'N',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N' };

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

index_t *
load_index
(
 char * file
)
{
   // Alloc index struct.
   index_t * index = malloc(sizeof(index_t));
   
   // Read chromosome index.
   char * chr_file = malloc(strlen(file)+5);
   strcpy(chr_file, file);
   strcpy(chr_file+strlen(file), ".chr");
   index->chr = read_CHRindex(chr_file);
   if (index->chr == NULL) {
      exit(EXIT_FAILURE);
   }
   free(chr_file);

   // Open OCC file.
   char * occ_file = malloc(strlen(file)+5);
   strcpy(occ_file, file);
   strcpy(occ_file+strlen(file), ".occ");
   int fd_occ = open(occ_file, O_RDONLY);
   if (fd_occ == -1) {
      fprintf(stderr, "error opening '%s': %s\n", occ_file, strerror(errno));
      exit(EXIT_FAILURE);
   }
   free(occ_file);

   // Open SA file.
   char * sa_file = malloc(strlen(file)+4);
   strcpy(sa_file, file);
   strcpy(sa_file+strlen(file), ".sa");
   int fd_sa = open(sa_file, O_RDONLY);
   if (fd_sa == -1) {
      fprintf(stderr, "error opening '%s': %s\n", sa_file, strerror(errno));
      exit(EXIT_FAILURE);
   }
   free(sa_file);

   // Open GEN file.
   char * gen_file = malloc(strlen(file)+5);
   strcpy(gen_file, file);
   strcpy(gen_file+strlen(file), ".gen");
   int fd_gen = open(gen_file, O_RDONLY);
   if (fd_gen == -1) {
      fprintf(stderr, "error opening '%s': %s\n", gen_file, strerror(errno));
      exit(EXIT_FAILURE);
   }
   free(gen_file);

   // Load OCC index.
   long idxsize = lseek(fd_occ, 0, SEEK_END);
   lseek(fd_occ, 0, SEEK_SET);
   index->occ_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_occ, 0);
   close(fd_occ);
   if (index->occ_file == NULL) {
      fprintf(stderr, "error mmaping .occ index file: %s.\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   // Load SA index.
   idxsize = lseek(fd_sa, 0, SEEK_END);
   lseek(fd_sa, 0, SEEK_SET);
   index->sa_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_sa, 0);
   close(fd_sa);
   if (index->sa_file == NULL) {
      fprintf(stderr, "error mmaping .sa index file: %s.\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   // Load GEN index.
   idxsize = lseek(fd_gen, 0, SEEK_END);
   lseek(fd_gen, 0, SEEK_SET);
   index->gen_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_gen, 0);
   close(fd_gen);
   if (index->occ_file == NULL) {
      fprintf(stderr, "error mmaping .gen index file: %s.\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Set C and gsize.
   index->c = ((uint64_t *) index->occ_file);
   index->gsize = index->c[NUM_BASES];
   // Load occ tables.
   uint64_t occ_size = index->c[NUM_BASES+1];
   index->occ   = malloc(NUM_BASES * sizeof(void *));
   for (int i = 0; i < NUM_BASES; i++) {
      index->occ[i] = ((uint64_t *) index->occ_file) + NUM_BASES + 2 + i*occ_size;
   }
   // Load SA.
   index->pos = (uint64_t *) index->sa_file;
   // Load genome.
   index->genome = (char *) index->gen_file;

   return index;
}

ssize_t
write_index
(
 char * filename
)
{
   uint64_t gsize;

   // Allocate structures.
   uint64_t * C, * pos;
   uint64_t ** occ = malloc(NUM_BASES * sizeof(vstack_t *));

   // Parse genome and convert to integers. (0..NBASES-1)
   char * genome = compact_genome(filename, &gsize);
   C = compute_c(genome, gsize);

   // Output files.
   char * safile = malloc(strlen(filename)+4);
   strcpy(safile, filename);
   strcpy(safile + strlen(filename), ".sa");
   char * occfile = malloc(strlen(filename)+5);
   strcpy(occfile, filename);
   strcpy(occfile + strlen(filename), ".occ");
   char * genfile = malloc(strlen(filename)+5);
   strcpy(genfile, filename);
   strcpy(genfile + strlen(filename), ".gen");
   // Open files.
   int fsa = open(safile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int foc = open(occfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fgn = open(genfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   // Error control.
   if (fsa == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", safile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   if (foc == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", occfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   if (fgn == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", genfile, strerror(errno));
      exit(EXIT_FAILURE);
   }

   free(safile);
   free(occfile);
   free(genfile);

   // FM index
   uint64_t occ_size;
   bwt_index(genome, gsize, &pos, occ, &occ_size);

   // Write index.
   size_t s = 0, stot = 0;
   
   // .OCC FILE
   // Write C
   while (s < NUM_BASES*sizeof(uint64_t)) s += write(foc, C + (NUM_SYMBOLS-NUM_BASES) + s/sizeof(uint64_t), NUM_BASES*sizeof(uint64_t) - s);
   stot += s;
   // Write gsize (last value of C).
   s = 0;
   while (s < sizeof(uint64_t)) s += write(foc, &gsize, sizeof(uint64_t));
   stot += s;
   // Write OCC size.
   s = 0;
   while (s < sizeof(uint64_t)) s += write(foc, &occ_size, sizeof(uint64_t));
   stot += s;
   // Write OCC.
   for (int i = 0; i < NUM_BASES; i++) {
      s = 0;
      while (s < occ_size * sizeof(uint64_t)) s += write(foc,occ[i] + s/sizeof(uint64_t), occ_size*sizeof(uint64_t) - s);
      stot += s;
   }

   // .GEN FILE
   // Write genome bases.
   s = 0;
   while (s < gsize*sizeof(char)) s += write(fgn, genome + s/sizeof(char), gsize*sizeof(char) - s);
   stot += s;

   // .SA FILE
   // Write Suffix Array.
   s = 0;
   while (s < gsize*sizeof(uint64_t)) s += write(fsa, pos + s/sizeof(uint64_t), gsize*sizeof(uint64_t) - s);
   stot += s;

   return stot;
}


char *
compact_genome
(
 char     * filename,
 uint64_t * genomesize
)
{
   char * fileindex = malloc(strlen(filename)+5);
   strcpy(fileindex, filename);
   strcpy(fileindex + strlen(filename), ".chr");

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
   uint64_t gbufsize = GENOME_SIZE;
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
   uint64_t end = *genomesize - 1;
   char tmp;
   for (int i = 0 ; i < *genomesize/2; i++) {
      tmp = genome[i];
      genome[i] = genome[end-i];
      genome[end-i] = tmp;
   }
   genome[*genomesize] = 0;

   fclose(input);
   fclose(output);

   // Free memory.
   free(buffer);
   free(fileindex);
   
   return genome;
}


void
bwt_index
(
 char      * genome,
 uint64_t    gsize,
 uint64_t ** pos,
 uint64_t ** occ,
 uint64_t  * occ_size
)
{
   // Sort prefixes
   uint64_t * values = malloc((gsize+3)*sizeof(long));
   for (uint64_t i = 0; i < gsize; i++) values[i] = (uint64_t)translate[(int)genome[i]];
   values[gsize] = values[gsize+1] = values[gsize+2] = 0;
   uint64_t * sa = malloc(gsize*sizeof(uint64_t));
   suffixArray((long *)values, (long *)sa, gsize, NUM_BASES);
   free(values);
   *pos = sa;

   // Compute OCC. (MSB FIRST encoding)
   uint64_t occ_abs[NUM_BASES];
   uint64_t occ_intervals = (gsize+OCC_MARK_BITS-1)/OCC_MARK_BITS;
   uint64_t occ_words = occ_intervals * OCC_MARK_INTERVAL;
   uint64_t occ_marks = occ_intervals + 1;
   // Initial values.
   for (int i = 0; i < NUM_BASES; i++) {
      occ_abs[i]  = 0;
      occ[i] = calloc((occ_words+occ_marks),sizeof(uint64_t));
   }
   // Compute occ for the full genome.
   // Word starts at 1, this is the first mark.
   uint64_t word = 1;
   for (uint64_t i = 0, interval = 0; i < gsize; i++) {
      // Set occ bit.
      int base = translate[(uint8_t)genome[(sa[i] == 0 ? gsize-1 : sa[i]-1)]];
      if (base) {
         occ[base-1][word] |= 1;
         occ_abs[base-1]++;
      }
      // Next word.
      if (i % OCC_WORD_SIZE == OCC_WORD_SIZE-1) {
         word++; interval++;
         // Write Mark.
         if (interval == OCC_MARK_INTERVAL) {
            for (int j = 0; j < NUM_BASES; j++) occ[j][word] = occ_abs[j];
            word++;
            interval = 0;
         }
      }
      // Shift words.
      for (int j = 0; j < NUM_BASES; j++) occ[j][word] <<= 1;
   }
   // Shift last word.
   if (gsize%OCC_WORD_SIZE) {
      for (int j = 0; j < NUM_BASES; j++)
         occ[j][word] <<= OCC_WORD_SIZE - 1 - gsize%OCC_WORD_SIZE;
   }
   // Add last interval.
   for (int j = 0; j < NUM_BASES; j++) occ[j][occ_words+occ_marks-1] = occ_abs[j];
   // Set occ size.
   *occ_size = occ_words + occ_marks;
}


uint64_t *
compute_c
(
 char * genome,
 long gsize
)
{
   uint64_t * cnt = calloc(NUM_SYMBOLS, sizeof(long));
   for (uint64_t i = 0; i < gsize; i++) cnt[(int)translate[(int)genome[i]]]++;
   uint64_t * c = malloc(NUM_SYMBOLS * sizeof(long));
   c[0] = 0;
   for (int i = 1; i < NUM_SYMBOLS; i++) c[i] = c[i-1] + cnt[i-1];
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
   const int   reverse
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
      format = FASTA;
      break;
   case '@':
      format = FASTQ;
      break;
   default:
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
   char *q    = NULL;
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
            q = temp;
         } else continue;
         break;
      }

      /*size_t seqlen = strlen(seq);
      if (seqlen > MAXSEQLEN) {
         fprintf(stderr, "max sequence length exceeded (%d)\n", MAXSEQLEN);
         fprintf(stderr, "offending sequence:\n%s\n", seq);
         continue;
         }*/
      
      
      if (seq_push(&seqstack, tag, seq, q, reverse)) return NULL;
   }

   free(line);
   free(subline);
   free(temp);
   return seqstack;
}
