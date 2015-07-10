#include "bwmapper.h"

static const char sa_value[256] = {[0 ... 255] = 5,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['t'] = 4, ['n'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4, ['N'] = 5 };

static const char translate[256] = {[0 ... 255] = 4,
                           ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['n'] = 4,
                           ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['N'] = 4 };

static const char revert[256]  = {[0 ... 255] = 'N', [0] = 'A', [1] = 'C', [2] = 'G', [3] = 'T'};

static const char revcomp[256] = {[0 ... 255] = 'N',
                   ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A' };

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
   char * sa_file = malloc(strlen(file)+5);
   strcpy(sa_file, file);
   strcpy(sa_file+strlen(file), ".sar");
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
      fprintf(stderr, "error mmaping .sar index file: %s.\n", strerror(errno));
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
   index->c = (uint64_t *) index->occ_file;
   // Genome size (Forward strand only).
   index->gsize = index->c[NUM_BASES]/2;
   // Load occ tables.
   index->occ = index->c + NUM_BASES + 1;
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
   uint64_t * occ;

   // Parse genome and convert to integers. (0..NBASES-1)
   char * genome = compact_genome(filename, &gsize);

   // Output files.
   char * sarfile = malloc(strlen(filename)+5);
   strcpy(sarfile, filename);
   strcpy(sarfile + strlen(filename), ".sar");
   char * occfile = malloc(strlen(filename)+5);
   strcpy(occfile, filename);
   strcpy(occfile + strlen(filename), ".occ");
   char * genfile = malloc(strlen(filename)+5);
   strcpy(genfile, filename);
   strcpy(genfile + strlen(filename), ".gen");
   // Open files.
   int fsa = open(sarfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int foc = open(occfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fgn = open(genfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   // Error control.
   if (fsa == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", sarfile, strerror(errno));
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

   free(sarfile);
   free(occfile);
   free(genfile);

   // FM index (the size of the index is 2*gsize)
   uint64_t occ_size;
   bwt_index(genome, 2*gsize, &pos, &occ, &occ_size, &C);

   // Write index.
   size_t s = 0, stot = 0;
   
   // .OCC FILE
   // Write C
   while (s < (NUM_BASES+1)*sizeof(uint64_t)) s += write(foc, C + s/sizeof(uint64_t), (NUM_BASES+1)*sizeof(uint64_t) - s);
   stot += s;
   // Write OCC.
   s = 0;
   while (s < occ_size * sizeof(uint64_t)) s += write(foc,occ + s/sizeof(uint64_t), occ_size*sizeof(uint64_t) - s);
   stot += s;

   // .GEN FILE
   // Write genome bases (forward and reverse strand).
   s = 0;
   while (s < 2*gsize*sizeof(char)) s += write(fgn, genome + s/sizeof(char), 2*gsize*sizeof(char) - s);
   stot += s;

   // .SAR FILE
   // Write Suffix Array.
   int bits = 0;
   while (2*gsize < (1 << bits)) bits++;
   uint64_t sa_size = compact_array(pos, 2*gsize, bits);
   s = 0;
   while (s < sa_size*sizeof(uint64_t)) s += write(fsa, pos + s/sizeof(uint64_t), sa_size*sizeof(uint64_t) - s);
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

   // Initialize size
   *genomesize = 0;

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
   genome = realloc(genome, 2*(*genomesize)+1);
   if (genome == NULL) {
      fprintf(stderr, "error in 'compact_genome' (realloc): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Reverse complement.
   uint64_t div = *genomesize;
   for (int i = 0 ; i < div; i++)
      genome[div + i] = revcomp[(int)genome[div-i-1]];

   // Insert wilcard at the end.
   genome[2*div] = 0;

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
 uint64_t    idxsize,
 uint64_t ** pos,
 uint64_t ** occp,
 uint64_t  * occ_size,
 uint64_t ** Cp
)
{
   // Sort prefixes
   uint64_t * values = malloc((idxsize+3)*sizeof(long));
   for (uint64_t i = 0; i < idxsize; i++) values[i] = (uint64_t)sa_value[(int)genome[i]];
   values[idxsize] = values[idxsize+1] = values[idxsize+2] = 0;
   uint64_t * sa = malloc(idxsize*sizeof(uint64_t));
   suffixArray((long *)values, (long *)sa, idxsize, NUM_BASES);
   free(values);
   *pos = sa;

   // Words, marks and intervals.
   uint64_t occ_intervals = (idxsize+OCC_MARK_BITS-1)/OCC_MARK_BITS;
   uint64_t occ_words = occ_intervals * OCC_MARK_INTERVAL;
   uint64_t occ_marks = occ_intervals + 1;

   // Alloc variables.
   uint64_t * occ = *occp = malloc((occ_words+occ_marks)*NUM_BASES*sizeof(uint64_t));
   uint64_t occ_abs[NUM_BASES];
   uint64_t occ_tmp[NUM_BASES];

   // Initial values.
   for (int i = 0; i < NUM_BASES; i++) {
      occ_abs[i] = 0;
      occ_tmp[i] = 0;
      occ[i] = 0;
   }

   // Compute OCC. (MSB FIRST encoding)
   // Word starts at NUM_BASES, after the first mark.
   uint64_t word = NUM_BASES;
   uint64_t interval = 0;
   for (uint64_t i = 0; i < idxsize; i++) {
      // Set occ bit.
      int base = translate[(uint8_t)genome[(sa[i] == 0 ? idxsize-1 : sa[i]-1)]];
      occ_tmp[base] |= 1;
      occ_abs[base-1]++;
      // Next word.
      if (i % OCC_WORD_SIZE == OCC_WORD_SIZE-1) {
         for (int j = 0; j < NUM_BASES; j++) {
            occ[word++] = occ_tmp[j];
            occ_tmp[j] = 0;
         }
         interval++;
         // Write Mark.
         if (interval == OCC_MARK_INTERVAL) {
            for (int j = 0; j < NUM_BASES; j++) occ[word++] = occ_abs[j];
            interval = 0;
         }
      }
      // Shift words.
      for (int j = 0; j < NUM_BASES; j++) occ_tmp[j] <<= 1;
   }
   // Shift last word.
   if (idxsize%OCC_WORD_SIZE) {
      interval++;
      for (int j = 0; j < NUM_BASES; j++)
         occ[word++] = occ_tmp[j] << (OCC_WORD_SIZE - 1 - idxsize%OCC_WORD_SIZE);
   }
   // Add last mark.
   if (interval > 0) {
      // Fill the last interval with 0.
      for (int i = interval; i < OCC_MARK_INTERVAL; i++)
         for (int j = 0; j < NUM_BASES; j++) occ[word++] = 0;
      // Add mark.
      for (int j = 0; j < NUM_BASES; j++) occ[word++] = occ_abs[j];
   }
   // Set occ size.
   //   *occ_size = occ_words + occ_marks;
   *occ_size = word;

   // Fill C.
   uint64_t * C = *Cp = malloc((NUM_BASES+1)*sizeof(uint64_t));
   C[0] = 0;
   for (int i = 1; i <= NUM_BASES; i++) C[i] = C[i-1] + occ_abs[i-1];
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
