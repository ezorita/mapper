#include "indexbuild.h"
#include <signal.h>
#include <execinfo.h>

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

int
write_index
(
 char * filename,
 int    def_kmer,
 int    def_tau,
 int    def_rthr,
 int    threads
)
{
   // Output files.
   char * sarfile = add_suffix(filename, ".sar");
   char * bwtfile = add_suffix(filename, ".bwt");
   char * genfile = add_suffix(filename, ".gen");
   char * annfile = add_suffix(filename, ".ann");
   char * shtfile = add_suffix(filename, ".sht");

   // Open files.
   int fsa = open(sarfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fbw = open(bwtfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fgn = open(genfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fan = open(annfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fsh = open(shtfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);

   // Error control.
   if (fsa == -1) {
      fprintf(stderr, "[error] opening '%s' file: %s.\n", sarfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   if (fbw == -1) {
      fprintf(stderr, "[error] opening '%s' file: %s.\n", bwtfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   if (fgn == -1) {
      fprintf(stderr, "[error] opening '%s' file: %s.\n", genfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   if (fan == -1) {
      fprintf(stderr, "[error] opening '%s' file: %s.\n", annfile, strerror(errno));
      return EXIT_FAILURE;
   }
   if (fsh == -1) {
      fprintf(stderr, "[error] opening '%s' file: %s.\n", shtfile, strerror(errno));
      return EXIT_FAILURE;
   }
   
   free(sarfile);
   free(bwtfile);
   free(genfile);
   free(annfile);
   free(shtfile);


   clock_t tstart;
   size_t s = 0, stot = 0, bytes;
   // Parse genome and reverse complement.
   uint64_t gsize;
   fprintf(stderr, "[proc] reading      genome file  ..."); tstart = clock();
   char * genome = compact_genome(filename, &gsize);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);

   // Compute suffix array.
   fprintf(stderr, "[proc] computing    suffix array ..."); tstart = clock();
   uint64_t * sa = (uint64_t *)compute_sa(genome, gsize);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);

   // Compute OCC.
   uint64_t occ_size;
   fprintf(stderr, "[proc] computing    occ table    ..."); tstart = clock();
   uint64_t * occ = compute_occ(genome, sa, gsize, &occ_size);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);

   // Compress suffix array.
   fprintf(stderr, "[proc] compressing  suffix array ..."); tstart = clock();
   uint64_t sa_bits = 0;
   while (gsize > ((uint64_t)1 << sa_bits)) sa_bits++;
   uint64_t sa_size = compact_array(sa, gsize, sa_bits);
   sa = realloc(sa, sa_size*sizeof(uint64_t));
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);

   // Compute C.
   fprintf(stderr, "[proc] computing    C table      ..."); tstart = clock();
   uint64_t * c = compute_c(occ + occ_size - NUM_BASES);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);

   // Write .BWT file
   fprintf(stderr, "[proc] writing bwt...");
   // mark interval
   s = 0;
   uint64_t mark_int = OCC_MARK_INTERVAL;
   while (s < sizeof(uint64_t)) s += write(fbw, &mark_int, sizeof(uint64_t));
   stot += s;
   // Write C
   s  = 0;
   while (s < (NUM_BASES+1)*sizeof(uint64_t)) s += write(fbw, c + s/sizeof(uint64_t), (NUM_BASES+1)*sizeof(uint64_t) - s);
   stot += s;
   // Write OCC.
   s = 0;
   while (s < occ_size * sizeof(uint64_t)) s += write(fbw,occ + s/sizeof(uint64_t), occ_size*sizeof(uint64_t) - s);
   stot += s;
   close(fbw);
   fprintf(stderr, " %ld bytes written.\n",stot);
   free(c); free(occ);
   // Write .GEN FILE
   fprintf(stderr, "[proc] writing gen...");
   // Write genome bases (forward and reverse strand).
   s = 0;
   while (s < gsize*sizeof(char)) s += write(fgn, genome + s/sizeof(char), gsize*sizeof(char) - s);
   stot += s;
   close(fgn);
   free(genome);
   fprintf(stderr, " %ld bytes written.\n",s);
   // .SAR FILE
   fprintf(stderr, "[proc] writing sar...");
   // Write sa word width.
   bytes = s = 0;
   while (s < sizeof(uint64_t)) s += write(fsa, &sa_bits, sizeof(uint64_t));
   bytes += s;
   s = 0;
   while (s < sa_size*sizeof(uint64_t)) s += write(fsa, sa + s/sizeof(uint64_t), sa_size*sizeof(uint64_t) - s);
   bytes += s;
   stot += bytes;
   close(fsa);
   free(sa);
   fprintf(stderr, " %ld bytes written.\n",bytes);
   fprintf(stderr, "[info] index base created successfully. %ld bytes written.\n", stot);
   fprintf(stderr, "[info] default genome annotation: k:%d, d:%d, repeat_thr:%d.\n", def_kmer, def_tau, def_rthr);
   fprintf(stderr, "[proc] reading index... ");
   index_t * index = index_load_base(filename);
   if (index == NULL)
      return EXIT_FAILURE;
   fprintf(stderr, "done.\n");
   fprintf(stderr, "[proc] creating annotation index... ");
   uint8_t count = 0;
   if (write(fan, &count, sizeof(uint8_t)) < 1) {
      fprintf(stderr, "\n[error] could not write annotation index.\n");
      return EXIT_FAILURE;
   }
   close(fan);
   fprintf(stderr, "done.\n[proc] creating seed table index... ");
   if(write(fsh, &count, sizeof(uint8_t)) < 1) {
      fprintf(stderr, "\n[error] could not write seed table index.\n");
      return EXIT_FAILURE;
   }
   close(fsh);
   fprintf(stderr, "done.\n");
   index_add_annotation(def_kmer,def_tau,def_rthr,3,threads,index,filename);
   return 0;
}

int64_t *
compute_sa
(
 char    * genome,
 uint64_t  gsize
)
{
   char * data = malloc(gsize);
   for (uint64_t i = 0; i < gsize; i++) data[i] = uppercase[(int)genome[i]];
   int64_t * sa = malloc(gsize*sizeof(int64_t));
   if (sa == NULL) return NULL;
   divsufsort((unsigned char *) data, sa, gsize);
   free(data);
   return sa;
}

uint64_t *
compute_occ
(
 char     * genome,
 uint64_t * sa,
 uint64_t   gsize,
 uint64_t * occ_size
)
{
   // Words, marks and intervals.
   uint64_t occ_intervals = (gsize+OCC_MARK_BITS-1)/OCC_MARK_BITS;
   uint64_t occ_words = occ_intervals * OCC_MARK_INTERVAL;
   uint64_t occ_marks = occ_intervals + 1;

   // Alloc variables.
   uint64_t * occ = malloc((occ_words+occ_marks)*NUM_BASES*sizeof(uint64_t));
   if (occ == NULL) return NULL;
   uint64_t occ_abs[NUM_BASES+1];
   uint64_t occ_tmp[NUM_BASES+1];

   // Initial values.
   for (int i = 0; i < NUM_BASES; i++) {
      occ_abs[i] = 0;
      occ_tmp[i] = 0;
      occ[i] = 0;
   }

   // Compute OCC. (MSB FIRST encoding)
   uint64_t word = NUM_BASES, interval = 0;

   for (uint64_t i = 0; i < gsize; i++) {
      // Set occ bit.
      int base = translate[(uint8_t)genome[(sa[i] > 0 ? sa[i] - 1 : gsize - 1)]];
      occ_tmp[base] |= 1;
      occ_abs[base]++;
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
   if (gsize%OCC_WORD_SIZE) {
      interval++;
      for (int j = 0; j < NUM_BASES; j++)
         occ[word++] = occ_tmp[j] << (OCC_WORD_SIZE - 1 - gsize%OCC_WORD_SIZE);
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
   *occ_size = word;

   return occ;
}

uint64_t *
compute_c
(
 uint64_t * occ
)
{
   // Fill C.
   uint64_t * C = malloc((NUM_BASES+1)*sizeof(uint64_t));
   if (C == NULL) return NULL;
   // Count the wildcard before 'A'.
   C[0] = 1;
   for (int i = 1; i <= NUM_BASES; i++) C[i] = C[i-1] + occ[i-1];
   return C;
}


uint64_t
compact_array
(
 uint64_t * array,
 uint64_t     len,
 int         bits
)
{
   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64-bits);
   uint64_t word = 0;
   int   lastbit = 0;
   
   // Clear upper bits of array[0].
   array[0] &= mask;

   for (uint64_t i = 0, current; i < len; i++) {
      // Save the current value.
      current = array[i];
      // Store the compact version.
      array[word] |= (current & mask) << lastbit;
      // Update bit offset.
      lastbit += bits;
      // Word is full.
      if (lastbit >= 64) {
         lastbit = lastbit - 64;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         array[++word] = (current & mask) >> (bits - lastbit);
      }
   }

   return word + (lastbit > 0);
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
   genome = realloc(genome, 2*(*genomesize)+2);
   if (genome == NULL) {
      fprintf(stderr, "error in 'compact_genome' (realloc): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // Reverse complement.
   uint64_t div = *genomesize;
   for (int i = 0 ; i < div; i++)
      genome[div + i] = revcomp[(int)genome[div-i-1]];

   *genomesize = *genomesize * 2 + 1;

   // Insert wildcard at the end.
   genome[2*div] = '$';
   genome[2*div+1] = 0;

   fclose(input);
   fclose(output);

   // Free memory.
   free(buffer);
   free(fileindex);
   
   return genome;
}
