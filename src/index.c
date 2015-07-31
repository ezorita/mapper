#include "index.h"
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


int main(int argc, char *argv[])
{
   //   signal(SIGSEGV, SIGSEGV_handler);
   if (argc != 2) {
      fprintf(stderr, "usage: buildindex <genome file>\n");
      exit(EXIT_FAILURE);
   }

   return write_index(argv[1]);
}


int
write_index
(
 char * filename
)
{
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
   char * lutfile = malloc(strlen(filename)+5);
   strcpy(lutfile, filename);
   strcpy(lutfile + strlen(filename), ".lut");
   char * lcpfile = malloc(strlen(filename)+5);
   strcpy(lcpfile, filename);
   strcpy(lcpfile + strlen(filename), ".lcp");

   // Open files.
   int fsa = open(sarfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int foc = open(occfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int fgn = open(genfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int flt = open(lutfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   int flc = open(lcpfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
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
   if (flt == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", lutfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   if (flc == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", lcpfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   
   free(sarfile);
   free(occfile);
   free(genfile);
   free(lutfile);
   free(lcpfile);

   clock_t tstart;
   // Parse genome and reverse complement.
   uint64_t gsize;
   fprintf(stderr, "reading      genome file  ..."); tstart = clock();
   char * genome = compact_genome(filename, &gsize);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute suffix array.
   fprintf(stderr, "computing    suffix array ..."); tstart = clock();
   uint64_t * sa = (uint64_t *)compute_sa(genome, 2*gsize+1);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute OCC.
   uint64_t occ_size, wilcard_bwt;
   fprintf(stderr, "computing    occ table    ..."); tstart = clock();
   uint64_t * occ = compute_occ(genome, sa, 2*gsize+1, &occ_size, &wilcard_bwt);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compress suffix array.
   fprintf(stderr, "compressing  suffix array ..."); tstart = clock();
   uint64_t sa_bits = 0;
   while (2*gsize > (1 << sa_bits)) sa_bits++;
   uint64_t sa_size = compact_array(sa, 2*gsize+1, sa_bits);
   sa = realloc(sa, sa_size*sizeof(uint64_t));
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute C.
   fprintf(stderr, "computing    C table      ..."); tstart = clock();
   uint64_t * c = compute_c(occ + occ_size - NUM_BASES);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute LUT.
   fprintf(stderr, "computing    lookup table ..."); tstart = clock();
   uint64_t * lut = compute_lut(c, occ, LUT_KMER_SIZE);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compress LUT.
   fprintf(stderr, "compressing  lookup table ..."); tstart = clock();
   uint64_t lut_kmers = 1 << (2*LUT_KMER_SIZE);
   uint64_t lut_size = compact_array(lut, lut_kmers, sa_bits);
   lut = realloc(lut, lut_size*sizeof(uint64_t));
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute LCP.
   //   index_t index = {.size = 2*gsize+1, .c = c, .occ = occ, .wilcard_bwt = wilcard_bwt};
   uint64_t * lcp_index, lcp_index_size;
   fprintf(stderr, "computing    LCP intervals..."); tstart = clock();
   lcpstack_t * lcp = compute_lcp(&lcp_index_size, &lcp_index, LCP_MIN_DEPTH, sa_bits, sa, 2*gsize+1, genome);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   
   // Write index.
   size_t s = 0, stot = 0, bytes;
   // .OCC FILE
   fprintf(stderr, "writing occ...");
   // Write C
   while (s < (NUM_BASES+1)*sizeof(uint64_t)) s += write(foc, c + s/sizeof(uint64_t), (NUM_BASES+1)*sizeof(uint64_t) - s);
   stot += s;
   // Write OCC.
   s = 0;
   while (s < occ_size * sizeof(uint64_t)) s += write(foc,occ + s/sizeof(uint64_t), occ_size*sizeof(uint64_t) - s);
   stot += s;
   fprintf(stderr, " %ld bytes written.\n",stot);

   // .GEN FILE
   fprintf(stderr, "writing gen...");
   // Write genome bases (forward and reverse strand).
   s = 0;
   while (s < 2*gsize*sizeof(char)) s += write(fgn, genome + s/sizeof(char), 2*gsize*sizeof(char) - s);
   stot += s;
   fprintf(stderr, " %ld bytes written.\n",s);

   // .SAR FILE
   fprintf(stderr, "writing sar...");
   // Write sa word width.
   bytes = s = 0;
   while (s < sizeof(uint64_t)) s += write(fsa, &sa_bits, sizeof(uint64_t));
   bytes += s;
   stot += s;
   s = 0;
   while (s < sa_size*sizeof(uint64_t)) s += write(fsa, sa + s/sizeof(uint64_t), sa_size*sizeof(uint64_t) - s);
   stot += s;
   bytes += s;
   fprintf(stderr, " %ld bytes written.\n",bytes);


   // .LUT FILE
   fprintf(stderr, "writing lut...");
   bytes = s = 0;
   while (s < lut_size*sizeof(uint64_t)) s += write(flt, lut + s/sizeof(uint64_t), lut_size*sizeof(uint64_t) - s);
   stot += s;
   bytes += s;
   fprintf(stderr, " %ld bytes written.\n",bytes);

   // .LCP FILE
   fprintf(stderr, "writing lcp...");
   // LCP index.
   bytes = s = 0;
   while(s < lcp_index_size*sizeof(uint64_t)) s += write(flc, lcp_index + s/sizeof(uint64_t), lcp_index_size*sizeof(uint64_t) - s);
   stot += s;
   bytes += s;
   // LCP table.
   s = 0;
   while(s < lcp->pos) s += write(flc, lcp->lcp + s, lcp->pos - s);
   stot += s;
   bytes += s;
   fprintf(stderr, " %ld bytes written.\n",bytes);

   fprintf(stderr, "done. %ld bytes written.\n", stot);

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
 uint64_t * occ_size,
 uint64_t * wilcard_bwt
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
   // Word starts at NUM_BASES, after the first mark.
   uint64_t word = NUM_BASES;
   uint64_t interval = 0;
   for (uint64_t i = 0; i < gsize; i++) {
      // Set occ bit.
      int base = translate[(uint8_t)genome[(sa[i] == 0 ? gsize-1 : sa[i]-1)]];
      if (base == 5) *wilcard_bwt = i;
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

uint64_t *
compute_lut
(
 uint64_t * c,
 uint64_t * occ,
 int        depth
)
{
   // Alloc LUT table.
   uint64_t * lut = malloc((1 << (2*depth))*sizeof(uint64_t));
   if (lut == NULL) return NULL;
   // Recursively compute k-mer start-end indices.
   int path[depth];
   recursive_lut(c,occ,lut,0,0,depth,path);

   return lut;
}

void
recursive_lut
(
 uint64_t * c,
 uint64_t * occ,
 uint64_t * lut,
 uint64_t   ptr,
 int        d,
 int        maxd,
 int      * path
)
{
   if (d == maxd) {
      // Compute position from path.
      uint64_t idx = 0;
      for (int i = 0; i < maxd; i++) {
         idx += path[i] * (1 << (2*i));
      }
      lut[idx] = ptr;
      return;
   }
   // 'N' is not included.
   for (int i = 0; i < NUM_BASES - 1; i++) {
      uint64_t newptr = c[i] + get_occ_nt(ptr-1,occ,i);
      path[d] = i;
      recursive_lut(c,occ,lut,newptr,d+1,maxd,path);
   }
}

int
lcp_stack_push8
(
 lcpstack_t ** lcp_stack,
 uint8_t       lcp,
 int8_t        offset
)
{
   lcpstack_t * stack = *lcp_stack;
   // Realloc stack if full.
   if (stack->pos + 1 >= stack->size) {
      uint64_t newsize = stack->size * 2;
      *lcp_stack = stack = realloc(stack, sizeof(lcpstack_t) + newsize);
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   // Save LCP value and offset to parent.
   stack->lcp[stack->pos++] = lcp;
   stack->lcp[stack->pos++] = offset;
   return 0;
}

int
lcp_stack_push32
(
 lcpstack_t ** lcp_stack,
 uint8_t      lcp,
 int32_t      offset
)
{
   lcpstack_t * stack = *lcp_stack;
   // Realloc stack if full.
   if (stack->pos + 4 >= stack->size) {
      uint64_t newsize = (stack->size*2 <= stack->pos+4 ? stack->pos+5 : stack->size * 2);
      *lcp_stack = stack = realloc(stack, sizeof(lcpstack_t) + newsize);
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   // Save LCP value.
   stack->lcp[stack->pos++] = lcp;
   // Save extended offset to parent as int32.
   stack->lcp[stack->pos++] = (offset >> 24) & 0x000000FF;
   stack->lcp[stack->pos++] = (offset >> 16) & 0x000000FF;
   stack->lcp[stack->pos++] = (offset >> 8)  & 0x000000FF;
   stack->lcp[stack->pos++] = offset & 0x000000FF;
   return 0;
}

void
lcp_index_sample
(
 uint64_t * lcp_index,
 uint64_t   pos,
 int        size32
)
{
   // Count number of marks.
   uint64_t marks = pos / LCP_MARK_BITS + 1;
   uint64_t word  = pos / LCP_WORD_SIZE;
   int      bit   = pos % LCP_WORD_SIZE;
   // LSB index 0
   // ...
   // MSB index 63
   // Set bit.
   if ((lcp_index[2*word + marks] >> bit) & 1)
      fprintf(stderr, "repeated LCP!\n");
   lcp_index[2*word + marks] |= ((uint64_t)(1) << bit);
   if (size32) lcp_index[2*word + marks + 1] |= ((uint64_t)(1) << bit);
}

int
seq_lcp
(
 char * seq_a,
 char * seq_b
)
{
   int lcp = 0;
   while (seq_a[lcp] == seq_b[lcp]) lcp++;
   return lcp;
}

lcpcorner_t
corner_pop
(
 cstack_t * stack
)
{
   if (stack->pos > 0)
      return stack->c[--stack->pos];
   else
      return (lcpcorner_t){-1,-1};
}

int
corner_push
(
 cstack_t    ** stackp,
 lcpcorner_t    c
)
{
   cstack_t * stack = *stackp;
   // Realloc if needed.
   if (stack->pos >= stack->size) {
      uint64_t newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(cstack_t) + newsize*sizeof(lcpcorner_t));
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   // Save corner.
   stack->c[stack->pos++] = c;

   return 0;
}

cstack_t *
cstack_new
(
 int size
)
{
   cstack_t * stack = malloc(sizeof(cstack_t) + size*sizeof(lcpcorner_t));
   if (stack == NULL) return NULL;
   stack->size = size;
   stack->pos  = 0;
   return stack;
}

int64_t
naive_lcp
(
 int           min_depth,
 uint64_t      idxsize,
 uint64_t    * lcpindex,
 lcpstack_t ** lcp,
 uint64_t    * sar,
 int           sar_bits,
 char        * genome
)
{
   cstack_t * top = cstack_new(STACK_LCP_INITIAL_SIZE);
   cstack_t * bottom = cstack_new(STACK_LCP_INITIAL_SIZE);
   // Add top corner.
   corner_push(&top, (lcpcorner_t){0,0});
   // Previous position and lcp (0).
   uint64_t pp = get_sa(0,sar,sar_bits);
   int pl = 0;
   for (uint64_t i = 1; i < idxsize; i++) {
      // Current position and lcp.
      uint64_t cp = get_sa(i,sar,sar_bits);
      int cl = seq_lcp(genome+pp, genome+cp);
      // If current LCP > previous LCP -> Top corner.
      if (cl > pl) {
         // Save sample on top of the top corner stack.
         lcpcorner_t tcorner = top->c[top->pos-1];
         int offset = tcorner.pos - (i-1);
         if (offset < -127) {
            lcp_index_sample(lcpindex, i-1, 1);
            if (lcp_stack_push32(lcp, pl, offset)) return -1;
         } else {
            lcp_index_sample(lcpindex, i-1, 0);
            if (lcp_stack_push8(lcp, pl, (char)offset)) return -1;
         }
         // Push corner to the stack.
         corner_push(&top, (lcpcorner_t){pl, i-1});
      } else if (cl < pl) {
         while (bottom->pos > 0) {
            if (bottom->c[bottom->pos-1].lcp <= cl) break;
            // Pop sample from bottom stack and save.
            lcpcorner_t bcorner = corner_pop(bottom);
            int offset = i-1 - bcorner.pos;
            if (offset > 127) {
               lcp_index_sample(lcpindex, i-1, 1);
               if (lcp_stack_push32(lcp, bcorner.lcp, offset)) return -1;
            } else {
               lcp_index_sample(lcpindex, i-1, 0);
               if (lcp_stack_push8(lcp, bcorner.lcp, (char)offset)) return -1;
            }
         }
         // Save current bottom corner.
         corner_push(&bottom, (lcpcorner_t){cl, i-1});
         // Remove top corners from stack.
         while (top->pos > 0) {
            if (top->c[top->pos-1].lcp < cl) break;
            top->pos--;
         }
      }
      // Assign previous.
      pp = cp;
      pl = cl;
   }
   // Save last sample.
   lcp_index_sample(lcpindex, idxsize-1, 0);
   lcp_stack_push32(lcp, pl, 0);
   
   return 0;

}


lcpstack_t *
compute_lcp
(
 uint64_t  * index_size,
 uint64_t ** lcp_index,
 int         min_depth,
 int         sar_bits,
 uint64_t  * sar,
 uint64_t    idxsize,
 char      * genome
)
{
   // Allocate LCP stack.
   lcpstack_t * stack = malloc(sizeof(lcpstack_t) + STACK_LCP_INITIAL_SIZE*sizeof(uint8_t));
   if (stack == NULL) return NULL;
   stack->pos = 0;
   stack->size = STACK_LCP_INITIAL_SIZE;
   // Allocate LCP index.
   uint64_t words = (idxsize + LCP_WORD_SIZE - 1)/LCP_WORD_SIZE;
   uint64_t intervals = (words + LCP_MARK_INTERVAL - 1) / LCP_MARK_INTERVAL;
   uint64_t marks = intervals + 1;
   uint64_t * lcpindex = *lcp_index = calloc(2*intervals*LCP_MARK_INTERVAL + marks,sizeof(uint64_t));
   *index_size = 2*intervals*LCP_MARK_INTERVAL + marks;
   /* Recursive
   // Run recursive search.
   fmdpos_t root = {.fp = 0, .rp = 0, .sz = index->size};
   // Insert sample at index 0.
   lcp_index_sample(lcpindex, 0, 0);
   lcp_stack_push8 (&stack, 0, (int8_t) 0);

   int64_t lcpval = recursive_lcp(root, 0, index, lcpindex, &stack, min_depth);
   */
   /* Naive */
   naive_lcp(min_depth, idxsize, lcpindex, &stack, sar, sar_bits, genome);
   /* Recursive
   // Insert sample at index size-1.
   lcp_index_sample(lcpindex, index->size - 1, 0);
   lcp_stack_push8 (&stack, lcpval, (int8_t) 0);
   */
   // Realloc stack.
   stack = realloc(stack, sizeof(lcpstack_t) + stack->pos);
   // LCP index marks.
   uint64_t abs_pos = 0;
   for (uint64_t i = 0, w = 0; i < words; i++) {
      if (i%LCP_MARK_INTERVAL == 0) lcpindex[w++] = abs_pos;
      abs_pos += 2*__builtin_popcountl(lcpindex[w++]);
      abs_pos += 3*__builtin_popcountl(lcpindex[w++]);
   }
   lcpindex[*index_size - 1] = abs_pos;
   // Return stack.
   return stack;
}

int64_t
recursive_lcp
(
 fmdpos_t      pos,
 int           depth,
 index_t     * index,
 uint64_t    * lcp_index,
 lcpstack_t ** lcp,
 int           mindepth
)
{
   if (pos.sz == 0) return -1;
   if (pos.sz == 1) return depth-1;
   // Maximum depth.
   if (depth == 255) return depth-1;
   int64_t lcpval = depth;
   for (int i = 0; i < NUM_BASES; i++) {
      // Forward search nucleotide 'i'.
      fmdpos_t newpos = extend_fw(i, pos, index);
      // No more Suffixes.
      if (newpos.sz == 0) continue;
      // One suffix exists.
      if (newpos.sz == 1) {
         lcpval = depth;
         continue;
      }
      // 2+ suffixes exist. (New interval)
      if ((newpos.fp > pos.fp) && (depth >= mindepth)) {
         //save sp at depth using offset sp - newsp.
         int64_t offset = (int64_t)pos.fp - (int64_t)newpos.fp;
         int ext_size = offset < -127;
         lcp_index_sample(lcp_index, newpos.fp, ext_size);
         if (ext_size) lcp_stack_push32(lcp, depth, (int32_t) offset);
         else          lcp_stack_push8 (lcp, depth, (int8_t)  offset);
      }
      // Extend suffix.
      lcpval = recursive_lcp(newpos, depth+1, index, lcp_index, lcp, mindepth);
      // Save ep if it does not coincide with last position.
      if ((newpos.fp + newpos.sz < pos.fp + pos.sz)) {
         if (lcpval >= mindepth) {
            int64_t offset = (int64_t)pos.fp + pos.sz - (int64_t)newpos.fp - newpos.sz;
            int ext_size = offset > 127;
            uint64_t p = newpos.fp + newpos.sz - 1;
            lcp_index_sample(lcp_index, p, ext_size);
            if (ext_size) lcp_stack_push32(lcp, lcpval, (int32_t) offset);
            else          lcp_stack_push8 (lcp, lcpval, (int8_t)  offset);
         }
      } else break;
   }
   return lcpval;
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

   // Insert wilcard at the end.
   genome[2*div] = '$';
   genome[2*div+1] = 0;

   fclose(input);
   fclose(output);

   // Free memory.
   free(buffer);
   free(fileindex);
   
   return genome;
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
         lastbit = 64 - lastbit;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         array[++word] = (current & mask) >> (bits - lastbit);
      }
   }

   return word + (lastbit > 0);
}
