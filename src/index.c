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
   //signal(SIGSEGV, SIGSEGV_handler);
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
   size_t s = 0, stot = 0, bytes;
   // Parse genome and reverse complement.
   uint64_t gsize;
   fprintf(stderr, "reading      genome file  ..."); tstart = clock();
   char * genome = compact_genome(filename, &gsize);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute suffix array.
   fprintf(stderr, "computing    suffix array ..."); tstart = clock();
   uint64_t * sa = (uint64_t *)compute_sa(genome, gsize);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compute OCC.
   uint64_t occ_size;
   fprintf(stderr, "computing    occ table    ..."); tstart = clock();
   uint64_t * occ = compute_occ(genome, sa, gsize, &occ_size);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Compress suffix array.
   fprintf(stderr, "compressing  suffix array ..."); tstart = clock();
   uint64_t sa_bits = 0;
   while (gsize > (1 << sa_bits)) sa_bits++;
   uint64_t sa_size = compact_array(sa, gsize, sa_bits);
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
   // Write .OCC file
   fprintf(stderr, "writing occ...");
   // mark interval
   s = 0;
   uint64_t mark_int = OCC_MARK_INTERVAL;
   while (s < sizeof(uint64_t)) s += write(foc, &mark_int, sizeof(uint64_t));
   stot += s;
   // Write C
   s  = 0;
   while (s < (NUM_BASES+1)*sizeof(uint64_t)) s += write(foc, c + s/sizeof(uint64_t), (NUM_BASES+1)*sizeof(uint64_t) - s);
   stot += s;
   // Write OCC.
   s = 0;
   while (s < occ_size * sizeof(uint64_t)) s += write(foc,occ + s/sizeof(uint64_t), occ_size*sizeof(uint64_t) - s);
   stot += s;
   fprintf(stderr, " %ld bytes written.\n",stot);
   free(c); free(occ);
   // Write .LUT file
   fprintf(stderr, "writing lut...");
   bytes = s = 0;
   uint64_t kmer_size = LUT_KMER_SIZE;
   while (s < sizeof(uint64_t)) s += write(flt, &kmer_size, sizeof(uint64_t));
   bytes += s;
   s = 0;
   while (s < lut_size*sizeof(uint64_t)) s += write(flt, lut + s/sizeof(uint64_t), lut_size*sizeof(uint64_t) - s);
   bytes += s;
   stot += bytes;
   fprintf(stderr, " %ld bytes written.\n",bytes);
   free(lut);
   // Compute LCP.
   fprintf(stderr, "computing    LCP intervals..."); tstart = clock();
   lcp_t lcp = compute_lcp(gsize, LCP_MIN_DEPTH, sa_bits, sa, genome);
   fprintf(stderr, " done [%.3fs]\n", (clock()-tstart)*1.0/CLOCKS_PER_SEC);
   // Write .GEN FILE
   fprintf(stderr, "writing gen...");
   // Write genome bases (forward and reverse strand).
   s = 0;
   while (s < gsize*sizeof(char)) s += write(fgn, genome + s/sizeof(char), gsize*sizeof(char) - s);
   stot += s;
   fprintf(stderr, " %ld bytes written.\n",s);
   // .SAR FILE
   fprintf(stderr, "writing sar...");
   // Write sa word width.
   bytes = s = 0;
   while (s < sizeof(uint64_t)) s += write(fsa, &sa_bits, sizeof(uint64_t));
   bytes += s;
   s = 0;
   while (s < sa_size*sizeof(uint64_t)) s += write(fsa, sa + s/sizeof(uint64_t), sa_size*sizeof(uint64_t) - s);
   bytes += s;
   stot += bytes;
   fprintf(stderr, " %ld bytes written.\n",bytes);
   // .LCP FILE
   fprintf(stderr, "writing lcp...");
   // LCP index.
   bytes = s = 0;
   mark_int = LCP_MARK_INTERVAL;
   while (s < sizeof(uint64_t)) s += write(flc, &mark_int, sizeof(uint64_t));
   bytes += s;
   s = 0;
   uint64_t min_depth = LCP_MIN_DEPTH;
   while (s < sizeof(uint64_t)) s += write(flc, &min_depth, sizeof(uint64_t));
   bytes += s;
   s = 0;
   while(s < sizeof(uint64_t)) s += write(flc, &(lcp.idx_size), sizeof(uint64_t));
   bytes += s;
   // Write sample index.
   s = 0;
   while(s < lcp.idx_size*sizeof(uint64_t)) s += write(flc, lcp.idx_sample + s/sizeof(uint64_t), lcp.idx_size*sizeof(uint64_t) - s);
   bytes += s;
   // Write extended index.
   s = 0;
   while(s < lcp.idx_size*sizeof(uint64_t)) s += write(flc, lcp.idx_extend + s/sizeof(uint64_t), lcp.idx_size*sizeof(uint64_t) - s);
   bytes += s;
   // Write LCP samples.
   s = 0;
   uint64_t nsamples = (lcp.lcp_sample)->pos;
   while(s < sizeof(uint64_t)) s+= write(flc, &nsamples, sizeof(uint64_t));
   bytes += s;
   s = 0;
   while(s < nsamples*sizeof(int8_t)) s += write(flc, (lcp.lcp_sample)->val + s/sizeof(int8_t), nsamples * sizeof(int8_t) - s);
   bytes += s;
   // Write ext samples.
   s = 0;
   nsamples = (lcp.lcp_extend)->pos;
   while(s < sizeof(uint64_t)) s+= write(flc, &nsamples, sizeof(uint64_t));
   bytes += s;
   s = 0;
   while(s < nsamples*sizeof(int32_t)) s += write(flc, (lcp.lcp_extend)->val + s/sizeof(int32_t), nsamples * sizeof(int32_t) - s);
   bytes += s;
   stot += bytes;
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
   for (uint64_t i = 0; i < gsize; i++) {
      // Set occ bit.
      int base = translate[(uint8_t)genome[(sa[i] == 0 ? gsize-1 : sa[i]-1)]];
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
   C[0] = 0;
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
stack_push_lcp
(
 stack8_t **  stackp,
 uint8_t      lcp,
 int32_t      offset
)
{
   stack8_t * stack = *stackp;
   // Realloc stack if full.
   if (stack->pos + 4 >= stack->size) {
      uint64_t newsize = (stack->size*2 <= stack->pos+4 ? stack->pos+5 : stack->size * 2);
      *stackp = stack = realloc(stack, sizeof(stack8_t) + newsize);
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   // Save LCP value.
   stack->val[stack->pos++] = lcp;
   // Save extended offset to parent as int32.
   memcpy(stack->val + stack->pos, &offset, sizeof(int32_t));
   stack->pos += 4;

   return 0;
}

void
lcp_index_set
(
 uint64_t * lcp_index,
 uint64_t   pos
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
   lcp_index[word + marks] |= ((uint64_t)(1) << bit);
}

int
seq_lcp
(
 char * seq_a,
 char * seq_b
)
{
   int lcp = 0;
   while (seq_a[lcp] == seq_b[lcp] && lcp < 255) lcp++;
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
      return (lcpcorner_t){-1,-1,-1};
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

int
sample_compar
(
 const void * a,
 const void * b,
 const int c
)
{
   lcpcorner_t * ca = (lcpcorner_t *)a;
   lcpcorner_t * cb = (lcpcorner_t *)b;
   return (ca->pos > cb->pos ? 1 : -1);
}

stack8_t *
stack8_new
(
 int size
)
{
   stack8_t * stack = malloc(sizeof(stack8_t) + size);
   if (stack == NULL) return NULL;
   stack->pos = 0;
   stack->size = size;
   return stack;
}

stack32_t *
stack32_new
(
 int size
)
{
   stack32_t * stack = malloc(sizeof(stack32_t) + size * sizeof(int));
   if (stack == NULL) return NULL;
   stack->pos = 0;
   stack->size = size;
   return stack;
}

int32_t
stack32_push
(
 stack32_t ** stackp,
 int32_t      val
)
{
   stack32_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      uint64_t newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(stack32_t) + newsize * sizeof(int));
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   stack->val[stack->pos++] = val;
   return 0;
}

int32_t
stack32_pop
(
 stack32_t * stack
)
{
   return stack->val[--stack->pos];
}


int64_t
naive_lcp
(
 int           min_depth,
 uint64_t      idxsize,
 uint64_t    * lcpsample,
 uint64_t    * lcpext,
 stack8_t   ** lcp,
 stack32_t  ** ext,
 uint64_t    * sar,
 int           sar_bits,
 char        * genome
)
{
   cstack_t * top = cstack_new(STACK_LCP_INITIAL_SIZE);
   cstack_t * bottom = cstack_new(STACK_LCP_INITIAL_SIZE);
   // Push topmost corner.
   corner_push(&top, (lcpcorner_t){-1, 0, 0});
   if (min_depth == 0) {
      lcp_index_set(lcpsample, 0);
      if (stack_push_lcp(lcp, 0, -1)) return -1;
   }
   // Compute LCP of i=1, store in previous.
   uint64_t cp = get_sa(0,sar,sar_bits);
   int cl = 0;
   uint64_t pp = get_sa(1,sar,sar_bits);
   int pl = seq_lcp(genome+pp, genome+cp);
   for (uint64_t i = 2; i < idxsize; i++) {
      // Current position and lcp.
      cp = get_sa(i,sar,sar_bits);
      cl = seq_lcp(genome+pp, genome+cp);
      // If current LCP > previous LCP -> Previous was top corner.
      if (cl > pl) {
         // Top corner - save sample.
         lcpcorner_t tcorner = top->c[top->pos-1];
         int offset = tcorner.pos - (i-1);
         if (cl >= min_depth) {
            lcp_index_set(lcpsample, i-1);
            if (offset < -127) lcp_index_set(lcpext, i-1);
            if (stack_push_lcp(lcp, pl, offset)) return -1;
         }
         // Push corner to the stack.
         corner_push(&top, (lcpcorner_t){pl, cl, i-1, 0});
      }
      // If current LCP < previous LCP -> Previous was bottom corner.
      else if (cl < pl) {
         // Bottom corner - save sample.
         if (pl >= min_depth) {
            if (stack_push_lcp(lcp, pl, 0xFFFFFFFF)) return -1;
         }
         while (bottom->pos > 0) {
            if (bottom->c[bottom->pos-1].lcp_next <= cl) break;
            // Pop sample from bottom stack and save.
            lcpcorner_t bcorner = corner_pop(bottom);
            // Bottom corner - update sample.
            if (bcorner.lcp >= min_depth) {
               int offset = i-1 - bcorner.pos;
               memcpy((*lcp)->val + bcorner.ptr, &offset, sizeof(int));
               lcp_index_set(lcpsample, bcorner.pos);
               if (offset > 127) lcp_index_set(lcpext, bcorner.pos);
            }
         }
         // Save current bottom corner.
         corner_push(&bottom, (lcpcorner_t){pl, cl, i-1, (*lcp)->pos - 4});
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

   // Add remaining bottom nodes.
   while (bottom->pos > 0) {
      // Pop sample from bottom stack and save.
      lcpcorner_t bcorner = corner_pop(bottom);
      // Bottom corner - update sample.
      if (bcorner.lcp >= min_depth) {
         int offset = idxsize-1 - bcorner.pos;
         memcpy((*lcp)->val + bcorner.ptr, &offset, sizeof(int));
         lcp_index_set(lcpsample, bcorner.pos);
         if (offset > 127) lcp_index_set(lcpext, bcorner.pos);
      }
   }

   if (min_depth == 0) {
      lcp_index_set(lcpsample, idxsize-1);
      if (stack_push_lcp(lcp, pl, 0)) return -1;
   }

   free(top);
   free(bottom);

   // Compact LCP.
   stack8_t * stack = *lcp;
   int32_t  val = 0;
   uint64_t p = 0, i = 0;
   while (i < stack->pos) {
      // Save LCP value.
      stack->val[p++] = stack->val[i++];
      // Compress offset.
      memcpy(&val, stack->val + i, sizeof(int32_t));
      if (val < -127) {
         stack->val[p++] = (int8_t)-128;
         stack32_push(ext, val);
      } else if (val > 127) {
         stack->val[p++] = (int8_t)127;
         stack32_push(ext, val);
      } else {
         stack->val[p++] = (int8_t)val - (val > 0);
      }
      i += 4;
   }

   stack->pos = p;
   *lcp = realloc(*lcp, sizeof(stack8_t) + stack->pos*sizeof(int8_t));
   if (*lcp == NULL) return -1;
   (*lcp)->size = p;

   *ext = realloc(*ext, sizeof(stack32_t) + (*ext)->pos*sizeof(int32_t));
   if (*ext == NULL) return -1;
   (*ext)->size = (*ext)->pos;
   return 0;
}


lcp_t
compute_lcp
(
 uint64_t    idxsize,
 int         min_depth,
 int         sar_bits,
 uint64_t  * sar,
 char      * genome
)
{
   // Allocate LCP stack.
   stack8_t * stack = stack8_new(STACK_LCP_INITIAL_SIZE);
   stack32_t * ext  = stack32_new(STACK_LCP_INITIAL_SIZE);
   // Allocate LCP index.
   uint64_t words = (idxsize + LCP_WORD_SIZE - 1)/LCP_WORD_SIZE;
   uint64_t intervals = (words + LCP_MARK_INTERVAL - 1) / LCP_MARK_INTERVAL;
   uint64_t marks = intervals + 1;
   uint64_t index_size = intervals*LCP_MARK_INTERVAL + marks;
   uint64_t * lcp_sample = calloc(index_size,sizeof(uint64_t));
   uint64_t * lcp_extend = calloc(index_size,sizeof(uint64_t));

   naive_lcp(min_depth, idxsize, lcp_sample, lcp_extend, &stack, &ext, sar, sar_bits, genome);

   // LCP sample marks.
   uint64_t abs_pos = 0;
   for (uint64_t i = 0, w = 0; i < words; i++) {
      if (i%LCP_MARK_INTERVAL == 0) lcp_sample[w++] = abs_pos;
      abs_pos += __builtin_popcountl(lcp_sample[w++]);
   }
   lcp_sample[index_size - 1] = abs_pos;

   // LCP ext marks.
   abs_pos = 0;
   for (uint64_t i = 0, w = 0; i < words; i++) {
      if (i%LCP_MARK_INTERVAL == 0) lcp_sample[w++] = abs_pos;
      abs_pos += __builtin_popcountl(lcp_sample[w++]);
   }
   lcp_sample[index_size - 1] = abs_pos;

   lcp_t lcp = {.idx_size = index_size, .idx_sample = lcp_sample, .idx_extend = lcp_extend, .lcp_sample = stack, .lcp_extend = ext};

   // Return stack.
   return lcp;
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

   *genomesize *= 2;

   // Insert wilcard at the end.
   genome[2*div] = 0;

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
         lastbit = lastbit - 64;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         array[++word] = (current & mask) >> (bits - lastbit);
      }
   }

   return word + (lastbit > 0);
}
