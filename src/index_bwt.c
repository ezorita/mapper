#include "index_bwt.h"

// Interface data types.
struct bwt_t {
   size_t      mmap_len;
   void      * mmap_ptr;
   uint64_t    occ_length;
   uint64_t    occ_mark_intv;
   uint64_t    occ_word_size;
   uint64_t    occ_mark_bits;
   uint64_t  * c;
   uint64_t  * occ;
   txt_t     * txt;
};

struct bwtquery_t {
   int64_t    fp;
   int64_t    rp;
   int64_t    sz;
   int64_t    dp;
   bwt_t    * bwt;
};


// Private function headers.
int64_t     get_occ       (int64_t ptr, bwt_t * bwt, int sym);
int         get_occ_all   (int64_t ptr, bwt_t * bwt, int64_t * occ);


// Interface function source.
bwtquery_t *
bwt_new_query
(
  bwt_t  * bwt
)
{
   // Declare variables.
   bwtquery_t * q = NULL;

   // Check arguments.
   error_test_msg(bwt == NULL, "argument 'bwt' is NULL.");
   
   // Alloc bwt query.
   q = malloc(sizeof(bwtquery_t));
   error_test_mem(q);

   // Fill query of depth 0.
   q->fp  = 0;
   q->rp  = 0;
   q->sz  = txt_length(bwt->txt);
   q->dp  = 0;
   q->bwt = bwt;

   return q;

 failure_return:
   free(q);
   return NULL;
}


bwtquery_t *
bwt_dup_query
(
  bwtquery_t  * q
)
{
   // Declare variables.
   bwtquery_t * qdup = NULL;

   // Check arguments.
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   error_test_msg(q->bwt == NULL, "argument 'q->bwt' is NULL.");
   
   qdup = malloc(sizeof(bwtquery_t));
   error_test_mem(qdup);

   // Copy and return pointer.
   memcpy(qdup, q, sizeof(bwtquery_t));
   
   return qdup;

 failure_return:
   free(qdup);
   return NULL;
}


bwtquery_t **
bwt_new_vec
(
  bwt_t  * bwt
)
{
   // Declare variables.
   bwtquery_t ** qv = NULL;
   int32_t sym_cnt = 0;

   // Check arguments.
   error_test_msg(bwt == NULL, "argument 'bwt' is NULL.");

   sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   
   qv = calloc(sym_cnt, sizeof(void *));
   error_test_mem(qv);

   // Initialize queries.
   for (int i = 0; i < sym_cnt; i++) {
      qv[i] = bwt_new_query(bwt);
      error_test(qv[i] == NULL);
   }

   return qv;

 failure_return:
   if (qv != NULL) {
      for (int i = 0; i < sym_cnt; i++) {
         free(qv[i]);
      }
   }
   free(qv);
   return NULL;
}


bwtquery_t **
bwt_dup_vec
(
  bwtquery_t  ** qv
)
{
   // Declare variables.
   bwtquery_t ** qvdup = NULL;
   int32_t sym_cnt = 0;

   // Check arguments.
   error_test_msg(qv == NULL, "argument 'qv' is NULL.");
   error_test_msg(qv[0] == NULL, "vector 'qv' is not allocated.");

   bwt_t * bwt = qv[0]->bwt;
   error_test_msg(bwt == NULL, "missing reference to bwt_t from bwtquery_t.");

   sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Alloc new vector.
   qvdup = calloc(sym_cnt, sizeof(void *));
   error_test_mem(qvdup);

   // Duplicate all queries.
   for (int i = 0; i < sym_cnt; i++) {
      qvdup[i] = bwt_dup_query(qv[i]);
      error_test(qvdup[i] == NULL);
   }

   return qvdup;

 failure_return:
   if (qvdup != NULL) {
      for (int i = 0; i < sym_cnt; i++) {
         free(qvdup[i]);
      }
   }
   free(qvdup);
   return NULL;
}


int
bwt_free_vec
(
  bwtquery_t ** qv
)
{
   // Check arguments.
   error_test_msg(qv == NULL, "argument 'qv' is NULL.");
   error_test_msg(qv[0] == NULL, "vector 'qv' is not allocated.");

   bwt_t * bwt = qv[0]->bwt;
   error_test_msg(bwt == NULL, "missing reference to bwt_t from bwtquery_t.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Free queries.
   for (int i = 0; i < sym_cnt; i++) {
      free(qv[i]);
   }
   
   // Free vector.
   free(qv);

   return 0;

 failure_return:
   return -1;
}


int
bwt_query_all
(
  int           end,
  bwtquery_t  * q,
  bwtquery_t ** qv
)
{
   // Declare variables.
   int64_t * occ_sp = NULL;
   int64_t * occ_ep = NULL;
   int64_t * fp     = NULL;
   int64_t * rp     = NULL;
   int64_t * sz     = NULL;

   // Param check.
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   error_test_msg(qv == NULL, "argument 'qv' is NULL.");

   error_test_msg(end != BWT_QUERY_SUFFIX && end != BWT_QUERY_PREFIX, "invalid 'end' flag.");

   bwt_t * bwt = q->bwt;
   error_test_msg(bwt == NULL, "missing reference to bwt_t from bwtquery_t.");

   sym_t * sym = txt_get_symbols(bwt_get_text(bwt));
   int32_t sym_cnt = sym_count(sym);

   for (int j = 0; j < sym_cnt; j++) {
      error_test_msg(qv[j] == NULL, "found NULL element in bwtquery vector.");
   }

   int64_t q_fp = q->fp;
   int64_t q_rp = q->rp;
   int64_t q_sz = q->sz;
   int64_t q_dp = q->dp;

   // Swap fp/rp for forward extension.
   if (end == BWT_QUERY_SUFFIX) {
      q_fp = q->rp;
      q_rp = q->fp;
   }

   // Alloc temp pointers for each symbol.
   occ_sp = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(occ_sp);

   occ_ep = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(occ_ep);

   fp     = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(fp);

   rp     = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(rp);

   sz     = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(sz);

   // Fetch occ values to update pointers.
   get_occ_all(q_fp - 1, bwt, occ_sp);
   get_occ_all(q_fp + q_sz - 1, bwt, occ_ep);

   // Update forward pointers and result sizes.
   int64_t tot = 0;
   for (int j = 0; j < sym_cnt; j++) {
      fp[j]  = bwt->c[j] + occ_sp[j];
      sz[j]  = occ_ep[j] - occ_sp[j];
      // Sum total span.
      tot += sz[j];
   }

   // Update reverse pointer in complement order.
   rp[sym_complement(0,sym)] = q_rp + (q_sz - tot);
   for (int j = 1; j < sym_cnt; j++) {
      rp[sym_complement(j,sym)] = rp[sym_complement(j-1,sym)] + sz[sym_complement(j-1,sym)];
   }

   if (end == BWT_QUERY_PREFIX) {
      // Store pointers to qv vector.
      for (int j = 0; j < sym_cnt; j++) {
         *(qv[j]) = (bwtquery_t){fp[j], rp[j], sz[j], q_dp+1, bwt};
      }
   } else {
      // Reverse again forward and backward pointers and reversible symbols.
      for (int j  = 0; j < sym_cnt; j++) {
         int32_t sym_comp = sym_complement(j,sym);
         *(qv[j]) = (bwtquery_t){rp[sym_comp], fp[sym_comp], sz[sym_comp], q_dp+1, bwt};
      }
   }

   // Free memory.
   free(occ_sp);
   free(occ_ep);
   free(fp);
   free(rp);
   free(sz);

   return 0;

 failure_return:
   free(occ_sp);
   free(occ_ep);
   free(fp);
   free(rp);
   free(sz);
   return -1;
}


int
bwt_query
(
  int          sym,
  int          end,
  bwtquery_t * q,
  bwtquery_t * qo
)
{
   // Declare variables.
   bwtquery_t ** qv = NULL;
   // Param check.
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   error_test_msg(qo == NULL, "argument 'qo' is NULL.");
   
   bwt_t * bwt = q->bwt;
   error_test_msg(bwt == NULL, "missing reference to bwt_t from bwtquery_t.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   error_test_msg(sym < 0, "sym index must be positive.");
   error_test_msg(sym >= sym_cnt, "sym index out of bounds.");

   // Alloc query vector.
   qv = bwt_new_vec(bwt);
   error_test_mem(qv);

   // Query all symbols (same algorithm complexity).
   error_test(bwt_query_all(end, q, qv) == -1);

   // Copy updated query to q.
   memcpy(qo, qv[sym], sizeof(bwtquery_t));

   // Free query vector.
   bwt_free_vec(qv);
   
   return 0;

 failure_return:
   bwt_free_vec(qv);
   return -1;
}


int
bwt_prefix
(
  int           sym,
  bwtquery_t  * q,
  bwtquery_t  * qo
)
{
   // Check arguments.
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   error_test_msg(qo == NULL, "argument 'qo' is NULL.");
   
   bwt_t * bwt = q->bwt;
   error_test_msg(bwt == NULL, "missing reference to bwt_t from bwtquery_t.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   error_test_msg(sym < 0, "sym index must be positive.");
   error_test_msg(sym >= sym_cnt, "sym index out of bounds.");

   // Update start pointer.
   int64_t sp;
   error_test((sp = get_occ(q->fp - 1, bwt, sym)) == -1);
   sp += bwt->c[sym];
   // Update end pointer.
   int64_t ep;
   error_test((ep = get_occ(q->fp + q->sz - 1, bwt, sym)) == -1);
   ep += bwt->c[sym] -1;

   qo->fp = sp;
   qo->sz = ep - sp + 1;
   qo->rp = -1;
   qo->dp += 1;
   qo->bwt = bwt;

   return 0;

 failure_return:
   return -1;
}

int
bwt_prefix_all
(
  bwtquery_t  * q,
  bwtquery_t ** qv
)
{
   // Declare variables.
   int64_t * occ_sp = NULL;
   int64_t * occ_ep = NULL;

   // Check arguments.
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   error_test_msg(qv == NULL, "argument 'qv' is NULL.");

   bwt_t * bwt = q->bwt;
   error_test_msg(bwt == NULL, "missing reference to bwt_t from bwtquery_t.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   
   for (int j = 0; j < sym_cnt; j++) {
      error_test_msg(qv[j] == NULL, "bwtquery_t vector component uninitialized.");
   }

   int64_t q_fp = q->fp;
   int64_t q_sz = q->sz;
   int64_t q_dp = q->dp;
   
   // Compute occs for all nt in one call (1 cache miss).
   occ_sp = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(occ_sp);
   occ_ep = malloc(sym_cnt * sizeof(int64_t));
   error_test_mem(occ_ep);

   // Compute occ.
   error_test(get_occ_all(q_fp - 1, bwt, occ_sp) == -1);
   error_test(get_occ_all(q_fp + q_sz - 1, bwt, occ_ep) == -1);

   // Update pointers.
   for (int j = 0; j < sym_cnt; j++) {
      qv[j]->fp = bwt->c[j] + occ_sp[j];
      qv[j]->sz = bwt->c[j] + occ_ep[j] - qv[j]->fp;
      qv[j]->sz = (qv[j]->sz < 0 ? 0 : qv[j]->sz);
      qv[j]->dp = q_dp + 1;
      qv[j]->rp = -1;
      qv[j]->bwt = bwt;
   }

   free(occ_sp);
   free(occ_ep);
   return 0;
   
 failure_return:
   free(occ_sp);
   free(occ_ep);
   return -1;
}


bwt_t *
bwt_build
(
  txt_t  * txt,
  sar_t  * sar
)
{
   return bwt_build_opt(txt, sar, BWT_OCC_MARK_INTV_DEF);
}


bwt_t *
bwt_build_opt
(
  txt_t     * txt,
  sar_t     * sar,
  uint64_t    mark_intv
)
{
   // Declare variables.
   bwt_t    * bwt     = NULL;
   uint64_t * occ_abs = NULL;
   uint64_t * occ_tmp = NULL;

   // Word and interval sizes.
   uint64_t word_size = BWT_OCC_WORD_SIZE_DEF;
   uint64_t mark_bits = mark_intv * word_size;

   // Check arguments.
   error_test_msg(txt == NULL, "argument 'txt' is NULL.");
   error_test_msg(sar == NULL, "argument 'sar' is NULL.");
   error_test_msg(mark_intv < 1, "argument 'mark_intv' must be greater than 0.");
   error_test_msg(word_size % 8 != 0, "argument 'word_size' must be multiple of 8.");
   error_test_msg(mark_bits < word_size, "argument 'mark_bits' must be >= 'word_size'.");

   // Alloc bwt structure.
   bwt = malloc(sizeof(bwt_t));
   error_test_mem(bwt);

   // Inherits info from text.
   bwt->mmap_len = 0;
   bwt->mmap_ptr = NULL;
   bwt->occ_mark_intv = mark_intv;
   bwt->occ_word_size = word_size;
   bwt->occ_mark_bits = mark_bits;
   bwt->txt           = txt;
   bwt->occ           = NULL;
   bwt->c             = NULL;

   int64_t text_len   = txt_length(txt);
   int32_t sym_cnt    = sym_count(txt_get_symbols(txt));

   // Words, marks and intervals.
   uint64_t n_intv = (text_len + mark_bits - 1) / mark_bits;
   uint64_t n_word = n_intv * mark_intv;
   uint64_t n_mark = n_intv + 1;

   // Alloc OCC and C structures.
   bwt->occ = malloc((n_word + n_mark) * sym_cnt * sizeof(uint64_t));
   error_test_mem(bwt->occ);
   bwt->c = malloc((sym_cnt+1)*sizeof(uint64_t));
   error_test_mem(bwt->c);

   // Alloc buffers.
   occ_abs = malloc((sym_cnt + 1) * sizeof(uint64_t));
   error_test_mem(occ_abs);
   occ_tmp = malloc((sym_cnt + 1) * sizeof(uint64_t));
   error_test_mem(occ_tmp);

   // Initial values.
   for (int i = 0; i < sym_cnt; i++) {
      occ_abs[i]  = 0;
      occ_tmp[i]  = 0;
      bwt->occ[i] = 0;
   }

   // Compute OCC. (MSB FIRST encoding)
   uint64_t word = sym_cnt, interval = 0;

   for (uint64_t i = 0; i < text_len; i++) {
      // Get symbol at Suffix Array predecessor.
      int64_t txt_ptr = sar_get(i, sar);
      int     sym     = txt_sym(txt_ptr > 0 ? txt_ptr - 1 : text_len - 1, txt);
      // Set bit and update total count.
      occ_tmp[sym] |= 1;
      occ_abs[sym]++;
      // Next word.
      if (i % word_size == word_size - 1) {
         for (int j = 0; j < sym_cnt; j++) {
            bwt->occ[word++] = occ_tmp[j];
            occ_tmp[j] = 0;
         }
         interval++;
         // Write Mark.
         if (interval == mark_intv) {
            for (int j = 0; j < sym_cnt; j++) bwt->occ[word++] = occ_abs[j];
            interval = 0;
         }
      }
      // Shift words.
      for (int j = 0; j < sym_cnt; j++) occ_tmp[j] <<= 1;
   }
   // Shift last word.
   if (text_len % word_size) {
      interval++;
      for (int j = 0; j < sym_cnt; j++)
         bwt->occ[word++] = occ_tmp[j] << (word_size - 1 - (text_len % word_size));
   }
   // Add last mark.
   if (interval > 0) {
      // Fill the last interval with 0.
      for (int i = interval; i < mark_intv; i++)
         for (int j = 0; j < sym_cnt; j++) bwt->occ[word++] = 0;
      // Add mark.
      for (int j = 0; j < sym_cnt; j++) bwt->occ[word++] = occ_abs[j];
   }
   // Set occ size.
   bwt->occ_length = word;

   // Free buffers.
   free(occ_tmp);
   free(occ_abs);

   // C array.

   // Count the wildcards before 'A'.
   bwt->c[0] = txt_wildcard_count(txt);

   // C is the cumulative sum of symbol counts.
   uint64_t * sym_count = bwt->occ + bwt->occ_length - sym_cnt;

   for (int j = 1; j <= sym_cnt; j++) {
      bwt->c[j] = bwt->c[j-1] + sym_count[j-1];
   }

   return bwt;

 failure_return:
   free(occ_tmp);
   free(occ_abs);
   bwt_free(bwt);
   return NULL;
}


void
bwt_free
(
  bwt_t  * bwt
)
{
   if (bwt != NULL) {
      if (bwt->mmap_ptr == NULL) {
         free(bwt->c);
         free(bwt->occ);
      } else {
         munmap(bwt->mmap_ptr, bwt->mmap_len);
      }
      free(bwt);
   }

   return;
}

// Helper functions.
txt_t *
bwt_get_text
(
  bwt_t  * bwt
)
{
   error_test_msg(bwt == NULL, "argument 'bwt' is NULL.");
   return bwt->txt;

 failure_return:
   return NULL;
}


int64_t
bwt_start
(
  bwtquery_t  * q
)
{
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   return q->fp;

 failure_return:
   return -1;

}


int64_t
bwt_rcstart
(
  bwtquery_t  * q
)
{
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   return q->rp;

 failure_return:
   return -1;
}


int64_t
bwt_size
(
  bwtquery_t  * q
)
{
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   return q->sz;

 failure_return:
   return -1;
}


int64_t
bwt_depth
(
  bwtquery_t  * q
)
{
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   return q->dp;

 failure_return:
   return -1;
}


bwt_t *
bwt_get_bwt
(
  bwtquery_t  * q
)
{
   error_test_msg(q == NULL, "argument 'q' is NULL.");
   return q->bwt;

 failure_return:
   return NULL;
}




// I/O functions.
int
bwt_file_write
(
  char   * filename,
  bwt_t  * bwt
)
{
   // Declare variables.
   int fd = -1;

   // Check arguments.
   error_test_msg(filename == NULL, "argument 'filename' is NULL.");
   error_test_msg(bwt == NULL, "argument 'bwt' is NULL.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   
   // Open file.
   fd = creat(filename, 0644);
   error_test_def(fd == -1);

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = BWT_FILE_MAGICNO;

   // Write magic.
   error_test_def(write(fd, &magic, sizeof(uint64_t)) == -1);
   
   // Write occ_len.
   error_test_def(write(fd, (int64_t *)&(bwt->occ_length), sizeof(int64_t)) == -1);

   // Write occ_mark_intv.
   error_test_def(write(fd, (int64_t *)&(bwt->occ_mark_intv), sizeof(int64_t)) == -1);
      
   // Write occ_word_size.
   error_test_def(write(fd, (int64_t *)&(bwt->occ_word_size), sizeof(int64_t)) == -1);

   // Write occ_mark_bits.
   error_test_def(write(fd, (int64_t *)&(bwt->occ_mark_bits), sizeof(int64_t)) == -1);

   // Write C array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint64_t *)bwt->c + e_cnt, (sym_cnt+1 - e_cnt)*sizeof(uint64_t));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(uint64_t);
   } while (e_cnt < sym_cnt+1);

   // Write OCC array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint64_t *)bwt->occ + e_cnt, (bwt->occ_length - e_cnt)*sizeof(uint64_t));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(uint64_t);
   } while (e_cnt < bwt->occ_length);
   
   close(fd);

   return 0;

 failure_return:
   if (fd != -1)
      close(fd);
   return -1;
}


bwt_t *
bwt_file_read
(
  char   * filename,
  txt_t  * txt
)
{
   // Declare variables.
   int fd = -1;
   uint64_t * data = NULL;
   bwt_t * bwt = NULL;

   // Check arguments.
   error_test_msg(filename == NULL, "argument 'filename' is NULL.");
   error_test_msg(txt == NULL, "argument 'txt' is NULL.");

   int32_t sym_cnt = sym_count(txt_get_symbols(txt));

   // Open file.
   fd = open(filename, O_RDONLY);
   error_test_def(fd == -1);

   // Alloc memory.
   bwt = malloc(sizeof(bwt_t));
   error_test_mem(bwt);

   // Set NULL pointers.
   bwt->occ = NULL;
   bwt->c = NULL;

   // Get file len and mmap file.
   struct stat sb;
   fstat(fd, &sb);
   error_test_msg(sb.st_size < 48, "'bwt' file is too small (sb.st_size < 48).");

   data = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
   error_test_def(data == NULL);

   bwt->mmap_len = sb.st_size;
   bwt->mmap_ptr = (void *) data;

   // Read file.
   // Read magic number.
   uint64_t magic = data[0];
   error_test_msg(magic != BWT_FILE_MAGICNO, "unrecognized 'bwt' file format (magicno).");

   bwt->occ_length    = data[1];
   bwt->occ_mark_intv = data[2];
   bwt->occ_word_size = data[3];
   bwt->occ_mark_bits = data[4];

   bwt->c   = data + 5;
   bwt->occ = bwt->c + sym_cnt + 1;
   bwt->txt = txt;

   close(fd);

   return bwt;

 failure_return:
   if (fd != -1)
      close(fd);
   if (data != NULL)
      munmap(data, sb.st_size);
   bwt_free(bwt);
   return NULL;
}



// Private function source.
int
get_occ_all
(
  int64_t    ptr,
  bwt_t    * bwt,
  int64_t  * occ
)
{
   // Declare variables.
   uint64_t * offset = NULL;

   // Check arguments.
   error_test_msg(bwt == NULL, "argument 'bwt' is NULL.");
   error_test_msg(occ == NULL, "argument 'occ' is NULL.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Depth 0 pointer, return full span.
   if (ptr == -1) {
      for (int j = 0; j < sym_cnt; j++) occ[j] = bwt->occ[j];
      return 0;
   }

   // Compute word, bit and absolute marker positions.
   int64_t wrdnum = ptr/bwt->occ_word_size;
   int64_t wrdptr = (wrdnum + wrdnum/bwt->occ_mark_intv + 1) * sym_cnt;
   int64_t mrkptr = (((wrdnum + bwt->occ_mark_intv/2)/bwt->occ_mark_intv) * (bwt->occ_mark_intv+1)) * sym_cnt;
   int64_t bit    = ptr % bwt->occ_word_size;

   offset = calloc(sym_cnt, sizeof(uint64_t));
   error_test_mem(offset);

   if (wrdptr > mrkptr) {
      // Sum bit offsets.
      uint64_t i = mrkptr + sym_cnt;
      while (i < wrdptr) 
         for (int j = 0; j < sym_cnt; j++)
            offset[j] += __builtin_popcountl(bwt->occ[i++]);
      // Sum partial word.
      for (int j = 0; j < sym_cnt; j++)
         offset[j] += __builtin_popcountl(bwt->occ[wrdptr++] >> (bwt->occ_word_size - 1 - bit));
      // Return sum.
      for (int j = 0; j < sym_cnt; j++)
         occ[j] = bwt->occ[mrkptr++] + offset[j];
   } else {
      // Sum partial word.
      if (bit < bwt->occ_word_size - 1)
         for (int j = 0; j < sym_cnt; j++)
            offset[j] += __builtin_popcountl(bwt->occ[wrdptr++] << (bit+1));
      else wrdptr += sym_cnt;
      // Sum bit offsets.
      while (wrdptr < mrkptr)
         for (int j = 0; j < sym_cnt; j++)
            offset[j] += __builtin_popcountl(bwt->occ[wrdptr++]);
      // Return subtraction.
      for (int j = 0; j < sym_cnt; j++)
         occ[j] = bwt->occ[mrkptr++] - offset[j];
   }

   free(offset);
   return 0;
   
 failure_return:
   free(offset);
   return -1;
}



int64_t
get_occ
(
  int64_t    ptr,
  bwt_t    * bwt,
  int        sym
)
{
   // Check arguments.
   error_test_msg(bwt == NULL, "argument 'bwt' is NULL.");

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   error_test_msg(sym < 0, "argument 'sym' must be positive.");
   error_test_msg(sym >= sym_cnt, "index 'sym' out of bounds.");

   // Depth 0 pointer, return full span.
   if (ptr == -1) return bwt->occ[sym];
   // Compute word, bit and absolute marker positions.
   int64_t wrdnum = ptr / bwt->occ_word_size;
   int64_t wrdptr = (wrdnum + wrdnum/bwt->occ_mark_intv + 1) * sym_cnt + sym;
   int64_t mrkptr = (((wrdnum + bwt->occ_mark_intv/2)/bwt->occ_mark_intv) * (bwt->occ_mark_intv+1)) * sym_cnt + sym;
   int64_t bit    = ptr % bwt->occ_word_size;

   uint64_t occ = bwt->occ[mrkptr];
   if (wrdptr > mrkptr) {
      int64_t  offset = 0;
      // Sum bit offsets.
      for (uint64_t i = mrkptr + sym_cnt; i < wrdptr; i+= sym_cnt)
         offset += __builtin_popcountl(bwt->occ[i]);
      // Sum partial word.
      offset += __builtin_popcountl(bwt->occ[wrdptr] >> (bwt->occ_word_size - 1 - bit));
      // Returm sum.
      occ += offset;
   } else {
      int64_t  offset = 0;
      // Sum partial word.
      if (bit < bwt->occ_word_size - 1)
         offset += __builtin_popcountl(bwt->occ[wrdptr] << (bit+1));
      // Sum bit offsets.
      for (uint64_t i = wrdptr + sym_cnt; i < mrkptr; i+= sym_cnt)
         offset += __builtin_popcountl(bwt->occ[i]);
      // Return subtraction.
      occ -= offset;
   }
   return (int64_t)occ;

 failure_return:
   return -1;
}
