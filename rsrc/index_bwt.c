#include "index_bwt.h"

// Interface data types.
struct bwt_t {
   uint64_t   occ_length;
   uint64_t   occ_mark_intv;
   uint64_t   occ_word_size;
   uint64_t   occ_mark_bits;
   uint64_t * c;
   uint64_t * occ;
   txt_t    * txt;
};

struct bwtquery_t {
   int64_t   fp;
   int64_t   rp;
   int64_t   sz;
   int64_t   dp;
   bwt_t   * bwt;
};


// Private function headers.
uint64_t    get_occ       (int64_t ptr, bwt_t * bwt, int sym);
int         get_occ_all   (int64_t ptr, bwt_t * bwt, int64_t * occ);


// Interface function source.
bwtquery_t *
bwt_new_query
(
  bwt_t  * bwt
)
{
   // Check arguments.
   if (bwt == NULL)
      return NULL;
   
   // Alloc bwt query.
   bwtquery_t * q = malloc(sizeof(bwtquery_t));
   if (q == NULL)
      return NULL;

   // Fill query of depth 0.
   q->fp  = 0;
   q->rp  = 0;
   q->sz  = txt_length(bwt->txt);
   q->dp  = 0;
   q->bwt = bwt;

   return q;
}


bwtquery_t *
bwt_dup_query
(
  bwtquery_t  * q
)
{
   // Check arguments.
   if (q == NULL)
      return NULL;

   if (q->bwt == NULL)
      return NULL;
   
   bwtquery_t * qdup = malloc(sizeof(bwtquery_t));
   if (qdup == NULL)
      return NULL;

   // Copy and return pointer.
   memcpy(qdup, q, sizeof(bwtquery_t));
   
   return qdup;
}


bwtquery_t **
bwt_new_vec
(
  bwt_t  * bwt
)
{
   // Check arguments.
   if (bwt == NULL)
      return NULL;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   
   bwtquery_t ** qv = malloc(sym_cnt*sizeof(void *));
   if (qv == NULL)
      return NULL;

   // Initialize queries.
   for (int i = 0; i < sym_cnt; i++) {
      qv[i] = bwt_new_query(bwt);
      if (qv[i] == NULL)
         return NULL;
   }

   return qv;
}


bwtquery_t **
bwt_dup_vec
(
  bwtquery_t  ** qv
)
{
   // Check arguments.
   if (qv == NULL || qv[0] == NULL)
      return NULL;

   bwt_t * bwt = qv[0]->bwt;

   if (bwt == NULL)
      return NULL;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Alloc new vector.
   bwtquery_t ** qvdup = malloc(sym_cnt*sizeof(void *));
   if (qvdup == NULL)
      return NULL;
   

   // Duplicate all queries.
   for (int i = 0; i < sym_cnt; i++) {
      qvdup[i] = bwt_dup_query(qv[i]);
      if (qvdup[i] == NULL)
         return NULL;
   }

   return qvdup;
}


int
bwt_free_vec
(
  bwtquery_t ** qv
)
{
   // Check arguments.
   if (qv == NULL)
      return -1;
   if (qv[0] == NULL)
      return -1;

   bwt_t * bwt = qv[0]->bwt;
   if (bwt == NULL)
      return -1;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Free queries.
   for (int i = 0; i < sym_cnt; i++) {
      free(qv[i]);
   }
   
   // Free vector.
   free(qv);

   return 0;
}


int
bwt_query_all
(
  int           end,
  bwtquery_t  * q,
  bwtquery_t ** qv
)
{
   // Param check.
   if (q == NULL || qv == NULL)
      return -1;

   if (end != BWT_QUERY_SUFFIX && end != BWT_QUERY_PREFIX)
      return -1;

   bwt_t * bwt = q->bwt;
   if (bwt == NULL)
      return -1;

   sym_t * sym = txt_get_symbols(bwt_get_text(bwt));
   int32_t sym_cnt = sym_count(sym);

   for (int j = 0; j < sym_cnt; j++) {
      if (qv[j] == NULL)
         return -1;
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
   int64_t * occ_sp = malloc(sym_cnt * sizeof(int64_t));
   int64_t * occ_ep = malloc(sym_cnt * sizeof(int64_t));
   int64_t * fp     = malloc(sym_cnt * sizeof(int64_t));
   int64_t * rp     = malloc(sym_cnt * sizeof(int64_t));
   int64_t * sz     = malloc(sym_cnt * sizeof(int64_t));

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
   // Param check.
   if (q == NULL || qo == NULL)
      return -1;
   
   bwt_t * bwt = q->bwt;
   if (bwt == NULL || sym < 0)
      return -1;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   if (sym >= sym_cnt)
      return -1;

   // Alloc query vector.
   bwtquery_t ** qv = bwt_new_vec(bwt);
   if (qv == NULL)
      return -1;

   // Query all symbols (same algorithm complexity).
   if (bwt_query_all(end, q, qv) < 0) {
      bwt_free_vec(qv);
      return -1;
   }

   // Copy updated query to q.
   memcpy(qo, qv[sym], sizeof(bwtquery_t));

   // Free query vector.
   bwt_free_vec(qv);
   
   return 0;
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
   if (q == NULL || qo == NULL)
      return -1;
   
   bwt_t * bwt = q->bwt;
   if (bwt == NULL || sym < 0)
      return -1;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   if (sym >= sym_cnt)
      return -1;

   // Update start pointer.
   int64_t sp;
   if ((sp = get_occ(q->fp - 1, bwt, sym)) < 0)
      return -1;
   sp += bwt->c[sym];
   // Update end pointer.
   int64_t ep;
   if ((ep = get_occ(q->fp + q->sz - 1, bwt, sym)) < 0)
      return -1;
   ep += bwt->c[sym] -1;

   qo->fp = sp;
   qo->sz = ep - sp + 1;
   qo->rp = -1;
   qo->dp += 1;
   qo->bwt = bwt;

   return 0;
}

int
bwt_prefix_all
(
  bwtquery_t  * q,
  bwtquery_t ** qv
)
{
   // Check arguments.
   if (q == NULL || qv == NULL)
      return -1;

   bwt_t * bwt = q->bwt;
   if (bwt == NULL)
      return -1;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   
   for (int j = 0; j < sym_cnt; j++) {
      if (qv[j] == NULL)
         return -1;
   }

   int64_t q_fp = q->fp;
   int64_t q_sz = q->sz;
   int64_t q_dp = q->dp;
   
   // Compute occs for all nt in one call (1 cache miss).
   int64_t * occ_sp = malloc(sym_cnt * sizeof(int64_t));
   int64_t * occ_ep = malloc(sym_cnt * sizeof(int64_t));

   // Compute occ.
   if (get_occ_all(q_fp - 1, bwt, occ_sp)) {
      free(occ_sp);
      free(occ_ep);
      return -1;
   }
   if (get_occ_all(q_fp + q_sz - 1, bwt, occ_ep)) {
      free(occ_sp);
      free(occ_ep);
      return -1;
   }
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
   uint64_t word_size = BWT_OCC_WORD_SIZE_DEF;
   uint64_t mark_bits = mark_intv * word_size;

   // Check arguments.
   if (!txt || !sar)
      return NULL;

   if (mark_intv < 1 || word_size % 8 != 0 || mark_bits < word_size)
      return NULL;

   // Alloc bwt structure.
   bwt_t * bwt = malloc(sizeof(bwt_t));
   if (bwt == NULL)
      return NULL;

   // Inherits info from text.
   bwt->occ_mark_intv = mark_intv;
   bwt->occ_word_size = word_size;
   bwt->occ_mark_bits = mark_bits;
   bwt->txt           = txt;

   int64_t text_len   = txt_length(txt);
   int32_t sym_cnt    = sym_count(txt_get_symbols(txt));

   // Words, marks and intervals.
   uint64_t n_intv = (text_len + mark_bits - 1) / mark_bits;
   uint64_t n_word = n_intv * mark_intv;
   uint64_t n_mark = n_intv + 1;

   // Alloc OCC and C structures.
   bwt->occ = malloc((n_word + n_mark) * sym_cnt * sizeof(uint64_t));
   if (bwt->occ == NULL) {
      free(bwt);
      return NULL;
   }

   bwt->c = malloc((sym_cnt+1)*sizeof(uint64_t));
   if (bwt->c == NULL) {
      free(bwt->occ);
      free(bwt);
      return NULL;
   }

   // Alloc buffers.
   uint64_t * occ_abs = malloc((sym_cnt + 1) * sizeof(uint64_t));
   uint64_t * occ_tmp = malloc((sym_cnt + 1) * sizeof(uint64_t));

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
}


void
bwt_free
(
  bwt_t  * bwt
)
{
   if (bwt != NULL) {
      free(bwt->c);
      free(bwt->occ);
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
   // Check arguments.
   if (bwt == NULL)
      return NULL;
   
   return bwt->txt;
}


int64_t
bwt_start
(
  bwtquery_t  * q
)
{
   if (q == NULL)
      return -1;
   else
      return q->fp;
}


int64_t
bwt_rcstart
(
  bwtquery_t  * q
)
{
   if (q == NULL)
      return -1;
   else
      return q->rp;
}


int64_t
bwt_size
(
  bwtquery_t  * q
)
{
   if (q == NULL)
      return -1;
   else
      return q->sz;
}


int64_t
bwt_depth
(
  bwtquery_t  * q
)
{
   if (q == NULL)
      return -1;
   else
      return q->dp;
}


bwt_t *
bwt_get_bwt
(
  bwtquery_t  * q
)
{
   if (q == NULL)
      return NULL;
   else
      return q->bwt;
}




// I/O functions.
int
bwt_file_write
(
  char   * filename,
  bwt_t  * bwt
)
{
   // Check arguments.
   if (filename == NULL || bwt == NULL)
      return -1;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   
   // Open file.
   int fd = creat(filename, 0644);
   if (fd == -1)
      return -1;

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = BWT_FILE_MAGICNO;

   // Write magic.
   if (write(fd, &magic, sizeof(uint64_t)) == -1)
      goto close_and_error;
   
   // Write occ_len.
   if (write(fd, (int64_t *)&(bwt->occ_length), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write occ_mark_intv.
   if (write(fd, (int64_t *)&(bwt->occ_mark_intv), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write occ_word_size.
   if (write(fd, (int64_t *)&(bwt->occ_word_size), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write occ_mark_bits.
   if (write(fd, (int64_t *)&(bwt->occ_mark_bits), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write C array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint64_t *)bwt->c + e_cnt, (sym_cnt+1 - e_cnt)*sizeof(uint64_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(uint64_t);
   } while (e_cnt < sym_cnt+1);

   // Write OCC array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint64_t *)bwt->occ + e_cnt, (bwt->occ_length - e_cnt)*sizeof(uint64_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(uint64_t);
   } while (e_cnt < bwt->occ_length);
   
   close(fd);
   return 0;

 close_and_error:
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
   /* WHAT ABOUT USING MMAP INSTEAD? */
   // Check arguments.
   if (filename == NULL || txt == NULL)
      return NULL;

   int32_t sym_cnt = sym_count(txt_get_symbols(txt));

   // Open file.
   int fd = open(filename, O_RDONLY);
   if (fd == -1)
      return NULL;

   // Alloc memory.
   bwt_t * bwt = malloc(sizeof(bwt_t));
   if (bwt == NULL)
      return NULL;
   // Set NULL pointers.
   bwt->occ = NULL;
   bwt->c = NULL;
   
   // Read file.
   uint64_t magic;
   ssize_t b_cnt;
   ssize_t e_cnt;

   // Read magic number.
   if (read(fd, &magic, sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;
   if (magic != BWT_FILE_MAGICNO)
      goto free_and_return;

   // Read 'occ_length'.
   if (read(fd, &(bwt->occ_length), sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;

   // Read 'occ_mark_intv'.
   if (read(fd, &(bwt->occ_mark_intv), sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;

   // Read 'occ_word_size'.
   if (read(fd, &(bwt->occ_word_size), sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;

   // Read 'occ_mark_bits'.
   if (read(fd, &(bwt->occ_mark_bits), sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;

   // Read C array.
   bwt->c = malloc((sym_cnt+1)*sizeof(uint64_t));
   if (bwt->c == NULL)
      goto free_and_return;

   e_cnt = 0;
   do {
      b_cnt = read(fd, (uint64_t *)bwt->c + e_cnt, (sym_cnt+1 - e_cnt) * sizeof(uint64_t));
      if (b_cnt == -1)
         goto free_and_return;
      e_cnt += b_cnt / sizeof(uint64_t);
   } while (e_cnt < sym_cnt+1);

   // Read OCC array.
   bwt->occ = malloc((bwt->occ_length)*sizeof(uint64_t));
   if (bwt->occ == NULL)
      goto free_and_return;

   e_cnt = 0;
   do {
      b_cnt = read(fd, (uint64_t *)bwt->occ + e_cnt, (bwt->occ_length - e_cnt) * sizeof(uint64_t));
      if (b_cnt == -1)
         goto free_and_return;
      e_cnt += b_cnt / sizeof(uint64_t);
   } while (e_cnt < bwt->occ_length);

   bwt->txt = txt;

   return bwt;

 free_and_return:
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
   // Check arguments.
   if (bwt == NULL || occ == NULL)
      return -1;

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

   uint64_t * offset = calloc(sym_cnt, sizeof(uint64_t));
   if (offset == NULL) {
      return 1;
   }

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
}



uint64_t
get_occ
(
  int64_t    ptr,
  bwt_t    * bwt,
  int        sym
)
{
   // Check arguments.
   if (bwt == NULL)
      return -1;

   int32_t sym_cnt = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   if (sym < 0 || sym >= sym_cnt)
      return -1;

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
   return occ;
}
