#include "index_bwt.h"

// Interface data types.
struct bwt_t {
   int32_t    sym_cnt;
   int32_t    rev_cnt;
   uint64_t   txt_length;
   uint64_t   occ_length;
   uint64_t   occ_mark_intv;
   uint64_t   occ_word_size;
   uint64_t   occ_mark_bits;
   uint64_t * c;
   uint64_t * occ;
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
   q->sz  = bwt->txt_length;
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
   
   bwtquery_t ** qv = malloc(bwt->sym_cnt*sizeof(void *));
   if (qv == NULL)
      return NULL;

   // Initialize queries.
   for (int i = 0; i < bwt->sym_cnt; i++) {
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

   // Alloc new vector.
   bwtquery_t ** qvdup = malloc(bwt->sym_cnt*sizeof(void *));
   if (qvdup == NULL)
      return NULL;
   

   // Duplicate all queries.
   for (int i = 0; i < bwt->sym_cnt; i++) {
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

   // Free queries.
   for (int i = 0; i < bwt->sym_cnt; i++) {
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

   bwt_t * bwt = q->bwt;
   if (bwt == NULL)
      return -1;

   for (int j = 0; j < bwt->sym_cnt; j++) {
      if (qv[j] == NULL)
         return -1;
   }
      
   // Alloc temp pointers for each symbol.
   int64_t * occ_sp = malloc(bwt->sym_cnt * sizeof(int64_t));
   int64_t * occ_ep = malloc(bwt->sym_cnt * sizeof(int64_t));
   int64_t * fp     = malloc(bwt->sym_cnt * sizeof(int64_t));
   int64_t * rp     = malloc(bwt->sym_cnt * sizeof(int64_t));
   int64_t * sz     = malloc(bwt->sym_cnt * sizeof(int64_t));

   // Fetch occ values to update pointers.
   get_occ_all(q->fp - 1, bwt, occ_sp);
   get_occ_all(q->fp + q->sz - 1, bwt, occ_ep);

   // Update forward pointers and result sizes.
   int64_t tot = 0;
   for (int j = 0; j < bwt->sym_cnt; j++) {
      fp[j]  = bwt->c[j] + occ_sp[j];
      sz[j]  = occ_ep[j] - occ_sp[j];
      // Sum total span.
      tot += sz[j];
   }

   // Since wildcards '$' are the lexicographically smallest they will reduce
   // the interval size of the first symbol of fp, hence impacting on the last reversible
   // element of rp.
   if (bwt->rev_cnt > 1) {
      rp[bwt->rev_cnt-1] = q->rp + (q->sz - tot);
      for (int j = bwt->rev_cnt-2; j >= 0; j--) rp[j] = rp[j+1] + sz[j+1];
      // The first reversible symbol goes after the first fp symbol (hence the last in rp).
      rp[bwt->rev_cnt] = rp[0] + sz[0];
   }
   // If there are no reversible symbols, the wildcards impact directly on the first symbol.
   else {
      rp[0] = q->rp + (q->sz - tot);
   }

   // The rest of non-reversible symbols are updated sequentially.
   if (bwt->sym_cnt > bwt->rev_cnt) {
      for (int j = bwt->rev_cnt + 1; j < bwt->sym_cnt; j++) {
         rp[j] = rp[j-1] + sz[j-1];
      }
   }

   if (end == BWT_QUERY_PREFIX) {
      // Store pointers to qv vector.
      for (int j = 0; j < bwt->sym_cnt; j++) {
         *(qv[j]) = (bwtquery_t){fp[j], rp[j], sz[j], q->dp+1};
      }
   } else {
      // Reverse again forward and backward pointers and reversible symbols.
      int nr = bwt->rev_cnt-1;
      int j  = 0;
      for (; j <= nr; j++) {
         *(qv[j]) = (bwtquery_t){rp[nr-j], fp[nr-j], sz[nr-j], q->dp+1};
      }
      for (; j < bwt->sym_cnt; j++) {
         *(qv[j]) = (bwtquery_t){rp[j], fp[j], sz[j], q->dp+1};
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
 bwtquery_t * q
)
{
   // Param check.
   if (q == NULL)
      return -1;
   
   bwt_t * bwt = q->bwt;
   if (bwt == NULL || sym < 0 || sym >= bwt->sym_cnt)
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
   memcpy(q, qv[sym], sizeof(bwtquery_t));

   // Free query vector.
   bwt_free_vec(qv);
   
   return 0;
}


int
bwt_prefix
(
 int           sym,
 bwtquery_t  * q
)
{
   // Check arguments.
   if (q == NULL)
      return -1;
   
   bwt_t * bwt = q->bwt;
   if (bwt == NULL || sym < 0 || sym >= bwt->sym_cnt)
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

   q->fp = sp;
   q->sz = ep - sp + 1;
   q->rp = -1;
   q->dp += 1;

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
   
   for (int j = 0; j < bwt->sym_cnt; j++) {
      if (qv[j] == NULL)
         return -1;
   }
   
   // Compute occs for all nt in one call (1 cache miss).
   int64_t * occ_sp = malloc(bwt->sym_cnt * sizeof(int64_t));
   int64_t * occ_ep = malloc(bwt->sym_cnt * sizeof(int64_t));

   // Compute occ.
   if (get_occ_all(q->fp - 1, bwt, occ_sp)) {
      free(occ_sp);
      free(occ_ep);
      return -1;
   }
   if (get_occ_all(q->fp + q->sz - 1, bwt, occ_ep)) {
      free(occ_sp);
      free(occ_ep);
      return -1;
   }
   // Update pointers.
   for (int j = 0; j < bwt->sym_cnt; j++) {
      qv[j]->fp = bwt->c[j] + occ_sp[j];
      qv[j]->sz = bwt->c[j] + occ_ep[j] - qv[j]->fp;
      qv[j]->sz = (qv[j]->sz < 0 ? 0 : qv[j]->sz);
      qv[j]->dp = q->dp + 1;
      qv[j]->rp = -1;
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
   return bwt_build_opt(txt, sar, BWT_OCC_MARK_INTV_DEF, BWT_OCC_WORD_SIZE_DEF, BWT_OCC_MARK_BITS_DEF);
}


bwt_t *
bwt_build_opt
(
  txt_t     * txt,
  sar_t     * sar,
  uint64_t    mark_intv,
  uint64_t    word_size,
  uint64_t    mark_bits
)
{
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
   bwt->sym_cnt       = txt_sym_count(txt);
   bwt->rev_cnt       = txt_rev_count(txt);
   bwt->txt_length    = txt_length(txt);
   bwt->occ_mark_intv = mark_intv;
   bwt->occ_word_size = word_size;
   bwt->occ_mark_bits = mark_bits;

   // Words, marks and intervals.
   uint64_t n_intv = (bwt->txt_length + mark_bits - 1) / mark_bits;
   uint64_t n_word = n_intv * mark_intv;
   uint64_t n_mark = n_intv + 1;

   // Alloc OCC and C structures.
   bwt->occ = malloc((n_word + n_mark) * bwt->sym_cnt * sizeof(uint64_t));
   if (bwt->occ == NULL) 
      return NULL;

   bwt->c = malloc((bwt->sym_cnt+1)*sizeof(uint64_t));
   if (bwt->c == NULL)
      return NULL;

   // Alloc buffers.
   uint64_t * occ_abs = malloc((bwt->sym_cnt + 1) * sizeof(uint64_t));
   uint64_t * occ_tmp = malloc((bwt->sym_cnt + 1) * sizeof(uint64_t));

   // Initial values.
   for (int i = 0; i < bwt->sym_cnt; i++) {
      occ_abs[i]  = 0;
      occ_tmp[i]  = 0;
      bwt->occ[i] = 0;
   }

   // Compute OCC. (MSB FIRST encoding)
   uint64_t word = bwt->sym_cnt, interval = 0;

   for (uint64_t i = 0; i < bwt->txt_length; i++) {
      // Get symbol at Suffix Array predecessor.
      int64_t txt_ptr = sar_get(i, sar);
      int     sym     = txt_get_sym(txt_ptr > 0 ? txt_ptr - 1 : bwt->txt_length - 1, txt);
      // Set bit and update total count.
      occ_tmp[sym] |= 1;
      occ_abs[sym]++;
      // Next word.
      if (i % word_size == word_size - 1) {
         for (int j = 0; j < bwt->sym_cnt; j++) {
            bwt->occ[word++] = occ_tmp[j];
            occ_tmp[j] = 0;
         }
         interval++;
         // Write Mark.
         if (interval == mark_intv) {
            for (int j = 0; j < bwt->sym_cnt; j++) bwt->occ[word++] = occ_abs[j];
            interval = 0;
         }
      }
      // Shift words.
      for (int j = 0; j < bwt->sym_cnt; j++) occ_tmp[j] <<= 1;
   }
   // Shift last word.
   if (bwt->txt_length % word_size) {
      interval++;
      for (int j = 0; j < bwt->sym_cnt; j++)
         bwt->occ[word++] = occ_tmp[j] << (word_size - 1 - (bwt->txt_length % word_size));
   }
   // Add last mark.
   if (interval > 0) {
      // Fill the last interval with 0.
      for (int i = interval; i < mark_intv; i++)
         for (int j = 0; j < bwt->sym_cnt; j++) bwt->occ[word++] = 0;
      // Add mark.
      for (int j = 0; j < bwt->sym_cnt; j++) bwt->occ[word++] = occ_abs[j];
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
   uint64_t * sym_count = bwt->occ + bwt->occ_length - bwt->sym_cnt;

   for (int j = 1; j <= bwt->sym_cnt; j++) {
      bwt->c[j] = bwt->c[j-1] + sym_count[j-1];
   }

   return bwt;
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

   // Depth 0 pointer, return full span.
   if (ptr == -1) {
      for (int j = 0; j < bwt->sym_cnt; j++) occ[j] = bwt->occ[j];
      return 0;
   }

   // Compute word, bit and absolute marker positions.
   int64_t wrdnum = ptr/bwt->occ_word_size;
   int64_t wrdptr = (wrdnum + wrdnum/bwt->occ_mark_intv + 1) * bwt->sym_cnt;
   int64_t mrkptr = (((wrdnum + bwt->occ_mark_intv/2)/bwt->occ_mark_intv) * (bwt->occ_mark_intv+1)) * bwt->sym_cnt;
   int64_t bit    = ptr % bwt->occ_word_size;

   uint64_t * offset = calloc(bwt->sym_cnt, sizeof(uint64_t));
   if (offset == NULL) {
      return 1;
   }

   if (wrdptr > mrkptr) {
      // Sum bit offsets.
      uint64_t i = mrkptr + bwt->sym_cnt;
      while (i < wrdptr) 
         for (int j = 0; j < bwt->sym_cnt; j++)
            offset[j] += __builtin_popcountl(bwt->occ[i++]);
      // Sum partial word.
      for (int j = 0; j < bwt->sym_cnt; j++)
         offset[j] += __builtin_popcountl(bwt->occ[wrdptr++] >> (bwt->occ_word_size - 1 - bit));
      // Return sum.
      for (int j = 0; j < bwt->sym_cnt; j++)
         occ[j] = bwt->occ[mrkptr++] + offset[j];
   } else {
      // Sum partial word.
      if (bit < bwt->occ_word_size - 1)
         for (int j = 0; j < bwt->sym_cnt; j++)
            offset[j] += __builtin_popcountl(bwt->occ[wrdptr++] << (bit+1));
      else wrdptr += bwt->sym_cnt;
      // Sum bit offsets.
      while (wrdptr < mrkptr)
         for (int j = 0; j < bwt->sym_cnt; j++)
            offset[j] += __builtin_popcountl(bwt->occ[wrdptr++]);
      // Return subtraction.
      for (int j = 0; j < bwt->sym_cnt; j++)
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
   
   if (sym < 0 || sym >= bwt->sym_cnt)
      return -1;

   // Depth 0 pointer, return full span.
   if (ptr == -1) return bwt->occ[sym];
   // Compute word, bit and absolute marker positions.
   int64_t wrdnum = ptr / bwt->occ_word_size;
   int64_t wrdptr = (wrdnum + wrdnum/bwt->occ_mark_intv + 1) * bwt->sym_cnt + sym;
   int64_t mrkptr = (((wrdnum + bwt->occ_mark_intv/2)/bwt->occ_mark_intv) * (bwt->occ_mark_intv+1))*bwt->sym_cnt + sym;
   int64_t bit    = ptr % bwt->occ_word_size;

   uint64_t occ = bwt->occ[mrkptr];
   if (wrdptr > mrkptr) {
      int64_t  offset = 0;
      // Sum bit offsets.
      for (uint64_t i = mrkptr + bwt->sym_cnt; i < wrdptr; i+= bwt->sym_cnt)
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
      for (uint64_t i = wrdptr + bwt->sym_cnt; i < mrkptr; i+= bwt->sym_cnt)
         offset += __builtin_popcountl(bwt->occ[i]);
      // Return subtraction.
      occ -= offset;
   }
   return occ;
}
