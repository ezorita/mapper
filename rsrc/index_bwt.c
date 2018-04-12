#include "index_bwt.h"

// Interface data types.
struct bwt_t {
   int32_t    n_symbols;
   int32_t    rev_symbols; // How many symbols are reversible. (4 in DNA)
   uint64_t   seq_length;
   uint64_t   occ_length;
   uint64_t   occ_mark_int;
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
   q->sz  = bwt->seq_length;
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
   
   bwtquery_t ** qv = malloc(bwt->n_symbols*sizeof(void *));
   if (qv == NULL)
      return NULL;

   // Initialize queries.
   for (int i = 0; i < bwt->n_symbols; i++) {
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
   bwtquery_t ** qvdup = malloc(bwt->n_symbols*sizeof(void *));
   if (qvdup == NULL)
      return NULL;
   

   // Duplicate all queries.
   for (int i = 0; i < bwt->n_symbols; i++) {
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
   for (int i = 0; i < bwt->n_symbols; i++) {
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

   for (int j = 0; j < bwt->n_symbols; j++) {
      if (qv[j] == NULL)
         return -1;
   }
      
   // Alloc temp pointers for each symbol.
   int64_t * occ_sp = malloc(bwt->n_symbols * sizeof(int64_t));
   int64_t * occ_ep = malloc(bwt->n_symbols * sizeof(int64_t));
   int64_t * fp     = malloc(bwt->n_symbols * sizeof(int64_t));
   int64_t * rp     = malloc(bwt->n_symbols * sizeof(int64_t));
   int64_t * sz     = malloc(bwt->n_symbols * sizeof(int64_t));

   // Fetch occ values to update pointers.
   get_occ_all(q->fp - 1, bwt, occ_sp);
   get_occ_all(q->fp + q->sz - 1, bwt, occ_ep);

   // Update forward pointers and result sizes.
   int64_t tot = 0;
   for (int j = 0; j < bwt->n_symbols; j++) {
      fp[j]  = bwt->c[j] + occ_sp[j];
      sz[j]  = occ_ep[j] - occ_sp[j];
      // Sum total span.
      tot += sz[j];
   }

   // Since wildcards '$' are the lexicographically smallest they will reduce
   // the interval size of the first symbol of fp, hence impacting on the last reversible
   // element of rp.
   if (bwt->rev_symbols > 1) {
      rp[bwt->rev_symbols-1] = q->rp + (q->sz - tot);
      for (int j = bwt->rev_symbols-2; j >= 0; j--) rp[j] = rp[j+1] + sz[j+1];
      // The first reversible symbol goes after the first fp symbol (hence the last in rp).
      rp[bwt->rev_symbols] = rp[0] + sz[0];
   }
   // If there are no reversible symbols, the wildcards impact directly on the first symbol.
   else {
      rp[0] = q->rp + (q->sz - tot);
   }

   // The rest of non-reversible symbols are updated sequentially.
   if (bwt->n_symbols > bwt->rev_symbols) {
      for (int j = bwt->rev_symbols + 1; j < bwt->n_symbols; j++) {
         rp[j] = rp[j-1] + sz[j-1];
      }
   }

   if (end == BWT_QUERY_PREFIX) {
      // Store pointers to qv vector.
      for (int j = 0; j < bwt->n_symbols; j++) {
         *(qv[j]) = (bwtquery_t){fp[j], rp[j], sz[j], q->dp+1};
      }
   } else {
      // Reverse again forward and backward pointers and reversible symbols.
      int nr = bwt->rev_symbols-1;
      int j  = 0;
      for (; j <= nr; j++) {
         *(qv[j]) = (bwtquery_t){rp[nr-j], fp[nr-j], sz[nr-j], q->dp+1};
      }
      for (; j < bwt->n_symbols; j++) {
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
   if (bwt == NULL || sym < 0 || sym >= bwt->n_symbols)
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
   if (bwt == NULL || sym < 0 || sym >= bwt->n_symbols)
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
   
   for (int j = 0; j < bwt->n_symbols; j++) {
      if (qv[j] == NULL)
         return -1;
   }
   
   // Compute occs for all nt in one call (1 cache miss).
   int64_t * occ_sp = malloc(bwt->n_symbols * sizeof(int64_t));
   int64_t * occ_ep = malloc(bwt->n_symbols * sizeof(int64_t));

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
   for (int j = 0; j < bwt->n_symbols; j++) {
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
      for (int j = 0; j < bwt->n_symbols; j++) occ[j] = bwt->occ[j];
      return 0;
   }

   // Compute word, bit and absolute marker positions.
   int64_t wrdnum = ptr/bwt->occ_word_size;
   int64_t wrdptr = (wrdnum + wrdnum/bwt->occ_mark_int + 1) * bwt->n_symbols;
   int64_t mrkptr = (((wrdnum + bwt->occ_mark_int/2)/bwt->occ_mark_int) * (bwt->occ_mark_int+1)) * bwt->n_symbols;
   int64_t bit    = ptr % bwt->occ_word_size;

   uint64_t * offset = calloc(bwt->n_symbols, sizeof(uint64_t));
   if (offset == NULL) {
      return 1;
   }

   if (wrdptr > mrkptr) {
      // Sum bit offsets.
      uint64_t i = mrkptr + bwt->n_symbols;
      while (i < wrdptr) 
         for (int j = 0; j < bwt->n_symbols; j++)
            offset[j] += __builtin_popcountl(bwt->occ[i++]);
      // Sum partial word.
      for (int j = 0; j < bwt->n_symbols; j++)
         offset[j] += __builtin_popcountl(bwt->occ[wrdptr++] >> (bwt->occ_word_size - 1 - bit));
      // Return sum.
      for (int j = 0; j < bwt->n_symbols; j++)
         occ[j] = bwt->occ[mrkptr++] + offset[j];
   } else {
      // Sum partial word.
      if (bit < bwt->occ_word_size - 1)
         for (int j = 0; j < bwt->n_symbols; j++)
            offset[j] += __builtin_popcountl(bwt->occ[wrdptr++] << (bit+1));
      else wrdptr += bwt->n_symbols;
      // Sum bit offsets.
      while (wrdptr < mrkptr)
         for (int j = 0; j < bwt->n_symbols; j++)
            offset[j] += __builtin_popcountl(bwt->occ[wrdptr++]);
      // Return subtraction.
      for (int j = 0; j < bwt->n_symbols; j++)
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
   
   if (sym < 0 || sym >= bwt->n_symbols)
      return -1;

   // Depth 0 pointer, return full span.
   if (ptr == -1) return bwt->occ[sym];
   // Compute word, bit and absolute marker positions.
   int64_t wrdnum = ptr / bwt->occ_word_size;
   int64_t wrdptr = (wrdnum + wrdnum/bwt->occ_mark_int + 1) * bwt->n_symbols + sym;
   int64_t mrkptr = (((wrdnum + bwt->occ_mark_int/2)/bwt->occ_mark_int) * (bwt->occ_mark_int+1))*bwt->n_symbols + sym;
   int64_t bit    = ptr % bwt->occ_word_size;

   uint64_t occ = bwt->occ[mrkptr];
   if (wrdptr > mrkptr) {
      int64_t  offset = 0;
      // Sum bit offsets.
      for (uint64_t i = mrkptr + bwt->n_symbols; i < wrdptr; i+= bwt->n_symbols)
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
      for (uint64_t i = wrdptr + bwt->n_symbols; i < mrkptr; i+= bwt->n_symbols)
         offset += __builtin_popcountl(bwt->occ[i]);
      // Return subtraction.
      occ -= offset;
   }
   return occ;
}
