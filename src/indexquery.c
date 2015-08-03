#include "indexquery.h"

uint64_t
get_sa
(
 uint64_t   pos,
 uint64_t * sa,
 int        bits
)
{
   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64-bits);
   uint64_t bit = pos*bits;
   uint64_t word = bit/64;
   bit %= 64;
   if (bit + bits > 64)
      return (((sa[word] >> bit) & mask) | (sa[word+1] & mask) << (64-bit)) & mask;
   else
      return (sa[word] >> bit) & mask;
}

int
get_occ
(
 int64_t   ptr,
 uint64_t * occ,
 int64_t * val
)
{
   if (ptr == -1) {
      for (int j = 0; j < NUM_BASES; j++) val[j] = 0;
      return 0;
   }
   int64_t wrdnum = ptr/OCC_WORD_SIZE;
   int64_t wrdptr = (wrdnum + wrdnum/OCC_MARK_INTERVAL + 1)*NUM_BASES;
   int64_t mrkptr = (((wrdnum + OCC_MARK_INTERVAL/2)/OCC_MARK_INTERVAL)*(OCC_MARK_INTERVAL+1))*NUM_BASES;
   int64_t bit    = ptr%OCC_WORD_SIZE;

   uint64_t * offset = calloc(NUM_BASES,sizeof(uint64_t));
   if (wrdptr > mrkptr) {
      // Sum bit offsets.
      uint64_t i = mrkptr + NUM_BASES;
      while (i < wrdptr) 
         for (int j = 0; j < NUM_BASES; j++)
            offset[j] += __builtin_popcountl(occ[i++]);
      // Sum partial word.
      for (int j = 0; j < NUM_BASES; j++)
         offset[j] += __builtin_popcountl(occ[wrdptr++] >> (OCC_WORD_SIZE - 1 - bit));
      // Returm sum.
      for (int j = 0; j < NUM_BASES; j++)
         val[j] = occ[mrkptr++] + offset[j];
   } else {
      // Sum partial word.
      if (bit < OCC_WORD_SIZE - 1)
         for (int j = 0; j < NUM_BASES; j++)
            offset[j] += __builtin_popcountl(occ[wrdptr++] << (bit+1));
      else wrdptr += NUM_BASES;
      // Sum bit offsets.
      while (wrdptr < mrkptr)
         for (int j = 0; j < NUM_BASES; j++)
            offset[j] += __builtin_popcountl(occ[wrdptr++]);
      // Return subtraction.
      for (int j = 0; j < NUM_BASES; j++)
         val[j] = occ[mrkptr++] - offset[j];
   }
   free(offset);
   return 0;
}


uint64_t
get_occ_nt
(
 int64_t   ptr,
 uint64_t * occ,
 int       nt
)
{
   if (ptr == -1) return 0;
   int64_t wrdnum = ptr/OCC_WORD_SIZE;
   int64_t wrdptr = (wrdnum + wrdnum/OCC_MARK_INTERVAL + 1)*NUM_BASES + nt;
   int64_t mrkptr = (((wrdnum + OCC_MARK_INTERVAL/2)/OCC_MARK_INTERVAL)*(OCC_MARK_INTERVAL+1))*NUM_BASES + nt;
   int64_t bit    = ptr%OCC_WORD_SIZE;

   int64_t offset = 0;
   if (wrdptr > mrkptr) {
      // Sum bit offsets.
      for (uint64_t i = mrkptr + NUM_BASES; i < wrdptr; i+= NUM_BASES)
            offset += __builtin_popcountl(occ[i]);
      // Sum partial word.
      offset += __builtin_popcountl(occ[wrdptr] >> (OCC_WORD_SIZE - 1 - bit));
      // Returm sum.
      return occ[mrkptr] + offset;
   } else {
      // Sum partial word.
      if (bit < OCC_WORD_SIZE - 1)
         offset += __builtin_popcountl(occ[wrdptr] << (bit+1));
      // Sum bit offsets.
      for (uint64_t i = wrdptr + NUM_BASES; i < mrkptr; i+= NUM_BASES)
         offset += __builtin_popcountl(occ[i]);
      // Return subtraction.
      return occ[mrkptr] - offset;
   }

}


fmdpos_t
extend_bw
(
 int        nt,
 fmdpos_t   pos,
 index_t  * index

)
{
   int64_t occ_sp[NUM_BASES];
   int64_t occ_ep[NUM_BASES];
   int64_t fp[NUM_BASES];
   int64_t rp[NUM_BASES];
   int64_t sz[NUM_BASES];
   get_occ(pos.fp - 1, index->occ, occ_sp);
   get_occ(pos.fp + pos.sz - 1, index->occ, occ_ep);
   for (int j = 0; j < NUM_BASES; j++) {
      fp[j] = index->c[j] + occ_sp[j];
      sz[j] = occ_ep[j] - occ_sp[j];
   }
   // Intervals of the reverse strand pointer.
   // T
   rp[3] = pos.rp;
   // G-C-A
   for (int j = 2; j >= 0; j--) rp[j] = rp[j+1] + sz[j+1];
   // N
   rp[4] = rp[0] + sz[0];
   
   return (fmdpos_t) {fp[nt], rp[nt], sz[nt]};
}


fmdpos_t
extend_fw
(
 int        nt,
 fmdpos_t   pos,
 index_t  * index
)
{
   // Reverse complement.
   if (nt < 4) nt = 3 - nt;
   // Lookup index using the reverse strand.
   fmdpos_t newpos = extend_bw(nt, (fmdpos_t) {pos.rp, pos.fp, pos.sz}, index);
   // Return position wrt forward strand.
   return (fmdpos_t){newpos.rp, newpos.fp, newpos.sz};
}

bwpos_t
bw_search
(
 int       nt,
 bwpos_t   pos,
 index_t * index
)
{
   bwpos_t newpos;
   newpos.sp = index->c[nt] + get_occ_nt(pos.sp - 1, index->occ, nt);
   newpos.ep = index->c[nt] + get_occ_nt(pos.ep, index->occ, nt) - 1;
   return newpos;
}


bwpos_t
bw_shrink
(
 bwpos_t    pos,
 index_t  * index
)
{
   return (bwpos_t){0,0};
}

