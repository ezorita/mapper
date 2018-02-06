#include "bwtquery.h"

uint64_t
get_sa
(
 uint64_t   pos,
 sar_t    * sa
)
{
   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64 - sa->bits);
   uint64_t bit = pos*sa->bits;
   uint64_t word = bit/64;
   bit %= 64;
   if (bit + sa->bits > 64)
      return (((sa->sa[word] >> bit) & mask) | (sa->sa[word+1] & mask) << (64-bit)) & mask;
   else
      return (sa->sa[word] >> bit) & mask;
}

int
get_sa_range
(
 uint64_t   start,
 uint64_t   size,
 uint64_t * out,
 sar_t    * sa
)
{
   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64 - sa->bits);
   uint64_t bit = start*sa->bits;
   uint64_t word = bit/64;
   // Process bits.
   bit %= 64;
   uint64_t w = sa->sa[word];
   for (uint64_t i = 0; i < size; i++) {
      if (bit + sa->bits > 64) {
         out[i] = w >> bit;
         w = sa->sa[++word];
         out[i] |= (w & mask) << (64-bit);
         out[i] &= mask;
      } else {
         out[i] = (w >> bit) & mask;
      }
      bit = (bit + sa->bits)%64;
      if (bit == 0) w = sa->sa[++word];
   }
   return 0;
}

int
get_occ
(
 int64_t    ptr,
 uint64_t * occ,
 int64_t  * val
)
{
   if (ptr == -1) {
      for (int j = 0; j < NUM_BASES; j++) val[j] = occ[j];
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
      // Return sum.
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
 int64_t    ptr,
 uint64_t * occ,
 int        nt
)
{
   if (ptr == -1) return occ[nt];
   // Get word and bit pointers.
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
 bwt_t    * bwt
)
{
   if (pos.sz == 0) {
      pos.dp += 1;
      return pos;
   }
   int64_t occ_sp[NUM_BASES];
   int64_t occ_ep[NUM_BASES];
   int64_t fp[NUM_BASES];
   int64_t rp[NUM_BASES];
   int64_t sz[NUM_BASES];
   get_occ(pos.fp - 1, bwt->occ, occ_sp);
   get_occ(pos.fp + pos.sz - 1, bwt->occ, occ_ep);
   int64_t tot = 0;
   for (int j = 0; j < NUM_BASES; j++) {
      fp[j] = bwt->c[j] + occ_sp[j];
      sz[j] = occ_ep[j] - occ_sp[j];
      tot += sz[j];
   }
   // Intervals of the reverse strand pointer.
   // T (wildcard '$' is inside this interval if tot < pos.sz)
   rp[3] = pos.rp + (pos.sz - tot);
   // G-C-A
   for (int j = 2; j >= 0; j--) rp[j] = rp[j+1] + sz[j+1];
   // N
   rp[4] = rp[0] + sz[0];
   
   return (fmdpos_t) {fp[nt], rp[nt], sz[nt], pos.dp + 1};
}


fmdpos_t
extend_fw
(
 int        nt,
 fmdpos_t   pos,
 bwt_t    * bwt
)
{
   // Reverse complement.
   if (nt < 4) nt = 3 - nt;
   // Lookup index using the reverse strand.
   fmdpos_t newpos = extend_bw(nt, (fmdpos_t) {pos.rp, pos.fp, pos.sz, pos.dp}, bwt);
   // Return position wrt forward strand.
   return (fmdpos_t){newpos.rp, newpos.fp, newpos.sz, newpos.dp};
}

int
extend_bw_all
(
 fmdpos_t   pos,
 fmdpos_t * new,
 bwt_t    * bwt
)
{
   int64_t occ_sp[NUM_BASES];
   int64_t occ_ep[NUM_BASES];
   get_occ(pos.fp - 1, bwt->occ, occ_sp);
   get_occ(pos.fp + pos.sz - 1, bwt->occ, occ_ep);
   int64_t tot = 0;
   for (int j = 0; j < NUM_BASES; j++) {
      new[j].fp = bwt->c[j] + occ_sp[j];
      new[j].sz = occ_ep[j] - occ_sp[j];
      new[j].dp = pos.dp + 1;
      tot += new[j].sz;
   }
   // Intervals of the reverse strand pointer.
   // T (wildcard '$' is inside this interval if tot < pos.sz)
   new[3].rp = pos.rp + (pos.sz - tot);
   // G-C-A
   for (int j = 2; j >= 0; j--) new[j].rp = new[j+1].rp + new[j+1].sz;
   // N
   new[4].rp = new[0].rp + new[0].sz;
   
   return 1;
}

int
extend_fw_all
(
 fmdpos_t   pos,
 fmdpos_t * newpos,
 bwt_t    * bwt
)
{
   fmdpos_t tmp[NUM_BASES];
   extend_bw_all((fmdpos_t) {pos.rp, pos.fp, pos.sz, pos.dp}, tmp, bwt);
   newpos[0] = (fmdpos_t) {tmp[3].rp, tmp[3].fp, tmp[3].sz, tmp[3].dp};
   newpos[1] = (fmdpos_t) {tmp[2].rp, tmp[2].fp, tmp[2].sz, tmp[2].dp};
   newpos[2] = (fmdpos_t) {tmp[1].rp, tmp[1].fp, tmp[1].sz, tmp[1].dp};
   newpos[3] = (fmdpos_t) {tmp[0].rp, tmp[0].fp, tmp[0].sz, tmp[0].dp};
   newpos[4] = (fmdpos_t) {tmp[4].rp, tmp[4].fp, tmp[4].sz, tmp[4].dp};
   return 0;
}


int
suffix_string
(
 char    * suf,
 int       slen,
 uint64_t  minloci,
 bwpos_t * newpos,
 bwt_t   * bwt
)
// Extends from suf[0] until suf[slen-1] is reached or loci < minloci.
{
   int64_t loci;
   int len = slen-1;
   do {
      *newpos = bwt->bwt_base;
      for (int i = len; i >= 0; i--)
         if(suffix_extend(translate[(int)suf[i]], *newpos, newpos, bwt)) return -1;
      loci = (int64_t)newpos->ep - (int64_t)newpos->sp + 1;
      len--;
   } while (loci < minloci);

   return 0;
}

int
suffix_extend
(
 int       nt,
 bwpos_t   pos,
 bwpos_t * newpos,
 bwt_t   * bwt
)
{
   newpos->sp = bwt->c[nt] + get_occ_nt(pos.sp - 1, bwt->occ, nt);
   newpos->ep = bwt->c[nt] + get_occ_nt(pos.ep, bwt->occ, nt) - 1;
   newpos->depth = pos.depth + 1;
   return 0;
}

int
suffix_extend_all
(
 bwpos_t   pos,
 bwpos_t * newpos,
 bwt_t   * bwt
)
{
   // Compute occs for all nt in one call (1 cache miss).
   int64_t occ_sp[NUM_BASES];
   int64_t occ_ep[NUM_BASES];
   if (get_occ((int64_t)pos.sp - 1, bwt->occ, occ_sp)) return -1;
   if (get_occ((int64_t)pos.ep, bwt->occ, occ_ep)) return -1;
   // Update pointers.
   for (int nt = 0; nt < NUM_BASES; nt++) {
      newpos[nt].sp = bwt->c[nt] + occ_sp[nt];
      newpos[nt].ep = bwt->c[nt] + occ_ep[nt] - 1;
      newpos[nt].depth = pos.depth + 1;
   }
   return 0;
}

