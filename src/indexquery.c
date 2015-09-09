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
   if (ptr == -1) return occ[nt];
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

int
suffix_extend
(
 int       nt,
 bwpos_t   pos,
 bwpos_t * newpos,
 index_t * index
)
{
   newpos->sp = index->c[nt] + get_occ_nt(pos.sp - 1, index->occ, nt);
   newpos->ep = index->c[nt] + get_occ_nt(pos.ep, index->occ, nt) - 1;
   if (newpos->ep >= index->size || newpos->sp >= index->size) return -1;
   newpos->depth = pos.depth + 1;
   return 0;
}


int
suffix_shrink
(
 bwpos_t    pos,
 bwpos_t  * newpos,
 index_t  * index
)
{
   *newpos = pos;
   if (pos.sp == pos.ep)
      return suffix_ssv_search(pos.sp, newpos, index);
   else if (pos.sp < pos.ep)
      return suffix_ssv(pos, newpos, index);
   else return -1;
}

int
suffix_ssv_search
(
 uint64_t   pos,
 bwpos_t  * newpos,
 index_t  * index
)
// TODO:
// Possible algorithm update:
// Find only the first top corner when searching BWD and the first
// bottom corner when searching FWD, then add the offsets obtained
// from the sample list to reach the parent corner.
{
   *newpos = (bwpos_t){.depth = 0, .sp = pos, .ep = pos};
   uint64_t word  = pos/LCP_WORD_SIZE;
   uint64_t marks = word/LCP_MARK_INTERVAL + 1;
   int32_t  bit   = pos%LCP_WORD_SIZE;
   word += marks;

   uint64_t mark = marks * (LCP_MARK_INTERVAL + 1);

   uint64_t ptr = index->lcp_sample_idx[mark];
   for (uint64_t i = mark-1; i > word; i--)
      ptr -= __builtin_popcountl(index->lcp_sample_idx[i]);
   ptr -= __builtin_popcountl(index->lcp_sample_idx[word] >> bit);

   int is_top = 0, is_bot = 0;
   if ((index->lcp_sample_idx[word] >> bit) & 1) {
      // Pointer coincides with a corner.
      // Top or bottom?
      if (index->lcp_sample->lcp[ptr].offset < 0) {
         is_top = 1;
         // New depth.
         newpos->depth = index->lcp_sample->lcp[ptr+1].lcp;
      } else {
         is_bot = 1;
         // New depth.
         newpos->depth = index->lcp_sample->lcp[ptr].lcp;
      }
   } else {
      // Pointer is not a corner.
      // New depth.
      newpos->depth = index->lcp_sample->lcp[ptr].lcp;
   }

   // Forward search.
   if (!is_bot) {
      uint64_t topc = 1, samples = 0;
      uint64_t pptr = ptr + is_top, pword = word, pbit = bit;
      while (topc > 0) {
         lcpval_t val = index->lcp_sample->lcp[pptr++];
         topc += (val.offset < 0 ? 1 : -1);
         samples++;
      }
      // Now find the sample position in the BF.
      int32_t offset = 0;
      if (++pbit == LCP_WORD_SIZE) {
         pbit = 0;
         pword++;
      }
      uint64_t w = index->lcp_sample_idx[pword++] >> pbit;
      while (samples > 0) {
         int cnt = __builtin_popcountl(w);
         if (cnt >= samples) {
            while (samples > 0) {
               samples -= w & 1;
               w >>= 1;
               offset++;
            }
         } else {
            samples -= cnt;
            offset += LCP_WORD_SIZE - pbit;
            pbit = 0;
         }
         w = index->lcp_sample_idx[pword++];
      }
      newpos->ep = pos + offset;
   }

   // Backward search.
   if (!is_top) {
      uint64_t botc = 1, samples = 0;
      uint64_t nptr = ptr-1, nword = word, nbit = bit;
      while (botc > 0) {
         lcpval_t val = index->lcp_sample->lcp[nptr--];
         botc += (val.offset < 0 ? -1 : 1);
         samples++;
      }
      // Now find the sample position in the BF.
      int32_t offset = 0;
      if (--nbit == -1) {
         nbit = LCP_WORD_SIZE-1;
         nword--;
      }
      uint64_t w = index->lcp_sample_idx[nword--] << (LCP_WORD_SIZE - 1 - nbit);
      while (samples > 0) {
         int cnt = __builtin_popcountl(w);
         if (cnt >= samples) {
            while (samples > 0) {
               samples -= (w >> (LCP_WORD_SIZE-1)) & 1;
               w <<= 1;
               offset--;
            }
         } else {
            samples -= cnt;
            offset -= nbit;
            nbit = LCP_WORD_SIZE;
         }
         w = index->lcp_sample_idx[nword--];
      }
      newpos->sp = pos + offset;
   }

   return 0;
}

int
suffix_ssv
// Find sampled smaller value.
(
 bwpos_t    pos,
 bwpos_t  * newpos,
 index_t  * index
)
{
   *newpos = pos;
   
   // Start LCP.
   uint64_t sword  = pos.sp/LCP_WORD_SIZE;
   uint64_t smarks = sword/LCP_MARK_INTERVAL + 1;
   int32_t  sbit   = pos.sp%LCP_WORD_SIZE;
   sword += smarks;
   uint64_t smark = smarks*(LCP_MARK_INTERVAL+1);
   if (((index->lcp_sample_idx[sword] >> sbit) & 1) == 0) return -1;

   uint64_t ptr = index->lcp_sample_idx[smark];
   for (uint64_t i = smark-1; i > sword; i--)
      ptr -= __builtin_popcountl(index->lcp_sample_idx[i]);
   ptr -= __builtin_popcountl(index->lcp_sample_idx[sword] >> sbit);

   uint8_t slcp = index->lcp_sample->lcp[ptr].lcp;
   int32_t soff = index->lcp_sample->lcp[ptr].offset;
   if (soff >= 0) return -1;

   // End LCP.
   uint64_t eword  = pos.ep/LCP_WORD_SIZE;
   uint64_t emarks = eword/LCP_MARK_INTERVAL + 1;
   int32_t  ebit   = pos.ep%LCP_WORD_SIZE;
   eword += emarks;
   uint64_t emark = emarks*(LCP_MARK_INTERVAL+1);
   if (((index->lcp_sample_idx[eword] >> ebit) & 1) == 0) return -1;

   ptr = index->lcp_sample_idx[emark];
   for (uint64_t i = emark-1; i > eword; i--)
      ptr -= __builtin_popcountl(index->lcp_sample_idx[i]);
   ptr -= __builtin_popcountl(index->lcp_sample_idx[eword] >> ebit);

   uint8_t elcp = index->lcp_sample->lcp[ptr+1].lcp;
   int32_t eoff = index->lcp_sample->lcp[ptr].offset;
   if (eoff < 0) return -1;

   // Compare LCP.
   if (slcp >= elcp) {
      if (soff == -128) {
         ptr = index->lcp_extend_idx[smark];
         for (uint64_t i = smark-1; i > sword; i--)
            ptr -= __builtin_popcountl(index->lcp_extend_idx[i]);
         ptr -= __builtin_popcountl(index->lcp_extend_idx[sword] >> sbit);

         soff = index->lcp_extend->val[ptr];
      }
      newpos->sp = pos.sp + soff;
      newpos->depth = slcp;
   }
   if (slcp <= elcp) {
      if (++eoff == 128) {
         ptr = index->lcp_extend_idx[emark];
         for (uint64_t i = emark-1; i > eword; i--)
            ptr -= __builtin_popcountl(index->lcp_extend_idx[i]);
         ptr -= __builtin_popcountl(index->lcp_extend_idx[eword] >> ebit);

         eoff = index->lcp_extend->val[ptr];
      }
      newpos->ep = pos.ep + eoff;
      newpos->depth = elcp;
   }
   return 0;
}
