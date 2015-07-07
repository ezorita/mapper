#include "align.h"

void
dbf_update
(
 uint64_t   Eq,
 uint64_t * Pv,
 uint64_t * Mv,
 uint8_t  * Phin,
 uint8_t  * Mhin
 )
{
   uint64_t Xv, Xh, Ph, Mh;
   // Get current match sequence.
   Eq |= *Mhin;

   // Compute output of Xh blocks (Mh, Ph).
   Xh = (((Eq & *Pv) + *Pv) ^ *Pv) | Eq;
   Ph = *Mv | ~(Xh | *Pv);
   Mh = *Pv & Xh;

   // Save hout.
   uint8_t Phout = (Ph & MASK_MSB) > 0;
   uint8_t Mhout = (Mh & MASK_MSB) > 0;

   // Feed Ph, Mh to Xv blocks.
   Xv = Eq | *Mv;
   Ph = (Ph << 1) | *Phin;
   Mh = (Mh << 1) | *Mhin;
   *Pv = Mh | ~(Xv | Ph);
   *Mv = Ph & Xv;
   *Phin = Phout;
   *Mhin = Mhout;
}

void
dbf_fill_eq
(
 uint64_t ** Peq,
 char      * text,
 int         len,
 int         dir
)
{
char translate[256] = {[0 ... 255] = 4,
                       ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
                       ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3};

   uint32_t w = (len+WORD_SIZE-1)/WORD_SIZE;
   for (int i = 0; i < w; i++)
      for (int b = 0; b < WORD_SIZE && i*WORD_SIZE+b < len; b++) {
         int c = (int)translate[(uint8_t)text[dir*(i*WORD_SIZE+b)]];
         if (c < 4)  Peq[i][c] |= ((uint64_t)1) << b;
      }
}

void
dbf_free_eq
(
 uint64_t ** Eq,
 int         words
)
{
   for (int i = 0; i < words; i++) free(Eq[i]);
   free(Eq);
}

path_t
dbf_find_min
(
 uint32_t idx,
 uint32_t score,
 uint32_t awords,
 uint64_t * Pr,
 uint64_t * Mr,
 uint64_t * Pb,
 uint64_t * Mb
 )
{
   uint32_t cur = score;
   uint32_t min = score;
   uint32_t row = 0, col = 0;
   // Right wing.
   for (int j = awords-1, word = 0; j >= 0; j--, word++) {
      uint32_t abs_min = cur - __builtin_popcountl(Pr[j]);
      if (abs_min < min) {
         for (int i = WORD_SIZE-1, bit = 1; i >= 0; i--, bit++) {
            cur += ((Mr[j] >> i) & 1) - ((Pr[j] >> i) & 1);
            if (cur < min) {
               min = cur;
               row = WORD_SIZE*word + bit;
               col = 0;
            }
         } 
      } else {
         cur = abs_min + __builtin_popcountl(Mr[j]);
      }
   }
   // Bottom wing.
   cur = score;
   for (int j = awords-1, word = 0; j >= 0; j--, word++) {
      uint32_t abs_min = cur - __builtin_popcountl(Pb[j]);
      if (abs_min < min) {
         for (int i = WORD_SIZE-1, bit = 1; i >= 0; i--, bit++) {
            cur += ((Mb[j] >> i) & 1) - ((Pb[j] >> i) & 1);
            if (cur < min) {
               min = cur;
               col = WORD_SIZE*word + bit;
               row = 0;
            }
         } 
      } else {
         cur = abs_min + __builtin_popcountl(Mb[j]);
      }
   }

   return (path_t){.score = min, .row = idx - row, .col = idx - col};
}

path_t
dbf_align
(
 int         qlen,
 char      * qry,
 int         rlen,
 char      * ref,
 int         min_qlen,
 int         dir_q,
 int         dir_r,
 alignopt_t  opt
 )
{

char translate[256] = {[0 ... 255] = 4,
                       ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
                       ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3 };

   if (qlen < 1 || rlen < 1) return (path_t){0,0,0};
   int a_len = align_max(qlen,rlen);
   int align_w = (int)(a_len*opt.width_ratio);
   if (align_w < (qlen > rlen ? qlen - rlen : rlen - qlen)) return (path_t){-1,0,0};
   // Reduce width if align_w is bigger than the available query words.
   if (align_w > a_len) align_w = a_len;
   int awords = (align_w + WORD_SIZE - 1)/WORD_SIZE;
   align_w = awords * WORD_SIZE;

   // Number of Eq words.
   int eqwords = awords + (a_len + WORD_SIZE - 1)/WORD_SIZE;
   // Precompute query eq values.
   uint64_t ** Peq = malloc((eqwords+1)*sizeof(uint64_t*)); // Alloc 1 more to avoid segfault when shifting bits.
   for (int i = 0; i <= eqwords; i++) Peq[i] = calloc(4, sizeof(uint64_t));
   dbf_fill_eq(Peq + awords, qry, qlen, dir_q);
   // Precompute ref eq values.
   uint64_t ** Req = malloc((eqwords+1)*sizeof(uint64_t*));
   for (int i = 0; i <= eqwords; i++) Req[i] = calloc(4, sizeof(uint64_t));
   dbf_fill_eq(Req + awords, ref, rlen, dir_r);

   // Initialize
   // Right bitfield.
   uint64_t * Pr = calloc(awords, sizeof(uint64_t));
   uint64_t * Mr = malloc(awords * sizeof(uint64_t));
   memset(Mr, 0xFF, awords * sizeof(uint64_t));
   // Bottom bitfield.
   uint64_t * Pb = calloc(awords, sizeof(uint64_t));
   uint64_t * Mb = malloc(awords * sizeof(uint64_t));
   memset(Mb, 0xFF, awords * sizeof(uint64_t));

   // Alloc words for shifted Eq values.
   uint64_t * Eqr = malloc(awords * sizeof(uint64_t));
   uint64_t * Eqb = malloc(awords * sizeof(uint64_t));

   // Breakpoint.
   path_t   * bp_path = malloc((a_len/opt.bp_resolution+1)*sizeof(path_t));
   int    period  = 0;
   int    bp_pos  = 0;
   int    bp_ref  = 0;
   int    bp_cnt  = 0;
   double logAe = log(opt.rand_error) - log(opt.read_error);
   double logBe = log(1-opt.rand_error) - log(1-opt.read_error);

   // Horizontal/vertical input deltas.
   uint8_t Phinr, Mhinr, Phinb, Mhinb;
   uint32_t score = 0;

   for (int i = 0; i < a_len; i++) {
      Phinr = Phinb = 1;
      Mhinr = Mhinb = 0;

      // Update Eq.
      int b = i%WORD_SIZE;
      int w = i/WORD_SIZE;
      if (b) {
         for (int j = 0; j < awords; j++)
            Eqr[j] = (Peq[w+j][(int)translate[(uint8_t)ref[dir_r*i]]] >> b) | (Peq[w+j+1][(int)translate[(uint8_t)ref[dir_r*i]]] << (WORD_SIZE - b));
         for (int j = 0; j < awords; j++)
            Eqb[j] = (Req[w+j][(int)translate[(uint8_t)qry[dir_q*i]]] >> b) | (Req[w+j+1][(int)translate[(uint8_t)qry[dir_q*i]]] << (WORD_SIZE - b));
      }
      else {
         for (int j = 0; j < awords; j++) Eqr[j] = Peq[w+j][(int)translate[(uint8_t)ref[dir_r*i]]];
         for (int j = 0; j < awords; j++) Eqb[j] = Req[w+j][(int)translate[(uint8_t)qry[dir_q*i]]];
      }

      // Update bitfields.
      for (int j = 0; j < awords; j++) dbf_update(Eqr[j], Pr+j, Mr+j, &Phinr, &Mhinr);
      for (int j = 0; j < awords; j++) dbf_update(Eqb[j], Pb+j, Mb+j, &Phinb, &Mhinb);

      // Compute the central cell.
      uint64_t Pvc, Mvc, Phc, Mhc;
      int cref = translate[(uint8_t)ref[dir_r*i]];
      int cqry = translate[(uint8_t)qry[dir_q*i]];
      if (cref == cqry && cref < 4) {
         Mhc = Phinb;
         Phc = Mhinb;
         Mvc = Phinr;
         Pvc = Mhinr;
      } else {
         Mhc = Phinb & Mhinr;
         Phc = (!Phinb) & (Mhinb | (!Mhinr));
         Mvc = Phinr & Mhinb;
         Pvc = (!Phinr) & (Mhinr | (!Mhinb));
      }

      score += Phinr - Mhinr + Pvc - Mvc;

      // Now shift the bitfields and insert center cell.
      int j;
      for (j=0; j < awords-1; j++) Pr[j] = (Pr[j] >> 1) | (Pr[j+1] << (WORD_SIZE-1));
      Pr[j] = (Pr[j] >> 1) | (Pvc << (WORD_SIZE-1));
      for (j=0; j < awords-1; j++) Mr[j] = (Mr[j] >> 1) | (Mr[j+1] << (WORD_SIZE-1));
      Mr[j] = (Mr[j] >> 1) | (Mvc << (WORD_SIZE-1));
      for (j=0; j < awords-1; j++) Pb[j] = (Pb[j] >> 1) | (Pb[j+1] << (WORD_SIZE-1));
      Pb[j] = (Pb[j] >> 1) | (Phc << (WORD_SIZE-1));
      for (j=0; j < awords-1; j++) Mb[j] = (Mb[j] >> 1) | (Mb[j+1] << (WORD_SIZE-1));
      Mb[j] = (Mb[j] >> 1) | (Mhc << (WORD_SIZE-1));

      if (i%opt.bp_resolution == 0) {
         bp_path[bp_pos] = dbf_find_min(i,score,awords,Pr,Mr,Pb,Mb);
         if ((period++)%opt.bp_period == 1 || opt.bp_period == opt.bp_resolution) {
            // Compute likelihood of current reference.
            int Ee = bp_path[bp_pos].score - bp_path[bp_ref].score;
            int Me = i - bp_ref*opt.bp_resolution - Ee;
            float Je = Ee*logAe + Me*logBe;
            if (Je > opt.bp_thr) {
               bp_cnt++;
               //               fprintf(stdout, "\tML algorithm (i=%d, score=%d, cnt=%d) {bp_ref = %d (Je = %.2f)}\n", i, bp_path[bp_pos].score, bp_cnt, bp_ref*opt.bp_resolution, Je);
               if ((bp_cnt >= opt.bp_repeats && i >= min_qlen) || Je > opt.bp_max_thr) {
                  bp_cnt = opt.bp_repeats;
                  break;
               }
            } else {
               double maxJe;
               maxJe = -INFINITY;
               bp_cnt = 0;
               for (int b = 0; b < bp_pos; b += opt.bp_period) {
                  Ee = bp_path[bp_pos].score - bp_path[b].score;
                  Me = i - b*opt.bp_resolution - Ee;
                  Je = Ee*logAe + Me*logBe;
                  if (Je > maxJe) { maxJe = Je; bp_ref = b; }
               }
               if (maxJe > opt.bp_max_thr) {
                  bp_cnt = opt.bp_repeats;
                  break;
               }
               if (maxJe > opt.bp_thr) {
                  bp_cnt = 1;
               }
               //               fprintf(stdout, "\tML algorithm (i=%d, score=%d, cnt=%d) {*bp_ref = %d (Je = %.2f)}\n", i, bp_path[bp_pos].score, bp_cnt, bp_ref*opt.bp_resolution, maxJe);
            }
         }
         // Increase bp list position.
         bp_pos++;
      }
   }

   // High resolution breakpoint computation.
   path_t best_bp;
   int start = 0, end = bp_pos;
   if (bp_cnt >= opt.bp_repeats) {
      start = align_max(bp_ref - opt.bp_period + 1, 0);
      end   = align_min(bp_ref + opt.bp_period, bp_pos);
   }
   // Recompute highest likelihood.
   double maxJe = -INFINITY;
   //   fprintf(stdout, "FINE ML (%d to %d):\n", start, end);
   for (int b = start; b < end; b++) {
      int Ee = score - bp_path[b].score;
      int Me = bp_pos*opt.bp_resolution - b*opt.bp_resolution - Ee;
      double Je = Ee*logAe + Me*logBe;
      //      fprintf(stdout, "\tJe(%d)=%.2f [score = %d, len = %d]\n", b*opt.bp_resolution, Je, bp_path[b].score, a_len-1);
      if (Je > maxJe) { maxJe = Je; bp_ref = b; }
   }
   //   fprintf(stdout, "FINE ML (%d to %d) bp_ref=%d (maxJe = %.2f)\n", start, end, bp_ref*opt.bp_resolution, maxJe);
   if (maxJe >= opt.bp_thr) best_bp = bp_path[bp_ref];
   else best_bp = bp_path[bp_pos-1];

   // Free memory.
   free(Eqr);
   free(Eqb);
   dbf_free_eq(Peq, eqwords);
   dbf_free_eq(Req, eqwords);
   free(Mr); free(Pr); free(Mb); free(Pb);
   free(bp_path);

   if (best_bp.row >= qlen) best_bp.row = qlen-1;
   if (best_bp.col >= rlen) best_bp.row = rlen-1;


   return best_bp;
}
