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
 uint64_t  ** Peq,
 const char * text,
 const int    len,
 const int    dir
)
{
static const char translate[256] = {[0 ... 255] = 4,
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
align_bp
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
   int w = (int)(align_max(rlen,qlen)*opt.width_ratio);
   if (w > MAX_NAIVE_WIDTH) {
      return dbf_align_bp(qlen, qry, rlen, ref, min_qlen, dir_q, dir_r, opt);
   } else {
      return naive_align_bp(qlen, qry, rlen, ref, min_qlen, dir_q, dir_r, opt);
   }
}

path_t
dbf_align_bp
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
   // TODO:
   // To align until align_max(qlen,rlen) the ends of qry must be padded with 'N'.
   int a_len = align_min(qlen,rlen);
   int align_w = (int)(a_len*opt.width_ratio);
   if (align_w < (qlen > rlen ? qlen - rlen : rlen - qlen)) return (path_t){-1,0,0};
   // Reduce width if align_w is bigger than the available query words.
   if (align_w > a_len) align_w = a_len;
   int awords = align_max(1,(align_w + WORD_SIZE - 1)/WORD_SIZE);
   align_w = awords * WORD_SIZE;

   // Number of Eq words.
   int eqwords = awords + (a_len + WORD_SIZE - 1)/WORD_SIZE;
   // Precompute query eq values.
   uint64_t ** Peq = malloc((eqwords+1)*sizeof(uint64_t*)); // Alloc 1 more to avoid segfault when shifting bits.
   for (int i = 0; i <= eqwords; i++) Peq[i] = calloc(5, sizeof(uint64_t)); // 5 bases (A,C,G,T,N)
   dbf_fill_eq(Peq + awords, qry, qlen, dir_q);
   // Precompute ref eq values.
   uint64_t ** Req = malloc((eqwords+1)*sizeof(uint64_t*));
   for (int i = 0; i <= eqwords; i++) Req[i] = calloc(5, sizeof(uint64_t)); // 5 bases (A,C,G,T,N)
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
   double logAe = opt.logAe;//log(opt.rand_error) - log(opt.read_error);
   double logBe = opt.logBe;//log(1-opt.rand_error) - log(1-opt.read_error);

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
      int j = 0;
      for (j=0; j < awords-1; j++) Pr[j] = (Pr[j] >> 1) | (Pr[j+1] << (WORD_SIZE-1));
      Pr[j] = (Pr[j] >> 1) | (Pvc << (WORD_SIZE-1));
      for (j=0; j < awords-1; j++) Mr[j] = (Mr[j] >> 1) | (Mr[j+1] << (WORD_SIZE-1));
      Mr[j] = (Mr[j] >> 1) | (Mvc << (WORD_SIZE-1));
      for (j=0; j < awords-1; j++) Pb[j] = (Pb[j] >> 1) | (Pb[j+1] << (WORD_SIZE-1));
      Pb[j] = (Pb[j] >> 1) | (Phc << (WORD_SIZE-1));
      for (j=0; j < awords-1; j++) Mb[j] = (Mb[j] >> 1) | (Mb[j+1] << (WORD_SIZE-1));
      Mb[j] = (Mb[j] >> 1) | (Mhc << (WORD_SIZE-1));

      if (i%opt.bp_resolution == 0) {
         if (opt.bp_diagonal) {
            bp_path[bp_pos] = (path_t) {score, i, i};
         } else {
            bp_path[bp_pos] = dbf_find_min(i,score,awords,Pr,Mr,Pb,Mb);
         }
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
   dbf_free_eq(Peq, eqwords+1);
   dbf_free_eq(Req, eqwords+1);
   free(Mr); free(Pr); free(Mb); free(Pb);
   free(bp_path);

   if (best_bp.row >= qlen) best_bp.row = qlen-1;
   if (best_bp.col >= rlen) best_bp.row = rlen-1;


   return best_bp;
}


path_t
naive_align_bp
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

   char t[256] = {[0 ... 255] = 4,
                       ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
                       ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3 };

   if (qlen < 1 || rlen < 1) return (path_t){0,0,0};
   // TODO:
   // To align until align_max(qlen,rlen) the ends of qry must be padded with 'N'.
   int a_len = align_min(qlen,rlen);
   int w = align_max(opt.width_min, (int)(a_len*opt.width_ratio));
   if (w < (qlen > rlen ? qlen - rlen : rlen - qlen)) return (path_t){-1,0,0};
   // Reduce width if align_w is bigger than the available query words.
   if (w > a_len) w = a_len;

   // Breakpoint.
   path_t   * bp_path = malloc((a_len/opt.bp_resolution+1)*sizeof(path_t));
   int    period  = 0;
   int    bp_pos  = 0;
   int    bp_ref  = 0;
   int    bp_cnt  = 0;
   double logAe = opt.logAe;//log(opt.rand_error) - log(opt.read_error);
   double logBe = opt.logBe;//log(1-opt.rand_error) - log(1-opt.read_error);

   // Alloc cells.
   int * cells = malloc((2*w+3)*sizeof(int));
   int * c = cells + w + 1;
   // Initialize scores.
   c[-w-1] = c[w+1] = a_len;
   for (int i = 1; i <= w; i++) {
      c[-i] = c[i] = i;
   }
   c[0] = 0;

   int score, col = 0, row = 0;
   for (int i = 0; i < a_len; i++) {
      score = a_len;
      // Horizontal.
      for (int32_t j = (i < w ? -i : -w); j < 0; j++) {
         c[j] = align_min(align_min(c[j+1], c[j-1]) + 1, c[j] + (t[(int)qry[dir_q*(i+j)]] != t[(int)ref[dir_r*i]]));
         if (c[j] <= score) {
            col = i+j;
            row = i;
            score = c[j];
         }
      }
      // Vertical and center.
      for (int32_t j = (i < w ? i : w); j >= 0; j--) {
         c[j] = align_min(align_min(c[j+1], c[j-1]) + 1, c[j] + (t[(int)qry[dir_q*i]] != t[(int)ref[dir_r*(i-j)]]));
         if (c[j] <= score) {
            col = i;
            row = i-j;
            score = c[j];
         }
      }

      // Check alignment width.
      if (col-row == w || row-col == w) {
         opt.width_min = w+2;
         return naive_align_bp(qlen,qry,rlen,ref,min_qlen,dir_q,dir_r, opt);
      }
      
      if (i%opt.bp_resolution == 0) {
         bp_path[bp_pos] = (path_t) {score, row, col};
         if ((period++)%opt.bp_period == 1 || opt.bp_period == opt.bp_resolution) {
            // Compute likelihood of current reference.
            int Ee = bp_path[bp_pos].score - bp_path[bp_ref].score;
            int Me = i - bp_ref*opt.bp_resolution - Ee;
            float Je = Ee*logAe + Me*logBe;
            if (BREAKPOINT_DEBUG) {
               fprintf(stdout, "\tML algorithm (i=%d, score=%d, cnt=%d) {bp_ref = %d (Je = Ee(%d)*logAe(%.2f)+Me(%d)*logBe(%.2f) = %.2f)}\n", i, bp_path[bp_pos].score, bp_cnt, bp_ref*opt.bp_resolution, Ee, logAe, Me, logBe, Je);
            }
            if (Je > opt.bp_thr) {
               bp_cnt++;
               if (BREAKPOINT_DEBUG) {
                  fprintf(stdout, "\tML algorithm (i=%d, score=%d, cnt=%d) {bp_ref = %d (Je = %.2f)}\n", i, bp_path[bp_pos].score, bp_cnt, bp_ref*opt.bp_resolution, Je);
               }
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
               if (BREAKPOINT_DEBUG) {
                  fprintf(stdout, "\tML algorithm (i=%d, score=%d, cnt=%d) {*bp_ref = %d (Je = %.2f)}\n", i, bp_path[bp_pos].score, bp_cnt, bp_ref*opt.bp_resolution, maxJe);
               }
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
      end   = align_min(bp_ref + opt.bp_period + 1, bp_pos);
   }
   // Recompute highest likelihood.
   double maxJe = -INFINITY;
   if (BREAKPOINT_DEBUG) {
      fprintf(stdout, "FINE ML (%d to %d):\n", start, end-1);
   }
   for (int b = start; b < end; b++) {
      int Ee = score - bp_path[b].score;
      int Me = bp_pos*opt.bp_resolution - b*opt.bp_resolution - Ee;
      double Je = Ee*logAe + Me*logBe;
      if (BREAKPOINT_DEBUG) {
         fprintf(stdout, "\tJe(%d)=%.2f [score = %d, len = %d]\n", b*opt.bp_resolution, Je, bp_path[b].score, a_len);
      }
      if (Je > maxJe) { maxJe = Je; bp_ref = b; }
   }
   if (BREAKPOINT_DEBUG) {
      fprintf(stdout, "FINE ML (%d to %d) bp_ref=%d (maxJe = %.2f)\n", start, end, bp_ref*opt.bp_resolution, maxJe);
   }
   if (maxJe >= opt.bp_thr) best_bp = bp_path[bp_ref];
   else best_bp = bp_path[bp_pos-1];

   // Free memory.
   free(cells);
   free(bp_path);

   if (best_bp.row >= qlen) best_bp.row = qlen-1;
   if (best_bp.col >= rlen) best_bp.row = rlen-1;

   return best_bp;
}

int
align
(
 int         qlen,
 char      * qry,
 int         rlen,
 char      * ref,
 int         dir_q,
 int         dir_r,
 int         max_score,
 alignopt_t  opt,
 path_t    * path
 )
{
   int w = (int)(align_max(rlen,qlen)*opt.width_ratio);
   if (w > MAX_NAIVE_WIDTH) {
      return dbf_align(qlen, qry, rlen, ref, dir_q, dir_r, max_score, opt, path);
   } else {
      return naive_align(qlen, qry, rlen, ref, dir_q, dir_r, max_score, opt, path);
   }
}


int
dbf_align
(
 int         qlen,
 char      * qry,
 int         rlen,
 char      * ref,
 int         dir_q,
 int         dir_r,
 int         max_score,
 alignopt_t  opt,
 path_t    * path
 )
{

   char translate[256] = {[0 ... 255] = 4,
                          ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
                          ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3 };

   if (qlen < 1 || rlen < 1) {
      *path = (path_t){0,0,0};
      return 0;
   }
   int a_len = align_max(qlen,rlen);
   int align_w = align_max(opt.width_min,(int)(a_len*opt.width_ratio));
   if (align_w < (qlen > rlen ? qlen - rlen : rlen - qlen)) return -1;
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

   // Horizontal/vertical input deltas.
   uint8_t Phinr, Mhinr, Phinb, Mhinb;
   uint32_t score = 0;
   int i = 0;
   for (; i < a_len; i++) {
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
      if (score > max_score) {
         return 1;
      }

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
   }

   *path = dbf_find_min(i-1,score,awords,Pr,Mr,Pb,Mb);
   // TODO: this will not detect partial indels that compensate back to 0.

   // Free memory.
   free(Eqr);
   free(Eqb);
   dbf_free_eq(Peq, eqwords);
   dbf_free_eq(Req, eqwords);
   free(Mr); free(Pr); free(Mb); free(Pb);

   if (path->row >= qlen) path->row = qlen-1;
   if (path->col >= rlen) path->col = rlen-1;


   return 0;
}

int
naive_align
(
 int         qlen,
 char      * qry,
 int         rlen,
 char      * ref,
 int         dir_q,
 int         dir_r,
 int         max_score,
 alignopt_t  opt,
 path_t    * path
 )
{

   char t[256] = {[0 ... 255] = 4,
                          ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
                          ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3 };

   if (qlen < 1 || rlen < 1) {
      *path = (path_t){0,0,0};
      return 0;
   }
   int a_len = align_max(qlen,rlen);
   int w = align_max(opt.width_min,(int)(a_len*opt.width_ratio));
   if (w < (qlen > rlen ? qlen - rlen : rlen - qlen)) return -1;
   // Reduce width if align_w is bigger than the available query words.
   if (w > a_len) w = a_len;

   // Alloc cells.
   int * cells = malloc((2*w+3)*sizeof(int));
   int * c = cells + w + 1;
   // Initialize scores.
   c[-w-1] = c[w+1] = a_len;
   for (int i = 1; i <= w; i++) {
      c[-i] = c[i] = i;
   }
   c[0] = 0;
   // Iterate.
   int score = a_len, col = 0, row = 0;
   for (int i = 0; i < a_len; i++) {
      score = a_len;
      // Horizontal.
      for (int32_t j = (i < w ? -i : -w); j < 0; j++) {
         c[j] = align_min(align_min(c[j+1], c[j-1]) + 1, c[j] + (t[(int)qry[dir_q*(i+j)]] != t[(int)ref[dir_r*i]]));
         if (c[j] <= score) {
            score = c[j];
            col = i+j;
            row = i;
         }
      }
      // Vertical and center.
      for (int32_t j = (i < w ? i : w); j >= 0; j--) {
         c[j] = align_min(align_min(c[j+1], c[j-1]) + 1, c[j] + (t[(int)qry[dir_q*i]] != t[(int)ref[dir_r*(i-j)]]));
         if (c[j] <= score) {
            score = c[j];
            col = i;
            row = i-j;
         }
      }
      if (score > max_score) return 1;
   }

   // Save path.
   *path = (path_t){score, row, col};
   if (path->row >= qlen) path->row = qlen-1;
   if (path->col >= rlen) path->col = rlen-1;

   return 0;
}



int
mapq_align
(
 const char * qry,
 const char * qlt,
 const char * ref,
 const int    len,
 const int    w,
       int  * dist
)
{

   static const char t[256] = {[0 ... 255] = 4,
                  ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3,
                  ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3 };
   // Mismatch penalty
   int m_p[256];
   for (int i = 0; i < 256; i++) m_p[i] = i-33;
   // Indel penalty is maxQ.
   int i_p = 0;
   for (int i = 0; i < len; i++) if (i_p < qlt[i]) i_p = qlt[i];
   //   i_p = (m_p[i_p]+1)/2;
   i_p = m_p[i_p];
   
   // Alloc cells.
   int * cells = malloc((2*w+3)*sizeof(int));
   int * mism  = malloc((2*w+3)*sizeof(int));
   int * m = mism + w + 1;
   int * c = cells + w + 1;
   // Initialize scores.
   c[0] = m[0] = 0;
   c[-w-1] = c[w+1] = i_p * len;
   m[-w-1] = m[w+1] = len;
   for (int i = 1; i <= w; i++) {
      c[-i] = c[i] = i * i_p;
      m[-i] = m[i] = i;
   }

   for (int32_t i = 0; i < len; i++) {
      // Horizontal.
      for (int32_t j = (i < w ? -i : -w); j < 0; j++) {
         int indel = align_min(c[j+1], c[j-1]) + i_p;
         int match = (t[(int)qry[i+j]] == t[(int)ref[i]] ? 0 : m_p[(int)qlt[i+j]]);
         c[j] = align_min(indel, c[j] + match);
         m[j] = align_min(align_min(m[j+1], m[j-1]) + 1, m[j] + (match > 0));
      }
      // Vertical and center.
      double mp = m_p[(int)qlt[i]];
      for (int32_t j = (i < w ? i : w); j >= 0; j--) {
         int indel = align_min(c[j+1], c[j-1]) + i_p;
         int match = (t[(int)qry[i]] == t[(int)ref[i-j]] ? 0 : mp);
         c[j] = align_min(indel, c[j] + match);
         m[j] = align_min(align_min(m[j+1], m[j-1]) + 1, m[j] + (match > 0));
      }
   }

   // Find minimum and return.
   int mapq = c[w];
   *dist = m[w];
   for (int i = -w; i < w; i++) {
      if (c[i] < mapq) mapq = c[i], *dist = m[i];
   }

   // Free memory.
   free(cells);
   free(mism);

   return mapq;
}
