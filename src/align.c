#include "align.h"

align_t
sw_align
(
 char * read,
 int    rdlen,
 char * ref,
 int    rflen
)
{

   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

   // TODO:
   // Add compatibility for GAP_OPEN/GAP_EXTEND scores.

   align_t align = {.start = 0, .max = 0, .score = 0};

   char ref_val[rflen];
   for (int i = 0; i < rflen; i++) ref_val[i] = translate[(int)ref[i]];

   if (rdlen < 1 || rflen < 1) return align;
   int score[rflen];
   int start[rflen];
   score[0] = 0;
   
   // 1st column.
   char read_val = translate[(int)read[0]];
   for (int i = 0; i < rflen; i++) {
      score[i] = (read_val == ref_val[i]);
      start[i]   = 0;
   }

   // Next columns.
   for (int c = 1; c < rdlen; c++) {
      int match_old = 0;
      int start_old = c;
      read_val = translate[(int)read[c]];
      for (int i = 0; i < rflen; i++) {
         // Match.
         int sc = align_max(0, match_old + SCORE_MATCH*(read_val == ref_val[i] ? 1 : -1));
         int st = (sc > 0 ? start_old : c);
         // Deletion.
         if (i > 0 && score[i-1] > 1 && score[i-1] > sc) {
            if (sc == score[i-1] - 1) {
               st = align_min(start[i-1], st);
            } else {
               sc = score[i-1] + SCORE_DELETE;
               st = start[i-1];
            }
         }
         // Insertion.
         if (score[i] > 1 &&  score[i] > sc) {
            if (sc == score[i] - 1) {
               st = align_min(start[i] , st);
            } else {
               sc = score[i] + SCORE_INSERT;
               st = start[i];
            }
         }

         match_old = score[i];
         start_old = start[i];
         start[i] = st;
         score[i] = align_max(sc, 0);

         if (score[i] > align.score) {
            align.score  = score[i];
            align.start  = start[i];
            align.max    = c;
         }
      }
   }

   return align;
}
