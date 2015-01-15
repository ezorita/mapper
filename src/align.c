#include "align.h"

align_t
nw_align
(
 char * query,
 char * ref,
 int len_q,
 int dir_q,
 int dir_r
)
{
   struct cell_t {
      short score;
      short inserts;
      int   path;
   };

   // maximum length.
   int max_len = 2*len_q;

   // Allocate path and matrix.
   struct cell_t m[len_q * (max_len+1)];
   struct cell_t path[len_q + max_len];

   // The leftmost column is just a score reference.
   struct cell_t * matrix = m + len_q;

   // Fill leftmost column.
   for (int j = 0; j < len_q; j++) m[j] = j+1;

   // Initialize path.
   for (int i = 0; i < len_q + max_len; i++) path[i] = {0, 0, 0};

   // We will compute column-wise, so better translate the query values first.
   int q_val[len_q];
   for (int j = 0; j < len_q; j++) q_val[j] = translate[(int)query[j]];

   // Compute first column and assign cell (0,0) as best path.
   /*
   int r_val    = translate[(int)ref[0]];
   for (int j = 0; j < ALIGN_WIDTH; j++) {
      matrix[j].score = j + (r_val != q_val[j]);
      matrix[j].inserts = j;
      matrix[j].path = j-1;
   }
   path[0] = matrix[0];
   */
   
   int match, delete, insert;

   // Compute columns.
   for (int i = 0; i < max_len; i++) {
      // Translate next genome base.
      r_val    = translate[(int)ref[i]];

      // Compute upper value. (special case)
      // Scores.
      match  = i + (r_val != q_val[0]);
      delete = m[i*len_q] + 1;
      // Connect path.
      if (match < delete) matrix[len_q*i] = {match, 0, -1};           // IF match has lower score, end of path.
      else                matrix[len_q*i] = {delete, 0, len_q*(i-1)}; // ELSE connect path on the left.
      
      // Update column.
      for (int j = align_max(1, i - ALIGN_WIDTH); j >= align_min(len_q - 1, i + ALIGN_WIDTH); j++) {
         int idx = i*len_q + j;
         // Scores.
         match  = m[idx-1] + (r_val != q_val);
         delete = m[idx] + 1;
         insert = matrix[idx-1] + 1;
         
         // Assume match path.
         matrix[idx] = {match, matrix[idx-len_q-1], idx-len_q-1};
         // Correct if insertion is better.
         if (insert < matrix[idx].match || (insert == matrix[idx].match && 
         match, align_min(insert, delete));
         
         // Connect path.
         if      (score == delete) matrix[idx] = {score, matrix[idx-len_q].insert, idx-len_q};
         else if (score == insert) matrix[idx] = {score, matrix[idx-1].insert + 1, idx - 1};
         else    (score == match)  matrix[idx] = {score, matrix[idx-len_q-1].insert, idx-len_q-1};



         

         if 
            
         matrix[] 
      }
      
      // After the minimum alignment length, start computing the breakpoint scores.
      if (i > MIN_ALIGNMENT_LEN) {
          
      }
   }

}

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

   // Allocate score/start vectors.
   int score[rflen];
   int start[rflen];
   int ident[rflen];
   score[0] = 0;
   
   // 1st column.
   char read_val = translate[(int)read[0]];
   for (int i = 0; i < rflen; i++) {
      score[i] = ident[i] = (read_val == ref_val[i]);
      start[i] = 0;
   }

   // Next columns.
   for (int c = 1; c < rdlen; c++) {
      int match_old = 0;
      int id_old    = 0;
      int start_old = c;
      read_val = translate[(int)read[c]];
      for (int i = 0; i < rflen; i++) {
         // Match.
         int sc = align_max(0, match_old + SCORE_MATCH*(read_val == ref_val[i] ? 1 : -1));
         int st = (sc > 0 ? start_old : c);
         int id = (sc > 0 ? id_old : 0) + (read_val == ref_val[i]);
         // Deletion.
         if (i > 0 && score[i-1] > 1 && score[i-1] > sc) {
            if (sc == score[i-1] - 1) {
               if (start[i-1] > st) { // This will generate the shortest maximum.
                  st = start[i-1];
                  id = ident[i-1];
               }
            } else {
               sc = score[i-1] + (SCORE_DELETE);
               st = start[i-1];
               id = ident[i-1];
            }
         }

         // Insertion.
         if (score[i] > 1 &&  score[i] > sc) {
            if (sc == score[i] - 1) {
               if (start[i] > st) { // This will generate the shortest maximum
                  st = start[i];
                  id = ident[i];
               }
            } else {
               sc = score[i] + (SCORE_INSERT);
               st = start[i];
               id = ident[i];
            }
         }

         match_old = score[i];
         start_old = start[i];
         id_old    = ident[i];
         start[i] = st;
         score[i] = align_max(sc, 0);
         ident[i] = id;

         if (score[i] > align.score) {
            align.score  = score[i];
            align.start  = start[i];
            align.max    = c;
            align.ident  = ident[i];
         }
      }
   }

   return align;
}
