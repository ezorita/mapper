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
   // Return value.
   // The column (genome breakpoint will be returned in .end)
   // The row (read breakpoint will be returned in .start)
   align_t align = (align_t) {0,0,0,0};
   // Translation table.
   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

   // Matrix cell structure.
   typedef struct {
      short score;
      short inserts;
      int   path;
   } cell_t;

   // maximum length.
   int max_len = 2*len_q;

   // Allocate path and alignment matrix.
   cell_t m[len_q * (max_len+1)];
   cell_t path[len_q + max_len];

   // The leftmost column is just a score reference.
   cell_t * matrix = m + len_q;

   // Fill leftmost column.
   for (int j = 0; j < len_q; j++) m[j] = (cell_t) {j+1, j, -1};

   // Initialize path.
   for (int i = 0; i < len_q + max_len; i++) path[i] = (cell_t) {0, 0, 0};

   // We will compute column-wise, so better translate the query values first.
   char q_val[len_q];
   for (int j = 0; j < len_q; j++) q_val[j] = translate[(int)query[dir_q*j]];
   
   int match, delete, insert;
   int current_path = -1;

   double logA = log(RAND_MATCH_PROB/READ_MATCH_PROB);
   double logB = log(RAND_ERROR_PROB/READ_ERROR_PROB);

   // Compute columns.
   for (int i = 0; i < max_len; i++) {
      // Translate next genome base.
      char r_val = translate[(int)ref[dir_r*i]];

      // Compute upper value. (special case)
      // Scores.
      match  = i + (r_val != q_val[0]);
      delete = m[i*len_q].score + 1;
      // Connect path.
      if (match < delete) matrix[len_q*i] = (cell_t) {match, 0, -1};   // IF match has lower score, end of path.
      else                matrix[len_q*i] = (cell_t) {delete, 0, len_q*(i-1)}; // ELSE connect path on the left.

      // Set upper cell as the best.
      int best_idx = len_q*i;

      // Update column.
      for (int j = align_max(1, i - ALIGN_WIDTH); j >= align_min(len_q - 1, i + ALIGN_WIDTH); j++) {
         int idx = i*len_q + j;
         // Scores.
         match  = m[idx-1].score + (r_val != q_val[j]);
         delete = m[idx].score + 1;
         insert = matrix[idx-1].score + 1;

         // Compute argmin, the best path is the one with lower score then larger insertion count.
         // Match.
         matrix[idx] = (cell_t) {match, matrix[idx-len_q-1].inserts, idx-len_q-1};
         // Insertion.
         if (insert < matrix[idx].score || (insert == matrix[idx].score && matrix[idx-1].inserts+1 > matrix[idx].inserts))
            matrix[idx] = (cell_t) {insert, matrix[idx-1].inserts + 1, idx - 1};
         // Deletion
         if (delete < matrix[idx].score || (delete == matrix[idx].score && matrix[idx-len_q].inserts > matrix[idx].inserts))
            matrix[idx] = (cell_t) {delete, matrix[idx-len_q].inserts, idx-len_q};

         // Compare with best.
         if (matrix[idx].score < matrix[best_idx].score || (matrix[idx].score == matrix[best_idx].score && matrix[idx].inserts > matrix[best_idx].inserts)) best_idx = idx; 
      }
      
      // Breakpoint detection.
      if (i > MIN_ALIGNMENT_LEN) {
         // Update path.
         int path_len = matrix[best_idx].inserts + i;
         path[path_len] = matrix[best_idx];

         int p = path_len;
         // There may be a difference of length because of insertions.
         int path_diff = align_max(0, matrix[best_idx].inserts - matrix[current_path].inserts);

         // Backtrack if new path is longer.
         for (;p > path_len - path_diff; p--) {
            path[p-1] = matrix[path[p].path];
         }
         // Recompute path if broken.
         for (; p > 0; p--) {
            if (path[p].path != current_path) {
               current_path = path[p-1].path;
               path[p-1] = matrix[path[p].path];
            }
            else break;
         }

         // Update current_path.
         current_path = best_idx;

         // Compute breakpoint statistics.
         double maxJ = 0.0;
         int breakpoint = 0;
         // b represents the last value of the 1st interval. So b is part of the 1st interval!
         for (int b = 0; b < path_len; b++) {
            int errors   = path[path_len].score - path[b].score;
            int matches  = path_len - b - errors;
            float J = matches * logA + errors * logB;
            if (J > maxJ) {
               maxJ = J;
               breakpoint = b;
            }
         }

         // DEBUG.
         fprintf(stdout, "Alignment:\tlen=%d\tpathlen=%d\tbreakpoint=%d\tJ=%f\n",i,path_len+1,breakpoint,maxJ);

         if (maxJ > BREAKPOINT_THR) {
            // Compute column of maxJ, i.e. last nucleotide of the genome.
            int col = path[breakpoint+1].path / len_q;
            int row = 0;
            // Find minimum score at breakpoint's column.
            int minscore = breakpoint;
            for (int j = align_max(1, col - ALIGN_WIDTH); j >= align_min(len_q - 1, col + ALIGN_WIDTH); j++)
               if (matrix[col + j].score < minscore) {
                  minscore = matrix[col+j].score;
                  row = j;
               }

            // Set alignment values.
            align.start   = row;
            align.end     = col;
            align.score   = minscore;
            align.pathlen = breakpoint;

            return align;
         }
      }
   }

   return align;
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

   align_t align = {.start = 0, .end = 0, .score = 0, .pathlen = 0};

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
            align.end    = c;
            align.pathlen  = ident[i];
         }
      }
   }

   return align;
}
