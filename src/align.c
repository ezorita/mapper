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
   typedef struct cell_t cell_t;
   struct cell_t {
      short    score;
      short    r_gaps;
      short    row;
      short    col;
      struct cell_t * path;
   };

   // Allocate path and alignment matrix.
   int VEC_SIZE = 2*ALIGN_WIDTH + 3;
   int allocated = align_min(len_q+ALIGN_WIDTH+1, ALLOC_BLOCK_SIZE);
   cell_t ** Ls = malloc((len_q+ALIGN_WIDTH+1)*sizeof(cell_t*));

   for (int i = 0; i < allocated; i++) {
      Ls[i] = malloc(VEC_SIZE * sizeof(cell_t));
   }
   cell_t ** path = malloc((2*len_q+ALIGN_WIDTH)*sizeof(cell_t *));

   // Translate sequences.
   char val_q[len_q];
   char val_r[len_q];

   for (int i = 0; i < len_q; i++) {
      val_q[i] = translate[(int)query[dir_q*i]];
      val_r[i] = translate[(int)ref[dir_r*i]];
   }

   // Aux vars.
   int match, r_gap, g_gap;

   double logA = log(READ_ERROR_PROB);
   double logB = log(READ_MATCH_PROB);
   double logC = log(RAND_ERROR_PROB);
   double logD = log(RAND_MATCH_PROB);

   // Initialize (0,0) value.
   cell_t * uL = Ls[0] + ALIGN_WIDTH + 1;
   uL[0] = (cell_t) {0, 0, 0};

   // Initialize path.
   for (int i = 0; i < 2*len_q; i++) path[i] = NULL;

   // Compute inverted Ls.
   for (int i = 0; i < len_q + ALIGN_WIDTH; i++) {
      int l = i + 1;
      // Expand matrix.
      if (l >= allocated) {
         allocated = align_min(len_q+ALIGN_WIDTH+1, allocated + ALLOC_BLOCK_SIZE);
         for (int k = l; k < allocated; k++) {
            Ls[k] = malloc(VEC_SIZE * sizeof(cell_t));
         }
      }
      cell_t * L = Ls[l] + ALIGN_WIDTH + 1;

      int width;
      // Matrix border elements.
      if (i <= ALIGN_WIDTH) {
         width = i;
         // Vertical L border.
         if (i < len_q) {
            L[width+1] = (cell_t) {uL[width].score + 1, uL[width].r_gaps + 1, i - width - 1, i, uL + width};
         }
         // Horizontal L border.
         L[-width-1] = (cell_t) {uL[-width].score + 1, uL[-width].r_gaps, i, i - width - 1, uL - width};
      } else {
         width = ALIGN_WIDTH;
         // Vertical L border.
         if (i < len_q) {
            r_gap = uL[width].score + 1;
            match = uL[width+1].score + (val_q[i-width-1] != val_r[i]);
            if (match < r_gap || (match == r_gap && uL[width+1].r_gaps > uL[width].r_gaps))
               L[width+1] = (cell_t) {match, uL[width+1].r_gaps, i - width - 1, i, uL + width + 1};
            else
               L[width+1] = (cell_t) {r_gap, uL[width].r_gaps + 1, i - width - 1, i, uL + width};
         }
         // Horizontal L border.
         g_gap = uL[-width].score + 1;
         match = uL[-width-1].score + (val_q[i] != val_r[i-width-1]);
         if (match < g_gap || (match == g_gap && uL[-width-1].r_gaps >= uL[-width].r_gaps))
            L[-width-1] = (cell_t) {match, uL[-width-1].r_gaps, i, i - width - 1, uL - width - 1};
         else
            L[-width-1] = (cell_t) {g_gap, uL[-width].r_gaps, i, i - width - 1, uL - width};

      }

      // Compute L elements.
      for (int j = width ; j > align_max(0, i - len_q); j--) {
         // Vertical wing.
         if (i < len_q) {
            match  = uL[j].score + (val_q[i-j] != val_r[i]);
            r_gap = uL[j-1].score + 1;
            g_gap = L[j+1].score + 1;

            // Assume r_gap.
            L[j] = (cell_t) {r_gap, uL[j-1].r_gaps + 1, i - j, i, uL + j - 1};
            // Check g_gap.
            if (g_gap < L[j].score || (g_gap == L[j].score && L[j+1].r_gaps > L[j].r_gaps))
               L[j] = (cell_t) {g_gap, L[j+1].r_gaps, i - j, i, L + j + 1};
            // Check match.
            if (match < L[j].score || (match == L[j].score && uL[j].r_gaps >= L[j].r_gaps))
               L[j] = (cell_t) {match, uL[j].r_gaps, i - j, i, uL + j};
         }

         // Horizontal wing.
         match = uL[-j].score + (val_q[i] != val_r[i-j]);
         r_gap = L[-j-1].score + 1;
         g_gap = uL[-j+1].score + 1;

         // Assume r_gap.
         L[-j] = (cell_t) {r_gap, L[-j-1].r_gaps + 1, i, i - j, L - j - 1};
         // Check g_gap.
         if (g_gap < L[-j].score || (g_gap == L[-j].score && uL[-j+1].r_gaps > L[-j].r_gaps))
            L[-j] = (cell_t) {g_gap, uL[-j+1].r_gaps, i, i - j, uL - j + 1};
         // Check match.
         if (match < L[-j].score || (match == L[-j].score && uL[-j].r_gaps >= L[-j].r_gaps))
            L[-j] = (cell_t) {match, uL[-j].r_gaps, i, i - j, uL - j};
      }

      // Center cell.
      if (i < len_q) {
         match = uL[0].score + (val_q[i] != val_r[i]);
         r_gap = L[-1].score + 1;
         g_gap = L[1].score + 1;

         // Assume r_gap.
         L[0] = (cell_t) {r_gap, L[-1].r_gaps + 1, i, i, L - 1};
         // Check g_gap.
         if (g_gap < L[0].score || (g_gap == L[0].score && L[1].r_gaps > L[0].r_gaps))
            L[0] = (cell_t) {g_gap, L[1].r_gaps, i, i, L + 1};
         // Check match.
         if (match < L[0].score || (match == L[0].score && uL[0].r_gaps >= L[0].r_gaps))
            L[0] = (cell_t) {match, uL[0].r_gaps, i, i, uL};
      }

      // Find highest identity.
      // Priorities: 1. Center value. 2. Vertical wing. 3. Horizontal wing.
      double minerr = 1.0;
      int best_idx = 0;
      for (int j = -width; j < align_min(0, len_q - i); j++) {
         double err = L[j].score*1.0/(L[j].r_gaps + L[j].row + 1);
         if (err <= minerr) {
            minerr = err;
            best_idx = j;
         }
      }
      if (i < len_q) {
         for (int j = width; j >= 0; j--) {
            double err = L[j].score*1.0/(L[j].r_gaps + L[j].row + 1);
            if (err <= minerr) {
               minerr = err;
               best_idx = j;
            }
         }
      }

      // Update path.
      int path_len = L[best_idx].r_gaps + L[best_idx].row;
      path[path_len] = L + best_idx;

      // Recompute path if broken.
      for (int p = path_len; p > 0; p--) {
         if (path[p-1] == NULL || path[p]->path != path[p-1]) {
            // Recover vector position from cellid.
            path[p-1] = path[p]->path;
         }
         else break;
      }

      double maxJ;
     
      // Breakpoint detection.
      if (i >= MIN_ALIGNMENT_LEN) {

         // Compute breakpoint statistics.
         maxJ = -INFINITY;
         int breakpoint = 0;
         // b represents the last value of the 1st interval. So b is part of the 1st interval!
         for (int b = 0; b < path_len; b++) {
            int E1 = path[path_len]->score - path[b]->score;
            int M1 = path_len - b - E1;
            float J = E1*(logC-logA) + M1*(logD-logB);
            if (J > maxJ) {
               maxJ = J;
               breakpoint = b;
            }
         }
         
         if (maxJ > BREAKPOINT_THR) {
            cell_t * bp = path[breakpoint];
            // Compute column of maxJ, i.e. last nucleotide of the genome.
            // Set alignment values.
            align.start   = align_max(0,bp->row);
            align.end     = align_max(0,bp->col);
            align.score   = bp->score;
            align.pathlen = breakpoint + 1;
            
            break;
         }
      }

      // Reached the end of the alignment without detecting breakpoint.
      if (i == len_q + ALIGN_WIDTH - 1) {
         cell_t * bp = path[path_len];
         // Set alignment values.
         align.start   = align_max(0,bp->row);
         align.end     = align_max(0,bp->col);
         align.score   = bp->score;
         align.pathlen = path_len + 1;
         
         break;
      }
      uL = L;
   }

   // Free path and matrix.
   free(path);
   for (int k = 0; k < allocated; k++) free(Ls[k]);
   free(Ls);

   return align;
}

/*
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

   double logA = log(READ_ERROR_PROB);
   double logB = log(READ_MATCH_PROB);
   double logC = log(RAND_ERROR_PROB);
   double logD = log(RAND_MATCH_PROB);

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

      // Compute extreme values.
      int start_value = align_max(1, i - ALIGN_WIDTH);
      int end_value = align_min(len_q - 1, i + ALIGN_WIDTH);
      if (start_value > 1) {
         matrix[i*len_q + start_value - 1] = (cell_t) {i, 0, 0};
      }
      if (end_value < len_q - 1) {
         matrix[i*len_q + end_value + 1] = (cell_t) {i, 0, 0};
      }

      // Set upper cell as the best.
      int best_idx = len_q*i;
      int minscore = len_q + max_len;

      // Update column.
      for (int j = start_value; j <= end_value; j++) {
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
         if (matrix[idx].score < minscore || (matrix[idx].score == minscore && matrix[idx].inserts > matrix[best_idx].inserts)) {
            best_idx = idx; 
            minscore = matrix[idx].score;
         }
      }
      
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

      // Breakpoint detection.
      if (i >= MIN_ALIGNMENT_LEN) {

         // Compute breakpoint statistics.
         double maxJ = -INFINITY;
         int breakpoint = 0;
         // b represents the last value of the 1st interval. So b is part of the 1st interval!
         for (int b = 0; b < path_len; b++) {
            int E1 = path[path_len].score - path[b].score;
            int M1 = path_len - b - E1;
            float J = E1*(logC-logA) + M1*(logD-logB);
            if (J > maxJ) {
               maxJ = J;
               breakpoint = b;
            }
         }
         
         // DEBUG.
         //         fprintf(stdout, "Alignment:\tlen=%d\tpathlen=%d\tbreakpoint=%d\tJ=%f\n",i,path_len+1,breakpoint,maxJ);

         if (maxJ > BREAKPOINT_THR) {
            // Compute column of maxJ, i.e. last nucleotide of the genome.
            int col = (path[breakpoint+1].path / len_q);
            int c_idx = col * len_q;
            int row = 0;
            // Find minimum score at breakpoint's column.
            int minscore = breakpoint;
            for (int j = align_max(1, col - ALIGN_WIDTH); j <= align_min(len_q - 1, col + ALIGN_WIDTH); j++) {
               if (matrix[c_idx + j].score < minscore) {
                  minscore = matrix[c_idx+j].score;
                  row = j;
               }
            }
            fprintf(stdout, "align resolved: pathlen=%d, score=%d, inserts=%d, row=%d, col=%d, J=%f\n", path_len, minscore, matrix[c_idx+row].inserts, row, col, maxJ);
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
*/

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
