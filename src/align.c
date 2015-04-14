#include "align.h"

align_t
nw_align
(
 const char * query,
 const char * ref,
 const int len_q,
 const int dir_q,
 const int dir_r,
 const alignopt_t opt
)
{

   if (len_q < 1) {
      // Return alignment at position 0.
      return (align_t) {0, 0, 0, 0};
   }
   // Return value.
   // The column (genome breakpoint will be returned in .end)
   // The row (read breakpoint will be returned in .start)
   align_t align = (align_t) {0,0,0,0};
   // Translation table.
   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

#define PEN_SUB 0
#define PEN_INS 1
#define PEN_DEL 2

   // Matrix cell structure.
   typedef struct cell_t cell_t;
   struct cell_t {
      int   score;
      short ins;
      short dels;
      int   row;
      int   col;
      cell_t * path;
   };

   // Compute alignment sizes.
   int align_max =  (int) (opt.border_y0 + len_q * opt.border_slope);
   int max_plen  = 2 * len_q + align_max;
   int * align_len = malloc((len_q + align_max + 1) * sizeof(int));

   for (int i = 0; i < len_q + align_max + 1; i++) {
      align_len[i] = (int) (opt.border_y0 + i * opt.border_slope);
   }

   // Allocate path and alignment matrix.
   int allocated = align_min(len_q + align_max + 1, ALLOC_BLOCK_SIZE);
   cell_t ** Ls = malloc((len_q + align_max + 1) * sizeof(cell_t*));
   cell_t ** path = malloc(max_plen * sizeof(cell_t*));

   for (int i = 0; i < allocated; i++) {
      Ls[i] = malloc((2*align_len[i]+3) * sizeof(cell_t));
   }
   for (int i = 0 ; i < max_plen; i++) path[i] = NULL;

   // Translated sequences buffer.
   char * val_q = malloc(len_q);
   char * val_r = malloc(len_q + align_max);

   for (int i = 0; i < len_q; i++) 
      val_q[i] = translate[(int)query[dir_q*i]];

   // Aux vars.
   int match, q_gap, r_gap;
   int ins, dels, min_score;
   int path_len = 0;

   int lastbs = 0, lastbe = 0, bp_count = 0, align_end = 0;
   double logAe = log(opt.rand_error) - log(opt.read_error);
   double logBe = log(1-opt.rand_error) - log(1-opt.read_error);

   double logAs = log(opt.rand_subst) - log(opt.read_subst);
   double logBs = log(1-opt.rand_subst) - log(1-opt.read_subst);


   // Initialize (0,0) value.
   cell_t * uL = Ls[0] + align_len[0] + 1;
   uL[0] = (cell_t) {0, 0, 0, 0};

   // Compute inverted Ls.
   for (int i = 0; i < len_q + align_max; i++) {

      // Translate.
      val_r[i] = translate[(int)ref[dir_r*i]];
      int l = i + 1;
      // Expand matrix.
      if (l >= allocated) {
         allocated = align_min(len_q + align_max + 1, allocated + ALLOC_BLOCK_SIZE);
         for (int k = l; k < allocated; k++) {
            Ls[k] = malloc((2*align_len[k]+3) * sizeof(cell_t));
         }
      }
      cell_t * L = Ls[l] + align_len[l] + 1;

      int width, best_idx, best_score;
      cell_t * p;

      // Matrix border elements.
      if (i <= align_len[l]) {
         width = i;
         // Vertical L border.
         p = uL + width;
         L[width+1] = (cell_t) {p->score + 1, p->ins, p->dels + 1, i - width - 1, i, p};
         // Horizontal L border.
         if (i < len_q) {
            p = uL - width;
            L[-width-1] = (cell_t) {p->score + 1, p->ins + 1, p->dels, i, i - width - 1, p};
         }
      } else {
         width = align_len[l];
         // Vertical L border.
         if (width > align_len[l-1]) {
            p = uL + width;
            L[width+1] = (cell_t) {p->score + 1, p->ins, p->dels + 1, i - width - 1, i, p};
         }
         else {
            // Compute score.
            q_gap = uL[width].score + 1;
            match = uL[width+1].score + (val_q[i-width-1] != val_r[i]);
            // Add path connectors.
            min_score = align_min(q_gap, match);
            if (match == min_score) { p = uL + width + 1; ins = p->ins; dels = p->dels; }
            else { p = uL + width; ins = p->ins; dels = p->dels + 1; }
            // Store cell.
            L[width+1] = (cell_t) {min_score, ins, dels, i - width - 1, i, p};
         }
         // Horizontal L border.
         if (i < len_q) {
            if (width > align_len[l-1]) {
               p = uL - width;
               L[-width-1] = (cell_t) {p->score + 1, p->ins + 1, p->dels, i, i - width - 1, p};
            }
            else {
               // Compute score.
               r_gap = uL[-width].score + 1;
               match = uL[-width-1].score + (val_q[i] != val_r[i-width-1]);
               // Add path connectors.
               min_score = align_min(r_gap, match);
               if (match == min_score) { p = uL - width - 1; ins = p->ins; dels = p->dels; }
               else { p = uL - width; ins = p->ins + 1; dels = p->dels; }
               // Store cell.
               L[-width-1] = (cell_t) {min_score, ins, dels, i, i - width - 1, p};
            }
         }
      }

      best_idx = width+1;
      best_score = L[best_idx].score;

      // Compute L elements.
      for (int j = width ; j > align_max(0, i - len_q); j--) {
         // Vertical wing.

         // Compute score.
         match  = uL[j].score + (val_q[i-j] != val_r[i]);
         q_gap = uL[j-1].score + 1;
         r_gap = L[j+1].score + 1;
         // Add path connectors.
         min_score = align_min(align_min(q_gap, r_gap), match);
         if (match == min_score) { p = uL + j; ins = p->ins; dels = p->dels; }
         else if (r_gap == min_score) { p = L + j + 1; ins = p->ins + 1; dels = p->dels; }
         else { p = uL + j - 1; ins = p->ins; dels = p->dels + 1; }
         // Store cell.
         L[j] = (cell_t) {min_score, ins, dels, i - j, i, p};
         // Update best score.
         if (min_score <= best_score) {
            best_score = min_score;
            best_idx = j;
         }

         // Horizontal wing.
         if (i < len_q) {
            // Compute score.
            match = uL[-j].score + (val_q[i] != val_r[i-j]);
            r_gap = uL[-j+1].score + 1;
            q_gap = L[-j-1].score + 1;
            // Add path connectors.
            min_score = align_min(align_min(q_gap, r_gap), match);
            if (match == min_score) { p = uL - j; ins = p->ins; dels = p->dels; }
            else if (r_gap == min_score) { p = uL - j + 1; ins = p->ins + 1; dels = p->dels; }
            else { p = L - j - 1; ins = p->ins; dels = p->dels + 1; }
            // Store cell.
            L[-j] = (cell_t) {min_score, ins, dels, i, i - j, p};
            // Update best score.
            if (min_score <= best_score) {
               best_score = min_score;
               best_idx = -j;
            }
         }
      }

      // Center cell.
      if (i < len_q) {
         // Compute score.
         match = uL[0].score + (val_q[i] != val_r[i]);
         r_gap = L[1].score + 1;
         q_gap = L[-1].score + 1;
         // Add path connectors.
         min_score = align_min(align_min(q_gap, r_gap), match);
         if (match == min_score) { p = uL; ins = p->ins; dels = p->dels; }
         else if (r_gap == min_score) { p = L + 1; ins = p->ins + 1; dels = p->dels; }
         else { p = L - 1; ins = p->ins; dels = p->dels + 1; }
         // Store cell.
         L[0] = (cell_t) {min_score, ins, dels, i, i, p};
         // Update best score.
         if (min_score <= best_score) {
            best_score = min_score;
            best_idx = 0;
         }
      }

      // Update path.
      path_len = L[best_idx].ins + L[best_idx].col;
      path[path_len] = L + best_idx;

      // Recompute path if broken.
      for (int p = path_len; p > 0; p--) {
         if (path[p-1] == NULL || path[p]->path != path[p-1]) {
            // Recover vector position from cellid.
            path[p-1] = path[p]->path;
         }
         else break;
      }

      // Breakpoint detection.
      if (i % opt.bp_period == 0) {
         double maxJe, maxJs;
         // Compute breakpoint statistics.
         maxJe = maxJs = -INFINITY;
         int be = 0, bs = 0;
         // b represents the last value of the 1st interval. So b is part of the 1st interval!
         for (int b = 0; b < path_len; b += opt.bp_period) {
            int Ee = path[path_len]->score - path[b]->score;
            int Me = path_len - b - Ee;
            float Je = Ee*logAe + Me*logBe;
            if (Je > maxJe) { maxJe = Je; be = b; }
            int Es = Ee - path[path_len]->ins + path[b]->ins - path[path_len]->dels + path[b]->dels;
            int Ms = path_len - b - Es;
            float Js = Es*logAs + Ms*logBs;
            if (Js > maxJs) { maxJs = Js; bs = b; }
         }

         fprintf(stdout, "ML algorithm {pos = %d, be = %d (Je = %.2f), bs = %d (Js = %.2f)}\n", path_len, be, maxJe, bs, maxJs);
         
         if (maxJs > opt.bp_thr && maxJe > opt.bp_thr) {
            if (bs == lastbs && be == lastbe) bp_count++;
            else {
               bp_count = 1;
               lastbs = bs;
               lastbe = be;
            }
            if (bp_count >= opt.bp_repeats) {
               align_end = 1;
            }
         }
      }

      // End of alignment, compute ML breakpoint.
      if (L[best_idx].row >= len_q - 1 || align_end) {
         float maxJ = -INFINITY;
         int bp = 0;
         for (int b = 0; b < path_len; b++) {
            int E = path[path_len]->score - path[b]->score;
            int M = path_len - b - E;
            float J = E*logAe + M*logBe;
            if (J > maxJ) { maxJ = J; bp = b; }
         }

         // Set alignment values.
         cell_t * bpc = path[bp];
         align.start   = align_max(0,bpc->row);
         align.end     = align_max(0,bpc->col);
         align.score   = bpc->score;
         align.pathlen = bp + 1;
         
         cell_t * last = path[path_len];
        
         fprintf(stdout, "Breakpoint: %ld\nRead(0-%d): %d (%.2f%%)\n- sub: %d\n- ins: %d\n- del: %d\nRandom(%d-%d): %d (%.2f%%)\n- sub: %d\n- ins: %d\n- del: %d\n", align.start, bp, bpc->score, bpc->score*100.0/(bp+1), bpc->score - bpc->ins - bpc->dels, bpc->ins, bpc->dels, bp, path_len, last->score - bpc->score, (last->score - bpc->score)*100.0/(path_len-bp+1), last->score - last->ins - last->dels - bpc->score + bpc->ins + bpc->dels, last->ins - bpc->ins, last->dels - bpc->dels);

         break;
      }

      uL = L;
   }

   // Free data.
   free(align_len);
   free(val_q);
   free(val_r);
   for (int k = 0; k < allocated; k++) free(Ls[k]);
   free(Ls);


   return align;
}


/*
align_t
nw_align_lite
(
 const char * query,
 const char * ref,
 const int len_q,
 const int len_g,
 const int dir_q,
 const int dir_r,
 const double border_slope,
 const double border_y0,
 const double likelihood_thr,
 const double read_match_prob,
 const double rand_match_prob
)
{

   if (len_q < 1) {
      // Return alignment at position 0.
      return (align_t) {0, 0, 0, 0};
   }
   // Return value.
   // The column (genome breakpoint will be returned in .end)
   // The row (read breakpoint will be returned in .start)
   align_t align = (align_t) {0,0,0,0};
   // Translation table.
   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

#define PEN_SUB 0
#define PEN_INS 1
#define PEN_DEL 2

   char penalty[8][3] = {
      {-34, -29, -29}, // 000 prev = match
      {-17, -36, -38}, // 001 prev = sub
      {-40, -18, -99}, // 010 prev = ins
      {-17, -18, -38}, // 011 prev = sub | ins
      {-43, -99, -22}, // 100 prev = del
      {-17, -36, -22}, // 101 prev = sub | del
      {-40, -18, -22}, // 110 prev = ins | del
      {-17, -18, -22}  // 111 prev = ins | del | sub
   };

   // Matrix cell structure.
   typedef struct cell_t cell_t;
   struct cell_t {
      int score;
      int prev; // 001 mismatch, 010 for insertion, 100 for deletion.
      int path;
   };

   cell_t best;

   // Compute alignment sizes.
   int align_max =  (int) (border_y0 + len_q * border_slope);
   int * align_len = malloc((len_q + align_max + 1) * sizeof(int));
   for (int i = 0; i < len_q + align_max + 1; i++) {
      align_len[i] = (int) (border_y0 + i * border_slope);
   }

   // Allocate path and alignment matrix.
   int max_plen = (3*(len_q+align_max))/2;
   cell_t ** path = malloc(max_plen*sizeof(cell_t *));

   cell_t * m = malloc((align_max + 1) * (len_q + 1) * sizeof(cell_t));

   // Translated sequences buffer.
   char val_q[len_q];
   char val_r[len_q + align_max];

   for (int i = 0; i < len_q; i++) 
      val_q[i] = translate[(int)query[dir_q*i]];

   // Aux vars.
   int match, r_gap, r_gap;

   // Initialize (0,0) value.
   m[0] = (cell_t) {0, 0, 0};

   // Initialize path.
   for (int i = 0; i < max_plen; i++) path[i] = NULL;

   // Compute inverted Ls.
   for (int i = 0; i < len_q + align_max; i++) {
      // Translate.
      val_r[i] = translate[(int)ref[dir_r*i]];

      cell_t * L = Ls[l] + align_len[l] + 1;
      
#define idx(r,c,w) (r*w + c)

      int width;
      if (i > align_len[i]) {
         width = align_len[i];
         m[idx(l-width-1, l, width)] = {len_q, 0, 0};
         m[idx(l, l-width-1, width)] = {len_q, 0, 0};
      } else {
         width = i;
         cell_t l, r;
         l = m[idx(0,width+1,width)];
         m[idx(0,width+1,width)] = (cell_t) {l.score + penalty[l.prev][PEN_DEL], 0x2, idx(0,l-1,width)}
         r = m[idx(i+1,0,width)];
         m[idx(l,0,width)] = (cell_t) {r.score + penalty[r.prev][PEN_INS], 0x4, idx(l-1,0,width)}
      }

      for (int j = 0; j < 

      // Matrix border elements.
      if (i <= align_len[l]) {
         width = i;
         // Vertical L border.
         L[width+1] = (cell_t) {uL[width].score + 1, uL[width].r_gaps + 1, uL[width].g_gaps, i - width - 1, i, uL + width};
         // Horizontal L border.
         if (i < len_q) {
            L[-width-1] = (cell_t) {uL[-width].score + 1, uL[-width].r_gaps, uL[-width].g_gaps + 1, i, i - width - 1, uL - width};
         }
      } else {
         width = align_len[l];
         if (width > align_len[l-1]) 
            L[width+1] = (cell_t) {uL[width].score + 1, uL[width].r_gaps + 1, uL[width].g_gaps, i - width - 1, i, uL + width};
         else {
            // Vertical L border.
            r_gap = uL[width].score + 1;
            match = uL[width+1].score + (val_q[i-width-1] != val_r[i]);
            if (match > r_gap)
               L[width+1] = (cell_t) {r_gap, uL[width].r_gaps + 1, uL[width].g_gaps, i - width - 1, i, uL + width};
            else
               L[width+1] = (cell_t) {match, uL[width+1].r_gaps, uL[width+1].g_gaps, i - width - 1, i, uL + width + 1};
         }
         // Horizontal L border.
         if (i < len_q) {
            if (width > align_len[l-1])
               L[-width-1] = (cell_t) {uL[-width].score + 1, uL[-width].r_gaps, uL[-width].g_gaps + 1, i, i - width - 1, uL - width};
            else {
               g_gap = uL[-width].score + 1;
               match = uL[-width-1].score + (val_q[i] != val_r[i-width-1]);
               if (match > g_gap)
                  L[-width-1] = (cell_t) {g_gap, uL[-width].r_gaps, uL[-width].g_gaps + 1, i, i - width - 1, uL - width};
               else
                  L[-width-1] = (cell_t) {match, uL[-width-1].r_gaps, uL[-width-1].g_gaps, i, i - width - 1, uL - width - 1};
       
            }
         }
      }

      // Compute L elements.
      for (int j = width ; j > align_max(0, i - len_q); j--) {
         // Vertical wing.
         match  = uL[j].score + (val_q[i-j] != val_r[i]);
         r_gap = uL[j-1].score + 1;
         g_gap = L[j+1].score + 1;

         // Assume match. (Prioritize diagonal transitions)
         L[j] = (cell_t) {match, uL[j].r_gaps, uL[j].g_gaps, i - j, i, uL + j};
         // Check g_gap.
         if (g_gap < L[j].score)
            L[j] = (cell_t) {g_gap, L[j+1].r_gaps, L[j+1].g_gaps + 1, i - j, i, L + j + 1};
         // Check r_gap.
         if (r_gap < L[j].score)
            L[j] = (cell_t) {r_gap, uL[j-1].r_gaps + 1, uL[j-1].g_gaps, i - j, i, uL + j - 1};

         // Horizontal wing.
         if (i < len_q) {
            match = uL[-j].score + (val_q[i] != val_r[i-j]);
            r_gap = L[-j-1].score + 1;
            g_gap = uL[-j+1].score + 1;

            // Assume match. (Prioritize diagonal transitions)
            L[-j] = (cell_t) {match, uL[-j].r_gaps, uL[-j].g_gaps, i, i - j, uL - j};
            // Check g_gap.
            if (g_gap < L[-j].score)
               L[-j] = (cell_t) {g_gap, uL[-j+1].r_gaps, uL[-j+1].g_gaps + 1, i, i - j, uL - j + 1};
            // Check r_gap.
            if (r_gap < L[-j].score)
               L[-j] = (cell_t) {r_gap, L[-j-1].r_gaps + 1, L[-j-1].g_gaps, i, i - j, L - j - 1};
         }
      }

      // Center cell.
      if (i < len_q) {
         match = uL[0].score + (val_q[i] != val_r[i]);
         r_gap = L[-1].score + 1;
         g_gap = L[1].score + 1;
         
         // Assume match. (Prioiritize diagonal transitions)
         L[0] = (cell_t) {match, uL[0].r_gaps, uL[0].g_gaps, i, i, uL};
         // Check g_gap.
         if (g_gap < L[0].score)
            L[0] = (cell_t) {g_gap, L[1].r_gaps, L[1].g_gaps + 1, i, i, L + 1};
         // Check r_gap.
         if (r_gap < L[0].score)
            L[0] = (cell_t) {r_gap, L[-1].r_gaps + 1, L[-1].g_gaps, i, i, L - 1};

      }

      // Find best score (choose the closest to the diagonal).
      int best_idx = width;
      int best_score = L[best_idx].score;

      if (i >= len_q) {
         for (int j = width; j >= i - len_q + 1; j--) {
            if (L[j].score <= best_score) {
               best_score = L[j].score;
               best_idx = j;
            }
         }
      } else {
         for (int j = width; j > 0; j--) {
            if (L[-j].score <= best_score) {
               best_score = L[-j].score;
               best_idx = -j;
            }
            if (L[j].score <= best_score) {
               best_score = L[j].score;
               best_idx = j;
            }

         }
         if (L[0].score <= best_score) best_idx = 0;
      }

      // Update path.
      int path_pos = L[best_idx].r_gaps + L[best_idx].row;
      path[path_pos] = L + best_idx;

      // Recompute path if broken.
      for (int p = path_pos; p > 0; p--) {
         if (path[p-1] == NULL || path[p]->path != path[p-1]) {
            // Recover vector position from cellid.
            path[p-1] = path[p]->path;
         }
         else break;
      }
      
      best = L[best_idx];

      if (L[best_idx].row >= len_q - 1) break;

      double maxJ;
      double logA = log(1-rand_match_prob) - log(1-read_match_prob);
      double logB = log(rand_match_prob) - log(read_match_prob);
     
      // Breakpoint detection.
      if (i >= align_min) {

         // Compute breakpoint statistics.
         maxJ = -INFINITY;
         int breakpoint = 0;
         // b represents the last value of the 1st interval. So b is part of the 1st interval!
         for (int b = 0; b < path_pos; b++) {
            int E1 = path[path_pos]->score - path[b]->score;
            int M1 = path_pos - b - E1;
            float J = E1*logA + M1*logB;
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
         cell_t * bp = path[path_pos];
         // Set alignment values.
         align.start   = align_max(0,bp->row);
         align.end     = align_max(0,bp->col);
         align.score   = bp->score;
         align.pathlen = path_pos + 1;
         
         break;
      }

      uL = L;
   }

   int path_len = best.r_gaps + best.row + 1;
   int matches = path_len - best.score;
   int errors  = best.score - best.g_gaps - best.r_gaps;
   int inserts = best.g_gaps;
   int deletes = best.r_gaps;

   fprintf(stdout, "Alignment:\n- Length: %d\n- Score: %d\n- Matches: %d/%d (%.2f%%)\n- Mismatches: %d/%d (%.2f%%)\n- Insertions: %d/%d (%.2f%%)\n- Deletions: %d/%d (%.2f%%)\n", path_len, best.score, matches, path_len, matches/(float)path_len * 100.0, errors, path_len, errors/(float)path_len*100.0, inserts, path_len, inserts/(float)path_len*100.0, deletes, path_len, deletes/(float)path_len*100.0);

   // Free path and matrix.
free(align_len);
   free(path);
   for (int k = 0; k < allocated; k++) free(Ls[k]);
   free(Ls);


   return align;
}
*/
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
/*
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
*/
