#include "filter.h"

int
align_seeds
(
 char         * read,
 matchlist_t  * seeds,
 matchlist_t ** seqmatches,
 index_t      * index,
 filteropt_t    opt,
 alignopt_t     alignopt
)
{
   if (seeds->pos == 0) return 0;
   int significant = 0;
   int slen = strlen(read);
   matchlist_t * matches = *seqmatches;
   

   // Sort by seeded nucleotides.
   mergesort_mt(seeds->match, seeds->pos, sizeof(match_t), 0, 1, compar_seedhits);

   for(long k = 0; k < seeds->pos; k++) {
      // Get next seed.
      match_t seed = seeds->match[k];

      // Compute alignment limits.
      long r_min  = seed.read_e - seed.read_s + 1;
      long l_min  = 0;
      long r_qstart = seed.read_s + 1;
      long l_qstart = align_max(r_qstart - 1,0);
      long r_rstart = seed.ref_s + 1;
      long l_rstart = seed.ref_s;
     
      // Extend existing mapping.
      double last_eexp = INFINITY;
      int extend = 0;
      int extend_score = 0;
      int cancel_align = 0;
      for (int i = 0; i < matches->pos; i++) {
         match_t match = matches->match[i];
         if (match.e_exp > last_eexp) continue;
         // Check whether the seed is contiguous.
         long r_distr = seed.read_e - match.read_e;
         long g_distr = match.ref_e - seed.ref_e;
         long r_distl = match.read_s - seed.read_s;
         long g_distl = seed.ref_s - match.ref_s;
         long d_accept = opt.dist_accept;
         double rg_ratio = opt.read_ref_ratio;
         double eexp_accept = opt.align_accept_eexp;
         int ext_r = g_distr > -d_accept && r_distr > -d_accept && ((g_distr < (r_distr * rg_ratio)) || (g_distr < d_accept && r_distr < d_accept));
         int ext_l = g_distl > -d_accept && r_distl > -d_accept && ((g_distl < (r_distl * rg_ratio)) || (g_distl < d_accept && r_distl < d_accept));
         if (ext_r || ext_l) {
            r_min  = (ext_r ? r_distr : 0);
            l_min  = (ext_l ? r_distl : 0);
            r_qstart = match.read_e + 1;
            r_rstart = match.ref_e + 1;
            l_qstart = match.read_s - 1;
            l_rstart = match.ref_s - 1;
            last_eexp = match.e_exp;
            extend = i;
            extend_score = match.score;
         } else if (match.e_exp < eexp_accept) {
            // Check overlap
            int span = align_min(seed.read_e - seed.read_s, match.read_e - match.read_s);
            int overlap = align_max(0, align_min(seed.read_e, match.read_e) - align_max(seed.read_s, match.read_s));
            overlap = overlap > span*opt.overlap_max_tolerance;
            int seed_ratio = seed.hits < match.hits*opt.align_seed_filter_thr;
            // If overlap is higher than the maximum overlap tolerance, cancel the alignment.
            if (overlap && seed_ratio) {
               cancel_align = 1;
               break;
            }
         } else {
            if (r_distr < 0) {r_distr = -r_distr; g_distr = -g_distr;}
            if (r_distl < 0) {r_distl = -r_distl; g_distl = -g_distl;}
            int same_align = (g_distr > -d_accept && (g_distr < (r_distr*rg_ratio) || g_distr < d_accept)) || (g_distl > -d_accept && (g_distl < (r_distl*rg_ratio) || g_distl < d_accept));
            if (same_align) {
               cancel_align = 1;
               break;
            }
         }
      }
      long r_qlen = slen - r_qstart;
      long l_qlen = l_qstart + 1;
      if (cancel_align || (!r_qlen && !l_qlen)) continue;
      long r_rlen = align_min((long)(r_qlen * (1 + alignopt.width_ratio)), index->size - r_rstart);
      long l_rlen = align_min((long)(l_qlen * (1 + alignopt.width_ratio)), r_qstart);
      char * r_qry = read + r_qstart;
      char * l_qry = read + l_qstart;
      char * r_ref = index->genome + r_rstart;
      char * l_ref = index->genome + l_rstart;
      
      // Align forward and backward starting from (read_s, ref_s).
      long read_s, read_e, ref_s, ref_e;
      path_t align_r = (path_t){0,0,0}, align_l = (path_t){0,0,0};
      // Forward alignment (right).
      if(r_qlen) {
         align_r = dbf_align_bp(r_qlen, r_qry, r_rlen, r_ref, r_min, ALIGN_FORWARD, ALIGN_FORWARD, alignopt);
         read_e = r_qstart + align_r.row;
         ref_e = r_rstart + align_r.col;
      }
      else {
         read_e = r_qstart - 1;
         ref_e  = r_rstart - 1;
      }
      // Backward alignment (left).
      if(l_qlen) {
         align_l = dbf_align_bp(l_qlen, l_qry, l_rlen, l_ref, l_min, ALIGN_BACKWARD, ALIGN_BACKWARD, alignopt);
         read_s =  l_qstart - align_l.row;
         ref_s = l_rstart - align_l.col;
      }
      else {
         read_s = l_qstart + 1;
         ref_s  = l_rstart + 1;
      }
      
      // Compute significance.
      long score = extend_score + align_l.score + align_r.score;
      double ident = 1.0 - (score*1.0)/(align_max(read_e-read_s+1,ref_e-ref_s+1));
      double e_exp = INFINITY;
      if (ident > opt.align_filter_ident) 
         e_exp = e_value(ref_e - ref_s + 1, extend_score + align_l.score + align_r.score, index->size);

      // If significant, store.
      if (e_exp < opt.align_filter_eexp) {
         match_t hit;
         if (extend) {
            hit = matches->match[extend];
            hit.hits += seed.hits;
         }
         else {
            hit.flags  = 0;
            hit.hits   = seed.hits;
         }

         // Fill/Update hit.
         hit.score  = score;
         hit.ident  = ident;
         hit.ref_e  = ref_e;
         hit.ref_s  = ref_s;
         hit.read_e = read_e;
         hit.read_s = read_s;
         hit.e_exp  = e_exp;
         hit.mapq   = 0;

         // Add to significant matchlist.
         significant = 1;
         if (extend) {
            matches->match[extend] = hit;
         } else {
            matchlist_add(seqmatches, hit);
            matches = *seqmatches;
         }
      }
   }

   return significant;
}


int
compute_mapq
(
 matchlist_t ** intervals,
 int n_int,
 double e_ratio,
 seq_t seq,
 index_t * index
)
{
   for (int i = 1; i < n_int; i++) {
      matchlist_t * interval = intervals[i];
      if (interval->pos < 2) {
         double e_exp = -10*interval->match[0].e_exp - 3;
         interval->match[0].mapq = (e_exp < 60 ? e_exp : 60);
         continue;
      }
      match_t match = interval->match[0];
      double exp_thr = match.e_exp * e_ratio;
      int beg = match.read_s;
      int len = match.read_e - match.read_s + 1;

      double psum = 0.0;
      if (interval->match[1].e_exp < exp_thr) {
         double mapq = mapq_align(seq.seq + beg, seq.q + beg, index->genome + match.ref_s, len, len*0.02);
         // Assign mapq.
         interval->match[0].mapq = pow(10.0, -mapq);
         psum += interval->match[0].mapq;
      } else {
         double e_exp = -10*interval->match[0].e_exp - 3;
         interval->match[0].mapq = (e_exp < 60 ? e_exp : 60);
         continue;
      }
      
      int j = 1;
      for (; j < interval->pos; j++) {
         match_t match = interval->match[j];         
         if (interval->match[j].e_exp < exp_thr) {
            char * ref_ptr = index->genome + match.ref_s + beg - match.read_s;
            double mapq = mapq_align(seq.seq + beg, seq.q + beg, ref_ptr, len, len*0.02);
            interval->match[j].mapq = pow(10.0, -mapq);
            psum += interval->match[j].mapq;
         } else break;
      }
      j--;
      for (; j >= 0; j--) {
         double p_exp = (int)(-10*log10(1-interval->match[j].mapq/psum));
         double e_exp = (int)(-10*interval->match[j].e_exp) - 3;
         double quality = (p_exp < e_exp ? p_exp : e_exp);
         interval->match[j].mapq = quality > 60 ? 60 : quality;
      }
   }
   return 0;
}


matchlist_t **
merge_intervals
(
 matchlist_t * matches,
 double        overlap_ratio,
 int32_t     * nintervals
)
{
   // Sort candidate matches by size.
   mergesort_mt(matches->match, matches->pos, sizeof(match_t), 0, 1, compar_matcheexp);

   // Allocate interval list.
   int32_t n = 0;
   matchlist_t ** intervals = malloc((matches->pos+1)*sizeof(matchlist_t *));

   // Interval 0 is for overlapping matches.
   intervals[n++] = matchlist_new(matches->pos);

   // Compute intervals.
   for (int i = 0; i < matches->pos; i++) {
      match_t match = matches->match[i];
      int new = 1, total_overlap = 0;
      for (int j = 1; j < n; j++) {
         // Compare with the heads of other intervals.
         match_t ref = intervals[j]->match[0];
         // Continue if there is no overlap.
         if (ref.read_e <= match.read_s || match.read_e <= ref.read_s) continue;
         // Match is contained in ref.
         else if (ref.read_s <= match.read_s && match.read_e <= ref.read_e) {
            new = 0;
            matchlist_add(intervals+j, match);
            matches->match[i].interval = j;
            break;
         }
         // Match is partially contained in ref.
         else {
            int overlap = min(match.read_e - ref.read_s, ref.read_e - match.read_s);
            int    span = min(match.read_e - match.read_s, ref.read_e - ref.read_s);
            if (overlap*1.0/span >= overlap_ratio) {
               new = 0;
               matchlist_add(intervals+j, match);
               matches->match[i].interval = j;
               break;
            }
            total_overlap += overlap;
         }
      }
      if (new) {
         if (total_overlap*1.0/(match.read_e - match.read_s) > (1 - overlap_ratio)) {
            // Overlapped interval.
            matchlist_add(intervals, match);
            matches->match[i].interval = 0;
         } else {
            // New interval.
            intervals[n] = matchlist_new(matches->pos-i);
            matchlist_add(intervals + n, match);
            matches->match[i].interval = n;
            n++;
         }
      }
   }

   *nintervals = n;
   // Sort intervals by start position.
   mergesort_mt(intervals+1, n-1, sizeof(matchlist_t *), 0, 1, compar_intvstart);

   return intervals;
}

double
e_value
(
 int L,
 int m,
 long gsize
)
{
   double E = log10(gsize) + log10(2) * (m*3.0-2*L);
   for (int i = 0; i < m; i++) E += log10((L-i)/(double)(m-i));
   return E;
}


int
compar_seedhits
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = (match_t *) a;
   match_t * mb = (match_t *) b;
   
   if (ma->hits < mb->hits) return 1;
   else return -1;
}

int
compar_matcheexp
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = (match_t *) a;
   match_t * mb = (match_t *) b;

   if (mb->e_exp < ma->e_exp) return 1;
   else return -1;
}

int
compar_intvstart
(
 const void * a,
 const void * b,
 const int   param
)
{
   matchlist_t * ma = *(matchlist_t **) a;
   matchlist_t * mb = *(matchlist_t **) b;
   if (ma->pos == 0) return 1;
   if (mb->pos == 0) return -1;
   if (mb->match[0].read_s < ma->match[0].read_s) return 1;
   else return -1;
}
