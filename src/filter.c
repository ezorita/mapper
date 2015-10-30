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

//   if (seeds->pos == 0) return 0;
   int significant = 0;
   int slen = strlen(read);
   
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

      long r_qlen = slen - r_qstart;
      long l_qlen = l_qstart + 1;
      if (!r_qlen && !l_qlen) continue;
      long r_rlen = align_min((long)(r_qlen * (1 + alignopt.width_ratio)),
            index->size - r_rstart);
      long l_rlen = align_min((long)(l_qlen * (1 + alignopt.width_ratio)),
            r_qstart);
      char * r_qry = read + r_qstart;
      char * l_qry = read + l_qstart;
      char * r_ref = index->genome + r_rstart;
      char * l_ref = index->genome + l_rstart;
      
      // Align forward and backward starting from (read_s, ref_s).
      long read_s, read_e, ref_s, ref_e;
      path_t align_r = (path_t){0,0,0}, align_l = (path_t){0,0,0};
      // Forward alignment (right).
      if(r_qlen) {
         align_r = dbf_align_bp(r_qlen, r_qry, r_rlen, r_ref, r_min,
               ALIGN_FORWARD, ALIGN_FORWARD, alignopt);
         read_e = r_qstart + align_r.row;
         ref_e = r_rstart + align_r.col;
      }
      else {
         read_e = r_qstart - 1;
         ref_e  = r_rstart - 1;
      }
      // Backward alignment (left).
      if(l_qlen) {
         align_l = dbf_align_bp(l_qlen, l_qry, l_rlen, l_ref, l_min,
               ALIGN_BACKWARD, ALIGN_BACKWARD, alignopt);
         read_s =  l_qstart - align_l.row;
         ref_s = l_rstart - align_l.col;
      }
      else {
         read_s = l_qstart + 1;
         ref_s  = l_rstart + 1;
      }
      
      // Compute significance.
//      long score = extend_score + align_l.score + align_r.score;

      // If significant, store.
      match_t hit;
      hit.flags  = 0;
      hit.hits   = seed.hits;

      // Fill/Update hit.
      hit.ref_e  = ref_e;
      hit.ref_s  = ref_s;
      hit.read_e = read_e;
      hit.read_s = read_s;
      hit.mapq   = 0;

      // Add to significant matchlist.
      significant = 1;
      matchlist_add(seqmatches, hit);
   }

   return significant;
}


int
compute_mapq
(
 matchlist_t ** intervals,
 int n_int,
 double eval_ratio,
 double err_rate,
 seq_t seq,
 index_t * index
)
{
   for (int i = 1; i < n_int; i++) {
      matchlist_t * interval = intervals[i];
      match_t match = interval->match[0];
      double exp_thr = match.e_exp * eval_ratio;
      int beg = match.read_s;
      int len = match.read_e - match.read_s + 1;
      int err = (int) ((1-match.ident) * len);
      if (interval->pos < 2 || (interval->pos >= 2 && interval->match[1].e_exp > exp_thr)) {
         double e_exp = -10*match.e_exp - 3;
         double ident = 0;
         if (1-match.ident > err_rate)
            ident = -10*(log10(binom(len,err)) + (len-err)*log10(1-err_rate) + err*log10(err_rate));
         double mapq = (e_exp < 60 ? e_exp : 60) - ident;
         interval->match[0].mapq = (mapq < 0 ? 0 : mapq);
         continue;
      }

      if (VERBOSE_DEBUG) {
         fprintf(stdout, "mapq:\n");
      }

      /*
      int errs;
      int mapq = mapq_align(seq.seq + beg, seq.q + beg, index->genome + match.ref_s, len, len*0.02, &errs);
      // Assign mapq.
      interval->match[0].mapq = pow(10.0, -(mapq/10.0)) * binom(len,errs);
      if (VERBOSE_DEBUG) {
         fprintf(stdout, "[0] errors = %d, score = %d, prob = %f\n", errs, mapq, interval->match[0].mapq);
      }
      psum += interval->match[0].mapq;
      */
      double psum = 0.0;
      double * binoms = malloc(interval->pos*sizeof(double));
      int    * n_err = malloc(interval->pos*sizeof(int));
      int j;
      for (j = 0; j < interval->pos; j++) {
         match_t match = interval->match[j];         
         if (match.e_exp < exp_thr) {
            char * ref_ptr = index->genome + match.ref_s + beg - match.read_s;
            int mapq = mapq_align(seq.seq + beg, seq.q + beg, ref_ptr, len, len*0.02, n_err+j);
            binoms[j] = binom(len,n_err[j]);
            interval->match[j].mapq = pow(10.0, -(mapq/10.0)) * binoms[j];
            if (VERBOSE_DEBUG) {
               fprintf(stdout, "[%d] errors = %d, score = %d, prob = %f\n", j, n_err[j], mapq, interval->match[j].mapq);
            }
            psum += interval->match[j].mapq;
         } else break;
      }
      j--;
      for (; j >= 0; j--) {
         double p_div = interval->match[j].mapq/psum;
         double p_exp = (int)(-10*log10(max(1-p_div, 0.0000001)));
         double e_exp = (int)(-10*interval->match[j].e_exp) - 3;
         double ident = 0;
         if (n_err[j]*1.0/len > err_rate)
            ident = (int)(-10*(log10(binoms[j]) + (len-n_err[j])*log10(1-err_rate) + n_err[j]*log10(err_rate)));
         double quality = min(p_exp, e_exp);
         quality = (quality > 60 ? 60 : quality) - ident;
         interval->match[j].mapq = (quality < 0 ? 0 : quality);
      }
   }
   return 0;
}


matchlist_t **
single_interval
(
 matchlist_t * matches
)
{
   // Sort matches by e value.
   mergesort_mt(matches->match, matches->pos, sizeof(match_t), 0, 1, compar_matcheexp);
   
   matchlist_t ** intervals = malloc(2*sizeof(matchlist_t *));
   intervals[0] = matchlist_new(1);
   intervals[1] = matchlist_new(matches->pos);
   
   intervals[1]->pos = matches->pos;
   memcpy(intervals[1]->match, matches->match, matches->pos * sizeof(match_t));

   return intervals;
}

matchlist_t **
merge_intervals
(
 matchlist_t * matches,
 double        overlap_ratio,
 int32_t     * nintervals
)
{
   // Sort candidate matches by e-value.
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

int
filter_repeats
(
 matchlist_t ** intervals,
 matchlist_t  * repeats,
 int32_t        n_ints
)
{
   // Sort repeats by position.
   mergesort_mt(repeats->match, repeats->pos, sizeof(match_t), 0, 1, compar_match_read_beg);
   // Iterate over all intervals.
   for (int i = 1; i < n_ints; i++) {
      matchlist_t * itv = intervals[i];
      match_t ref = itv->match[0];
      int32_t ilen = ref.read_e - ref.read_s + 1;
      // Iterate over all repeats, find total overlap.
      int b, e, ov = 0;
      b = ref.read_s;
      e = b - 1;
      for (int j = 0; j < repeats->pos; j++) {
         int32_t r_beg = repeats->match[j].read_s, r_end = repeats->match[j].read_e;
         if (r_end < ref.read_s) continue;
         if (r_beg > ref.read_e) break;
         if (r_beg > e) {
            ov += e-b+1;
            b = r_beg;
            e = min(r_end, ref.read_e);
         } else {
            e = r_end;
         }
      }
      ov += e-b+1;
      // Apply penalty.
      double p = 1 - ov*1.0/ilen;
      for (int i = 0; i < itv->pos; i++) {
         if (VERBOSE_DEBUG) {
            fprintf(stdout, "Applying repeat penalty (len=%d, ov=%d): %f (mapq=%f, new mapq=%d)\n", ilen, ov, p, itv->match[i].mapq, (int)(itv->match[i].mapq*p));
         }
         itv->match[i].mapq *= p;
         itv->match[i].flags |= FLAG_REPEAT;

      }
   }
   return 0;
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
compar_match_read_beg
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = (match_t *) a;
   match_t * mb = (match_t *) b;

   if (mb->read_s < ma->read_s) return 1;
   else if (mb->read_s > ma->read_s) return -1;
   else return (mb->read_e < ma->read_e);
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
