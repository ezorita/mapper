#include "filter.h"

int
extend_align_bp
(
 match_t    * aln,
 char       * seq,
 size_t       slen,
 index_t    * index,
 alignopt_t   opt
)
{
#if VERBOSE_DEBUG == 1
   char * refseq = malloc(slen+1);
   char * qryseq = malloc(slen+1);
   // L extension.
   int nts = aln->read_s+1;
   memcpy(refseq, index->genome+aln->ref_s-aln->read_s, nts);
   memcpy(qryseq, seq, nts);
   refseq[nts] = 0;
   qryseq[nts] = 0;
   fprintf(stdout, "align L (%d nt):\nqry: %s\nref: %s\n", nts, qryseq, refseq);
   // seed.
   nts = aln->read_e - aln->read_s + 1;
   memcpy(refseq, index->genome+aln->ref_s, nts);
   memcpy(qryseq, seq + aln->read_s, nts);
   refseq[nts] = 0;
   qryseq[nts] = 0;
   fprintf(stdout, "seed (%d nt):\nqry: %s\nref: %s\n", nts, qryseq, refseq);
   // R extension.
   nts = slen-aln->read_e;
   memcpy(refseq, index->genome+aln->ref_e, nts);
   memcpy(qryseq, seq + aln->read_e, nts);
   refseq[nts] = 0;
   qryseq[nts] = 0;
   fprintf(stdout, "align R (%d nt):\nqry: %s\nref: %s\n", nts, qryseq, refseq);
   free(refseq);
   free(qryseq);
#endif

   aln->score = aln->read_e - aln->read_s + 1;

   // Extend right.
   if (aln->read_e < slen-1) {
      // Alignment limits.
      int64_t q_len = slen - aln->read_e;
      int64_t r_len = align_min(q_len, index->size - 1 - aln->ref_e);
      int seed_end_error = translate[(int)seq[aln->read_e]] != translate[(int)index->genome[aln->ref_e]];
      path_t r = align_bp(q_len, seq+aln->read_e, r_len, index->genome+aln->ref_e, 0, 1, 1, opt);
      // Add hits, avoid counting twice the mismatches at the ends of the seed.
      int hits = r.col - max(0,r.score - seed_end_error);
      aln->hits += hits;
      aln->read_e += r.col;
      aln->ref_e += r.row;
      aln->score += max(0, hits - opt.mismatch_penalty*r.score);

#if VERBOSE_DEBUG == 1
      fprintf(stdout, "R: qry:%d, ref:%d, mismatches:%d\n", r.col, r.row, r.score);
#endif

      // Extend deletion gaps.
      if (r.col < r.row && r.row == r_len-1) {
         q_len -= r.col;
         r = align_bp(q_len, seq+aln->read_e, q_len, index->genome+aln->ref_e, 0, 1, 1, opt);
         int hits = r.col - r.score;
         aln->hits += hits;
         aln->read_e += r.col;
         aln->ref_e += r.row;
         aln->score += max(0, hits - opt.mismatch_penalty*r.score);
#if VERBOSE_DEBUG == 1
         fprintf(stdout, "R (del gap): qry:%d, ref:%d, mismatches:%d\n", r.col, r.row, r.score);
#endif
      }
   }
   // Extend left.
   if (aln->read_s > 0) {
      int64_t q_len = aln->read_s + 1;
      int64_t r_len = align_min(q_len, aln->ref_s + 1);
      int seed_end_error = translate[(int)seq[aln->read_s]] != translate[(int)index->genome[aln->ref_s]];
      path_t l = align_bp(q_len, seq+aln->read_s, r_len, index->genome+aln->ref_s, 0, -1, -1, opt);
      int hits = l.col - max(0, l.score - seed_end_error);
      aln->hits += hits;
      aln->read_s -= l.col;
      aln->ref_s -= l.row;
      aln->score += max(0, hits - opt.mismatch_penalty*l.score);

#if VERBOSE_DEBUG == 1
      fprintf(stdout, "L: qry:%d, ref:%d, mismatches:%d\n", l.col, l.row, l.score);
#endif
      // Extend deletion gaps.
      if (l.row > l.col && l.row == r_len-1) {
         q_len -= l.col;
         l = align_bp(q_len, seq+aln->read_s, q_len, index->genome+aln->ref_s, 0, -1, -1, opt);
#if VERBOSE_DEBUG == 1
         fprintf(stdout, "L (del gap): qry:%d, ref:%d, mismatches:%d\n", l.col, l.row, l.score);
#endif
         hits = l.col - l.score;
         aln->read_s -= l.col;
         aln->ref_s -= l.row;
         aln->hits += hits;
         aln->score += max(0, hits - opt.mismatch_penalty*l.score);
      }
   }
#if VERBOSE_DEBUG == 1
   fprintf(stdout, "alignment: beg:%d, end:%d, hits:%d, score:%d\n", aln->read_s, aln->read_e, aln->hits, aln->score);
#endif

   // Return.
   return 0;
}


int
align_seeds
(
 char         * read,
 seedstack_t  * seeds,
 matchlist_t ** seqmatches,
 double       * cumlog,
 index_t      * index,
 filteropt_t    opt,
 alignopt_t     alignopt
)
{
    // Convert seeds to hits.
   uint64_t hit_count;
   // Forward seeds.
   hit_t * hits = compute_hits(seeds, index, &hit_count);
   int rval = align_hits(read,hits,hit_count,seqmatches,cumlog,index,opt,alignopt);
   free(hits);
   return rval;
}

int
add_align_info
(
 match_t        aln,
 int            slen,
 matchlist_t ** matchlist,
 double       * cumlog,
 index_t      * index,
 filteropt_t    opt
)
{
   int ov_beg = -1, ov_end = -1;
   matchlist_t * list = *matchlist;
   // Find overlapping intervals.
   int best = 1;
   for (size_t j = 0; j < list->pos; j++) {
      match_t m = list->match[j];
      if (m.read_e <= aln.read_s) continue;
      if (m.read_s >= aln.read_e) break;
      int overlap = min(aln.read_e - m.read_s, m.read_e - aln.read_s);
      int    span = min(aln.read_e - aln.read_s, m.read_e - m.read_s);
      if (overlap >= span*opt.overlap_max_tolerance) {
         if (aln.hits > m.hits) {
            // aln is better, extend overlap region.
            if (ov_beg < 0) {
               ov_beg = ov_end = j;
            } else {
               ov_end = j;
            }
            // Set second alignment.
            if (m.hits > aln.s_hits) {
               aln.s_hits = m.hits;
               aln.s_hits_cnt  = (m.hits == m.s_hits ? m.s_hits+1 : 1);
            } else if (m.hits == aln.s_hits) {
               aln.s_hits_cnt += 1;
            }
            if (m.score >= m.s_score) {
               aln.s_score = m.score;
               aln.s_score_cnt = 1 + (m.s_score == m.score ? m.s_score_cnt : 0);
            } else {
               aln.s_score = m.s_score;
               aln.s_score_cnt = m.s_score_cnt;
            }
            // Flag mapQ to recompute.
            aln.mapq = -1;
         } else {
            // Conditions to accept a second best hit:
            // 1. s_hits >= best_hits - annotation
            // 2. s_hits >= 'min_best_to_second'*best_hits.
            best = 0;
            if (aln.hits > m.s_hits && aln.hits >= m.hits-m.ann_d && aln.hits > m.hits*opt.min_best_to_second) {
               list->match[j].s_hits = aln.hits;
               list->match[j].s_hits_cnt  = 1;
               // Recompute neighbor mapQ.
               int mapq, maxq;
               mapq = neighbor_mapq(list->match[j], cumlog, &maxq);
               list->match[j].mapq = min(mapq,list->match[j].mapq);
               list->match[j].maxq = min(maxq,list->match[j].maxq);
            } else if (aln.hits == m.s_hits) {
               list->match[j].s_hits_cnt += 1;
               // Recompute neighbor mapQ.
               int mapq, maxq;
               mapq = neighbor_mapq(list->match[j], cumlog, &maxq);
               list->match[j].mapq = min(mapq,list->match[j].mapq);
               list->match[j].maxq = min(maxq,list->match[j].maxq);
            }
         }
         if (aln.score > m.s_score) {
            list->match[j].s_score = aln.score;
            list->match[j].s_score_cnt = 1;
            // Recompute score mapQ.
            int mapq;
            mapq = score_mapq(list->match[j]);
            list->match[j].mapq = min(mapq,list->match[j].mapq);
            list->match[j].maxq = min(mapq,list->match[j].maxq);
         } else if (aln.score == m.s_score) {
            list->match[j].s_score_cnt += 1;
            // Recompute score mapQ.
            int mapq;
            mapq = score_mapq(list->match[j]);
            list->match[j].mapq = min(mapq,list->match[j].mapq);
            list->match[j].maxq = min(mapq,list->match[j].maxq);
         }
      }
   }
   // Insert interval.
   if (best) {
      // Add sequence neighborhood info.
      check_neighbor_info(&aln,index);
      // Double check second best alignment.
      if (aln.s_hits < aln.hits - aln.ann_d) aln.s_hits = aln.s_hits_cnt =  0;
      // Compute tentative mapQ.
      tentative_mapq(&aln, cumlog);
      // Insert interval in empty gap.
      if (ov_beg < 0) {
         // Insert position.
         int ins = 0;
         for (;ins < list->pos && list->match[ins].read_s < aln.read_s; ins++);
         // Add match.
         matchlist_add(matchlist, aln);
         list = *matchlist;
         // Move memory.
         if (ins < list->pos - 1) {
            memmove(list->match+ins+1, list->match+ins, (list->pos - 1 - ins)*sizeof(match_t));
         }
         list->match[ins] = aln;
      }
      // Substitute interval(s).
      else {
         list->match[ov_beg] = aln;
         if (ov_end > ov_beg) {
            if (ov_end < list->pos - 1)
               memmove(list->match+ov_beg + 1, list->match+ov_end+1, (list->pos-1-ov_end)*sizeof(match_t));
            list->pos += ov_beg - ov_end;
         }
      }
      // Recompute gap identities.
      int beg = 0;
      for (int i = 0; i < list->pos; i++) {
         match_t m = list->match[i];
         int gap_span = (i == list->pos - 1 ? slen : list->match[i+1].read_s) - beg;
         list->match[i].gap_id = m.hits*1.0/gap_span;
         beg = m.read_e + 1;
      }
   }

   return 0;
}

int
check_neighbor_info
(
 match_t * aln,
 index_t * index
)
{
   // Undo reverse complement.
   int64_t span = aln->ref_e - aln->ref_s + 1;
   int64_t ref_start = aln->ref_s;
   if (ref_start >= (index->size >> 1)) {
      ref_start = index->size - 2 - aln->ref_e;
   }

   // Use only the largest available k-mer.
   int idx = -1, k = 0;
   for (int a = 0; a < index->ann->count; a++) {
      ann_t ann = index->ann->ann[a];
      if (span - ann.data->kmer > 0 && ann.data->kmer > k) {
         idx = a;
         k = ann.data->kmer;
      }
   }

   // Add closest neighbor distance,
   int max_d = 0, cnt_d = 0, last = -k, span_d = 0;
   if (idx >= 0) {
      ann_t ann = index->ann->ann[idx];
      for (int64_t n = 0; n <= span - ann.data->kmer; n++) {
         int log10cnt;
         int d = ann_read(ann, ref_start + n, &log10cnt); 
         if (d > max_d) {
            max_d = d;
            cnt_d = log10cnt;
            span_d = ann.data->kmer;
            last = n;
         } else if (d == max_d) {
            if (log10cnt > cnt_d) cnt_d = log10cnt;
            span_d += min(n - last, ann.data->kmer);
            last = n;
         }
      }
      aln->ann_d   = (int)max_d*(span_d*1.0/ann.data->kmer);
      aln->ann_cnt = 1;
      for (int n = 0; n < cnt_d; n++) aln->ann_cnt *= 10;
   } else {
      aln->ann_d   = 0;
      aln->ann_cnt = 0;
   }
   
   return 0;
}

int
align_hits
(
 char         * read,
 hit_t        * hits,
 uint64_t       hit_count,
 matchlist_t ** seqmatches,
 double       * cumlog,
 index_t      * index,
 filteropt_t    opt,
 alignopt_t     alignopt
)
{
   int slen = strlen(read);
   for(int i = 0; i < hit_count; i++) {
      int done = 0;
      for (size_t j = 0; j < (*seqmatches)->pos; j++) {
         match_t m = (*seqmatches)->match[j];
         // Filter out repeated alignments.
         // TODO:
         //  - Needs improvement, keep track of all the alignments.
         int64_t d = ((int64_t)m.ref_s) - ((int64_t)hits[i].locus);
         if (d < slen && d > -slen) {
            if (VERBOSE_DEBUG) fprintf(stdout, "repeated alignment\n");
            done = 1;
            break;
         }
         // Filter out redundant alignments.
         /*
         if (hits[i].beg >= m.read_s && hits[i].end <= m.read_e) {
            if ((m.mapq > opt.target_q || m.maxq < opt.target_q) && m.gap_id >= opt.min_gap_id) {
               if (VERBOSE_DEBUG) fprintf(stdout, "redundant alignment\n");
               done = 1;
               break;
            }
         }
         */
      }
      if (done) continue;

      match_t m = {
         .read_s      = hits[i].s_beg,
         .read_e      = hits[i].s_end,
         .ref_s       = hits[i].locus,
         .ref_e       = hits[i].locus + hits[i].s_end - hits[i].s_beg,
         .hits        = hits[i].s_end - hits[i].s_beg + 1 - hits[i].s_errors,
         .s_hits      = 0,
         .s_hits_cnt  = 0,
         .score       = 0,
         .s_score     = 0,
         .s_score_cnt = 0,
         .mapq        = -1,
         .maxq        = -1
      };
      // Extend seed.
      extend_align_bp(&m,read,slen,index,alignopt);

      // DEBUG.
      if (VERBOSE_DEBUG) fprintf(stdout, "seed extend: ref_s: %ld, read_s: %d, read_e: %d, hits:%d\n", m.ref_s, m.read_s, m.read_e, m.hits);

      // Filter out non-significant alignments.
      if (m.hits < opt.min_interval_hits) continue;

      // Now check overlap with any other interval.
      add_align_info(m,slen,seqmatches,cumlog,index,opt);
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
/*
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
*/
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
   //   return (ma->e_exp > mb->e_exp ? -1 : 1);
   return (ma->mapq > mb->mapq ? -1 : 1);
}
