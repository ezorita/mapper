#include "filter.h"

int
align_matches
(
 char * read,
 matchlist_t *  hitlist,
 matchlist_t ** seqmatches,
 index_t * index,
 filteropt_t opt,
 alignopt_t alignopt
)
{
   int slen = strlen(read);
   //   mergesort_mt(hitlist->match, hitlist->pos, sizeof(match_t), 0, 1, compar_seedhits);
   hit_t hits[20];
   int best = 0, second = 0;
   for(int i = 0; i < min(20,hitlist->pos); i++) {
      match_t hit = hitlist->match[i];
      hits[i] = (hit_t){.locus = hit.ref_s, .qrypos = hit.read_s, .depth = hit.read_e - hit.read_s + 1, .bulk = 0};
      // Compute alignment limits.
      //long r_min  = seed.read_e - seed.read_s + 1;
      long r_min  = hits[i].depth + 1;
      long l_min  = 0;
      //long r_qstart = seed.read_s + 1;
      long r_qstart = hits[i].qrypos + 1;
      long l_qstart = align_max(hits[i].qrypos,0);
      //long r_rstart = seed.ref_s + 1;
      long r_rstart = hits[i].locus + 1; 
      //long l_rstart = seed.ref_s;
      long l_rstart = hits[i].locus; 

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

      // Find whether alignment is safe in repeat annotation.
      int safe = 0;
      const int safe_tau = 1;
      ann_t ann = index->ann->ann[0];
      uint64_t ref_start = l_rstart - l_rlen + 1;
      if (ref_start >= index->size/2) {
         ref_start = index->size - 1 - ref_start - 25;
         for (uint64_t j = ref_start; j >= ref_start - slen + 25 && !safe; j--) {
            uint64_t bit = (ann.d+1)*j+safe_tau;
            safe = (ann.data[bit>>3] >> (bit&7)) & 1;
         }
      } else {
         for (uint64_t j = ref_start; j <= ref_start + slen - 25 && !safe; j++) {
            uint64_t bit = (ann.d+1)*j+safe_tau;
            safe = (ann.data[bit>>3] >> (bit&7)) & 1;
         }
      }
      
      // Align forward and backward starting from (read_s, ref_s).
      long read_s, read_e, ref_s, ref_e;
      path_t align_r = (path_t){0,0,0}, align_l = (path_t){0,0,0};
      // Forward alignment (right).
      if(r_qlen) {
         align_r = align_bp(r_qlen, r_qry, r_rlen, r_ref, r_min,
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
         align_l = align_bp(l_qlen, l_qry, l_rlen, l_ref, l_min,
               ALIGN_BACKWARD, ALIGN_BACKWARD, alignopt);
         read_s =  l_qstart - align_l.row;
         ref_s = l_rstart - align_l.col;
      }
      else {
         read_s = l_qstart + 1;
         ref_s  = l_rstart + 1;
      }
      int align_span = align_r.row + align_l.row;
      int identities = align_span - align_r.score - align_l.score;
      int ref_min = (*seqmatches)->match[0].ref_s - align_span;
      int ref_max = ref_min + 2*align_span;

      
      // Compute significance.
//      long score = extend_score + align_l.score + align_r.score;
      if (identities > best) {
         // If significant, store.
         match_t hit;
         second = best;
         hit.interval = (safe ? safe : (second ? -1 : 0));
         //         hit.interval = safe;
         hit.ident  = identities;

         // Fill/Update hit.
         hit.ref_e  = ref_e;
         hit.ref_s  = ref_s;
         hit.read_e = read_e;
         hit.read_s = read_s;

         // Mapping quality based on annotation table:
         double mapq = 0.0;
         size_t span = read_e - read_s;
         for (int a = 0; a < index->ann->count; a++) {
            ann = index->ann->ann[a];
            int bits = 0;
            while ((ann.d >> bits) > 0) bits++;
            bits += 2;
            uint16_t mask = 0xFFFF;
            mask >>= (16-bits);
            for (uint64_t j = ref_start; j <= ref_start + span - ann.k && !safe; j++) {
               int log10cnt;
               int d =ann_read(ann, j, &log10cnt); 
               /*
               uint64_t word = bit >> 3;
               uint8_t shift = bit & 7;
               uint16_t w = (((uint16_t)ann.data[word]) >> shift) | (((uint16_t)ann.data[word+1]) << (8-shift));
               int d        = w & (mask >> 2);
               int log10cnt = (w >> (bits-2)) & 3;
               */
               const double Perr = 0.01;
               // This is the first approximation, we are not taking into account the alignment
               // mismatches in this kmer, so it is a quite optimistic approximation.
               mapq = max(mapq, -10*(log10cnt+log10(d)*log10(Perr)-log10(binom(span,d))-d*log10(3)));
            }
         }
         hit.mapq = mapq;
         (*seqmatches)->match[0] = hit;
         (*seqmatches)->pos = 1;
         // update best.
         best = identities;
      } else if (identities >= second && (ref_min >= ref_s || ref_max <= ref_s)) {
         if (second == identities) {
            (*seqmatches)->match[0].s_cnt++;
         } else {
            second = identities;
            (*seqmatches)->match[0].s_hits = identities;
            (*seqmatches)->match[0].s_cnt = 1;
         }
      }
   }
   return 0;
}

int
align_hits
(
 char * read,
 hit_t * hits,
 uint64_t hit_count,
 matchlist_t ** seqmatches,
 index_t * index,
 filteropt_t opt,
 alignopt_t alignopt
)
{
   int slen = strlen(read);
   int best = 0, s_hits = 0, s_cnt = 1;
   if ((*seqmatches)->pos > 0) {
      best = (*seqmatches)->match[0].hits;
      s_hits = (*seqmatches)->match[0].s_hits;
      s_cnt = (*seqmatches)->match[0].s_cnt;
   }
   for(int i = 0; i < hit_count; i++) {
      // Compute alignment limits.
      //long r_min  = seed.read_e - seed.read_s + 1;
      long r_min  = hits[i].depth + 1;
      long l_min  = 0;
      //long r_qstart = seed.read_s + 1;
      long r_qstart = hits[i].qrypos + 1;
      long l_qstart = align_max(r_qstart - 1,0);
      //long r_rstart = seed.ref_s + 1;
      long r_rstart = hits[i].locus + 1; 
      //long l_rstart = seed.ref_s;
      long l_rstart = hits[i].locus; 

      // DEBUG.
      if (VERBOSE_DEBUG) {
         char * rf = malloc(slen+1);
         memcpy(rf, index->genome+hits[i].locus-hits[i].qrypos, slen);
         rf[slen] = 0;
         fprintf(stdout, "sequence alignment [qrypos:%d, locus:%ld]:\nq:%s\nr:%s\n",hits[i].qrypos, hits[i].locus,read,rf);
         free(rf);
      }

      long r_qlen = slen - r_qstart;
      long l_qlen = l_qstart + 1;
      if (!r_qlen && !l_qlen) continue;
      long r_rlen = align_min((long)(r_qlen * (1 + alignopt.width_ratio)),
            index->size - r_rstart);
      long l_rlen = align_min((long)(l_qlen * (1 + alignopt.width_ratio)),
            r_qstart);

      // Added this line to avoid duplicate match when
      // unique seed align is not safe (while reseeding).
      if ((*seqmatches)->pos) {
         int64_t d = ((int64_t)(*seqmatches)->match[0].ref_s) - ((int64_t)l_rstart - l_rlen);
         if (d < slen && d > -slen) {
            if (VERBOSE_DEBUG) fprintf(stdout, "repeated alignment, skipping\n");
            continue;
         }
      }

      char * r_qry = read + r_qstart;
      char * l_qry = read + l_qstart;
      char * r_ref = index->genome + r_rstart;
      char * l_ref = index->genome + l_rstart;

      // Align forward and backward starting from (read_s, ref_s).
      long read_s, read_e, ref_s, ref_e;
      path_t align_r = (path_t){0,0,0}, align_l = (path_t){0,0,0};
      // Forward alignment (right).
      if(r_qlen) {
         align_r = align_bp(r_qlen, r_qry, r_rlen, r_ref, r_min,
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
         align_l = align_bp(l_qlen, l_qry, l_rlen, l_ref, l_min,
               ALIGN_BACKWARD, ALIGN_BACKWARD, alignopt);
         read_s =  l_qstart - align_l.row;
         ref_s = l_rstart - align_l.col;
      }
      else {
         read_s = l_qstart + 1;
         ref_s  = l_rstart + 1;
      }
      //      int align_span = align_r.row + align_l.row;
      int align_span = read_e - read_s + 1;
      int identities = align_span - align_r.score - align_l.score;
      int ref_min = (*seqmatches)->match[0].ref_s - align_span;
      int ref_max = ref_min + 2*align_span;

      // DEBUG.
      if (VERBOSE_DEBUG) fprintf(stdout, "alignment: [%ld<--<%ld,%d>--%ld][%ld--<%ld,%d>-->%ld]\n",read_s, l_qlen, align_l.score,l_qstart,r_qstart,r_qlen,align_r.score,read_e);

      // Filter out low indentity alignments.
      //      if (identities*1.0/slen < 0.85) continue;
      // Hit.
      //      if (identities > best) {
      if (identities > best && identities*1.0/slen >= 0.85) {
         match_t hit;
         // Update scores.
         if (best > 0) hit.interval = 1;
         if (s_hits != best) s_cnt = 1;
         else s_cnt++;
         
         s_hits = best;
         best = identities;
         // Fill/Update hit.
         hit.ref_e  = ref_e;
         hit.ref_s  = ref_s;
         hit.read_e = read_e;
         hit.read_s = read_s;

         // Undo reverse complement.
         int64_t ref_start = ref_s;
         if (ref_s >= index->size/2) {
            ref_start = index->size - 2 - ref_e;
         }
         // TEST DEBUG FOR FEATURE ONLY.
         hit.e_exp[0] = hit.e_exp[1] = hit.e_exp[2] = 0;
         int max_dist = 0, max_cnt = 0;
         // Gather annotation scores.
         for (int a = 0; a < index->ann->count; a++) {
            ann_t ann = index->ann->ann[a];
            int kmers = ref_e - ref_s + 1 - ann.k + 1;
            int dis = 0, cnt = 0, span = 0, last = -ann.k;
            for (int n = 0; n < kmers; n++) {
               int log10cnt;
               int d = ann_read(ann, ref_start + n, &log10cnt); 
               if (d > dis) {
                  dis = d;
                  cnt = log10cnt;
                  span = ann.k;
                  last = n;
               } else if (d == dis) {
                  span += min(n - last, ann.k);
                  last = n;
               }
            }
            hit.e_exp[a] = dis*(span*1.0/ann.k);
            if (hit.e_exp[a] > max_dist) {
               max_dist = hit.e_exp[a];
               max_cnt = cnt;
            }
            if (VERBOSE_DEBUG) fprintf(stdout, "\nannotation (%d,%d): d=%d, cnt=%d, max_dist=(%d,%d)\n",ann.k, ann.d, dis, cnt,max_dist,max_cnt);
         }
         // Best and second best identities.
         if (s_hits < identities - max_dist) {
            s_hits = identities-max_dist;
            s_cnt = pow(10,max_cnt);
            hit.interval = 0;
         }
         hit.hits = best;
         hit.s_cnt = s_cnt;
         hit.s_hits = s_hits;
         hit.annotation = max_dist;
         if (VERBOSE_DEBUG) fprintf(stdout, "set BEST=%d (%.2f) [second=%d (%.2f), cnt=%d]\n", best, best*1.0/slen, s_hits, s_hits*1.0/slen,s_cnt);
         // Alignment flag.
         if ((*seqmatches)->pos == 1)
            hit.interval = (*seqmatches)->match[0].interval;
         else {
            hit.interval = 0;
            (*seqmatches)->pos = 1;
         }
         // Set best match.
         (*seqmatches)->match[0] = hit;

      } else if ((*seqmatches)->pos && identities >= s_hits && (ref_min >= ref_s || ref_max <= ref_s)) {
         if (identities == s_hits && (*seqmatches)->match[0].interval) {
            s_cnt++;
            (*seqmatches)->match[0].s_cnt = s_cnt;
         } else {
            s_hits = identities;
            s_cnt = 1;
            (*seqmatches)->match[0].s_hits = identities;
            (*seqmatches)->match[0].s_cnt = 1;
            (*seqmatches)->match[0].interval = 1;
         }
         if (VERBOSE_DEBUG) fprintf(stdout, "SECOND=%d (%.2f) cnt=%d [best=%d (%.2f)]\n", (int)s_hits, s_hits*1.0/slen, s_cnt, (int)best, best*1.0/slen);

      }
   }
   return 0;
}

int
align_seeds
(
 char         * read,
 seedstack_t  * seeds,
 seedstack_t  * rseeds,
 matchlist_t ** seqmatches,
 index_t      * index,
 filteropt_t    opt,
 alignopt_t     alignopt
)
{
   int slen = strlen(read);

   // Convert seeds to hits.
   uint64_t hit_count, rev_count;
   // Forward seeds.
   hit_t * hits = compute_hits(seeds, index, &hit_count);
   // Reverse seeds.
   if (rseeds != NULL && rseeds->pos > 0) {
      hit_t * rhits = compute_hits(rseeds, index, &rev_count);
      // Revert hits.
      for (int i = 0; i < rev_count; i++) {
         rhits[i].locus = (index->size - 1) - (rhits[i].locus + rhits[i].depth);
         rhits[i].qrypos= slen - (rhits[i].qrypos + rhits[i].depth);
      }
      // Merge F/R hits.
      hits = realloc(hits, (hit_count + rev_count)*sizeof(hit_t));
      memcpy(hits+hit_count, rhits, rev_count*sizeof(hit_t));
      hit_count += rev_count;
      free(rhits);
   }
   int rval = align_hits(read,hits,hit_count,seqmatches,index,opt,alignopt);
   free(hits);
   return rval;
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
