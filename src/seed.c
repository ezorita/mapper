#include "seed.h"

int
check_kmer
(
 uint8_t * seq,
 uint8_t * query,
 int       slen
)
{
   int strand = 0;
   for (int i = 0; i < slen; i++) {
      // Check 'N'
      if (seq[i] > 3) return -1;
      // Compare forward and reverse strands.
      if (!strand) {
         if (3 - seq[slen-1-i] < seq[i]) strand = -1;
         else if (3 - seq[slen-1-i] > seq[i]) strand = 1;
      }
   }
   if (strand < 0)
      for (int i = 0; i < slen; i++) query[i] = 3-seq[slen-i-1];
   else
      memcpy(query, seq, slen);
   return 0;
}

int
find_uniq_seeds
(
 int            slen,
 uint8_t      * query,
 uint8_t      * values,
 index_t      * index,
 seedstack_t ** seeds
)
{
   // Load seed table.
   sht_t sht = index->sht->sht[0];
   int k = sht.k;
   int d = sht.d;
   int r = sht.repeat_thr;
   int z = 0;
   int zp = -k;
   ann_t ann;
   ann_find(k,index->ann,&ann);
   // Query buffer.
   uint8_t * q = malloc(k);
   for (int i = 0; i <= slen - k; i++) {
      if (check_kmer(query+i, q, k)) {
         values[i] = -1;
         continue;
      }
      uint64_t key = XXH64(q, k, 0);
      uint8_t value = htable_get(sht.htable, key);
      values[i] = value;

      //DEBUG
      if (VERBOSE_DEBUG) {
         pathstack_t * ps = pathstack_new(64);
         blocksearch(q, k, d, index, &ps);
         int count = 0;
         for (int c = 0; c < ps->pos; c++) {
            count += ps->path[c].pos.sz;
         }
         fprintf(stdout, "[%d] ann:%d, loci:%d\n", i, value,count);
         free(ps);
      }

      // Mask out mismatched seeds.
      if (value > 0) {
         if (i-zp < k) {
            values[i] = 0;
            continue;
         }
         z = 1;
      } else if (z) {
         z = 0;
         zp = i;
      }
      // Check value
      if (value == 1) {
         // Mask out mismatched seeds.
         if (i-zp < k) {
            values[i] = 0;
            continue;
         }
         // Set last existing to current.
         z = i;
         // Find seed in index.
         bwpos_t pos = index->bwt->bwt_base;
         for (int j = i+k-1; j >= i && pos.ep >= pos.sp; j--) suffix_extend(query[j],pos,&pos,index->bwt);
         if (pos.ep == pos.sp) {
            uint64_t locus = get_sa(pos.sp, index->sar);
            uint64_t loc = locus;
            if (loc >= index->size/2) loc = index->size - 1 - loc - k;
             int ann_d = ann_read(ann, loc, NULL);
            if (ann_d > d) {
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = i, .ref_pos = pos};
               seedstack_push(seed, seeds);
            } else {
               values[i] = 2;
            }
         } else {
            values[i] = (pos.ep < pos.sp ? 0 : (pos.ep - pos.sp < r ? 2 : 3));
         }
      }
   }
   free(q);
   return 0;
}

hit_t
find_uniq_seed
(
 int            beg,
 int            slen,
 uint8_t      * query,
 uint8_t      * values,
 index_t      * index
)
{
   // Load seed table.
   sht_t sht = index->sht->sht[0];
   int     k = sht.k;
   int     d = sht.d;
   int     r = sht.repeat_thr;
   int     z = 0;
   int    zp = -k;
   // Find annotation table with same k.
   ann_t ann;
   ann_find(k, index->ann, &ann);
   // Query buffer.
   uint8_t * q = malloc(k);
   for (int i = beg; i <= slen - k; i++) {
      if (check_kmer(query+i, q, k)) {
         values[i] = -1;
         continue;
      }
      uint64_t key = XXH64(q, k, 0);
      int value = htable_get(sht.htable, key);
      values[i] = value;
      // DEBUG.
      /*
      if (VERBOSE_DEBUG) {
         fprintf(stdout, "[%d] ", i);
         for (int n = 0; n < k; n++) fprintf(stdout, "%c", bases[q[n]]);
         fprintf(stdout, "\t%d", value);
         pathstack_t * ps = pathstack_new(64);
         blocksearch(q, k, d, index, &ps);
         int count = 0;
         for (int c = 0; c < ps->pos; c++) {
            count += ps->path[c].pos.sz;
         }
         fprintf(stdout, "\t(%d)\t",count);
         for (int c = 0; c < ps->pos; c++) {
            uint64_t l = get_sa(ps->path[c].pos.fp, index->sar);
            fprintf(stdout, "[p:%ld,s:%d]", l, ps->path[c].score);
            for (int n = 0; n < k; n++)
               fprintf(stdout, "%c", index->genome[l+n]);
            fprintf(stdout, "\t");
         }
         fprintf(stdout, "\n");
         free(ps);
      }
      */
      if (value > 0) {
         // Mask out mismatched seeds.
         if (i-zp < k) {
            values[i] = 0;
            continue;
         }
         z = 1;
      }
      if (value == 1) {
         // Find seed in genome.
         bwpos_t pos = index->bwt->bwt_base;
         for (int j = i+k-1; j >= i && pos.ep >= pos.sp; j--) suffix_extend(query[j],pos,&pos,index->bwt);
         if (pos.ep == pos.sp) {
            uint64_t locus = get_sa(pos.sp, index->sar);
            uint64_t loc = locus;
            if (loc >= index->size/2) loc = index->size - 1 - loc - k;
            int ann_d = ann_read(ann, loc, NULL);
            if (ann_d > d) {
               free(q);
               return (hit_t) {.locus = locus, .qrypos = i, .depth = k, .bulk = 0};
            } else {
               values[i] = 2;
            }
         } else {
            values[i] = (pos.ep < pos.sp ? 0 : (pos.ep - pos.sp < r ? 2 : 3));
         }
      } else if (value == 0 && z) {
         zp = i;
         z = 0;
      } 
   }
   free(q);
   return (hit_t) {.bulk = 1};
}

int
find_thr_seed
(
 int            beg,
 int            slen,
 int            k,
 uint8_t      * query,
 index_t      * index
)
{
   sht_t sht = index->sht->sht[0];
   uint8_t * q = malloc(k);
   for (int i = beg; i <= slen - k; i++) {
      if (check_kmer(query+i, q, k)) continue;
      uint64_t key = XXH64(q, k, 0);
      // Test.
      int value = htable_get(sht.htable, key);
      if (value == 0 || value == 2) {
         free(q);
         return i;
      }
   }
   free(q);
   return -1;
}


int
count_mismatch
(
 char         * seq,
 size_t         slen,
 index_t      * index
)
{
   int count = 0;
   bwpos_t pos = {.sp = 0, .ep = index->size-1, .depth = 0};
   for (int i = slen-1; i >= 0; i--) {
      int nt = translate[(uint8_t)seq[i]];
      bwpos_t next;
      for (int j = 0; j < 4; j++) {
          if (j == nt) {
            suffix_extend(j, pos, &next, index->bwt);
         } else {
            bwpos_t tmp;
            suffix_extend(j, pos, &tmp, index->bwt);
            for (int k = i-1; k >= 0 && tmp.ep >= tmp.sp; k--) {
               suffix_extend(translate[(uint8_t)seq[k]], tmp, &tmp, index->bwt);
            }
            if (tmp.ep >= tmp.sp) {
               count += tmp.ep - tmp.sp + 1;
            }
         }
      }
      pos = next;
   }
   if (pos.ep >= pos.sp) count += pos.ep - pos.sp + 1;

   return count;
}

int
find_mismatch
(
 char         * seq,
 size_t         slen,
 int            qrypos,
 seedstack_t ** stack,
 index_t      * index
)
{
   bwpos_t pos = {.sp = 0, .ep = index->size-1, .depth = 0};
   for (int i = slen-1; i >= 0; i--) {
      int nt = translate[(uint8_t)seq[i]];
      bwpos_t next;
      for (int j = 0; j < 4; j++) {
          if (j == nt) {
            suffix_extend(j, pos, &next, index->bwt);
         } else {
            bwpos_t tmp;
            suffix_extend(j, pos, &tmp, index->bwt);
            for (int k = i-1; k >= 0 && tmp.ep >= tmp.sp; k--) {
               suffix_extend(translate[(uint8_t)seq[k]], tmp, &tmp, index->bwt);
            }
            if (tmp.ep >= tmp.sp) {
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = tmp};
               seedstack_push(seed, stack);
            }
         }
      }
      pos = next;
   }
   if (pos.ep >= pos.sp) {
      seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = pos};
      seedstack_push(seed, stack);
   }

   return 0;
}

int
seed_mismatch
(
 char         * seq,
 size_t         slen,
 int            qrypos,
 seedstack_t ** stack,
 index_t      * index
)
{
   int c = slen/2;
   fmdpos_t pos = {0,0,index->size};
   // forward extend.
   //  forward extend exact until c-1.
   int i = 0;
   for (; i < c; i++) {
      pos = extend_fw(translate[(int)seq[i]],pos,index->bwt);
   }
   //  mismatch all the next positions.
   for (; i < slen; i++) {
      int nt = translate[(int)seq[i]];
      fmdpos_t newpos[5];
      extend_fw_all(pos, newpos, index->bwt);
      for (int j = 0; j < 4; j++) {
         if (j == nt) {
            // Extend perfect seed.
            pos = newpos[j];
         } else {
            // Extend mismatched seed.
            for (int k = i+1; k < slen && newpos[j].sz; k++)
               newpos[j] = extend_fw(translate[(int)seq[k]],newpos[j],index->bwt);
            if (newpos[j].sz) {
               // Save seed.
               bwpos_t refpos = {.sp = newpos[j].fp, .ep = newpos[j].fp + newpos[j].sz - 1, .depth = slen};
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = refpos};
               seedstack_push(seed, stack);
            }
         }
      }
   }
   // reverse extend.
   pos = (fmdpos_t){0,0,index->size};
   i = slen-1;
   for (; i >= c; i--) {
      pos = extend_bw(translate[(int)seq[i]],pos,index->bwt);
   }
   //  mismatch all the next positions.
   for (; i >= 0; i--) {
      int nt = translate[(int)seq[i]];
      fmdpos_t newpos[5];
      extend_bw_all(pos, newpos, index->bwt);
      for (int j = 0; j < 4; j++) {
         if (j == nt) {
            // Extend perfect seed.
            pos = newpos[j];
         } else {
            // Extend mismatched seed.
            for (int k = i-1; k >= 0 && newpos[j].sz; k++)
               newpos[j] = extend_bw(translate[(int)seq[k]],newpos[j],index->bwt);
            if (newpos[j].sz) {
               // Save seed.
               bwpos_t refpos = {.sp = newpos[j].fp, .ep = newpos[j].fp + newpos[j].sz - 1, .depth = slen};
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = refpos};
               seedstack_push(seed, stack);
            }
         }
      }
   }
   return 0;
}



hit_t *
compute_hits
(
 seedstack_t * seeds,
 index_t     * index,
 uint64_t    * hit_cnt
)
{
   // Count loci.
   uint64_t cnt = 0;
   for (size_t i = 0; i < seeds->pos; i++) cnt += seeds->seed[i].ref_pos.ep - seeds->seed[i].ref_pos.sp + 1;
   *hit_cnt = cnt;

   // Find loci in Suffix Array.
   hit_t * hits = malloc(cnt * sizeof(hit_t));
   size_t l = 0;
   for (size_t i = 0; i < seeds->pos; i++) {
      seed_t seed = seeds->seed[i];
      for (size_t j = seed.ref_pos.sp; j <= seed.ref_pos.ep; j++) {
         int64_t loc = get_sa(j,index->sar);
         hits[l++] = (hit_t) {.locus = loc, .qrypos = seed.qry_pos, .depth = seed.ref_pos.depth, .bulk = seed.ref_pos.ep - seed.ref_pos.sp + 1};
      }
   }
   
   return hits;

}

int
chain_seeds
(
 seedstack_t * seeds,
 seedstack_t * rseeds,
 int           seqlen,
 matchlist_t * matchlist,
 index_t     * index,
 int           dist_accept,
 double        read_ref_ratio
)
{
   if (matchlist->size < 1) return -1;

   // Convert seeds to hits.
   uint64_t hit_count, rev_count;
   // Forward seeds.
   hit_t * hits = compute_hits(seeds, index, &hit_count);

   // DEBUG.
   if (VERBOSE_DEBUG) {
      for (int i = 0; i < hit_count; i++) {
         uint64_t g_start;
         int dir = 0;
         if (hits[i].locus >= index->size/2) {
            g_start = index->size - (hits[i].locus+hits[i].depth) - 2;
            dir = 1;
         } else {
            g_start = hits[i].locus + 1;
            dir = 0;
         }
         // Search chromosome name.
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "[%d,%d] (%d) %s,%ld%c\n",
                 hits[i].qrypos, hits[i].qrypos+hits[i].depth-1,
                 hits[i].depth,
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 dir ? '-' : '+'
                 );
      }
   }

   // Reverse seeds.
   if (rseeds->pos > 0) {
      hit_t * rhits = compute_hits(rseeds, index, &rev_count);
      // Revert hits.
      for (int i = 0; i < rev_count; i++) {
         rhits[i].locus = (index->size - 1) - (rhits[i].locus + rhits[i].depth);
         rhits[i].qrypos= seqlen - (rhits[i].qrypos + rhits[i].depth);
      }
      // Merge F/R hits.
      hits = realloc(hits, (hit_count + rev_count)*sizeof(hit_t));
      memcpy(hits+hit_count, rhits, rev_count*sizeof(hit_t));
      hit_count += rev_count;
      free(rhits);
   }
   
   // Merge-sort hits.
   //   mergesort_mt(hits, hit_count, sizeof(hit_t), 0, 1, compar_hit_locus);
   qsort(hits, hit_count, sizeof(hit_t), compar_hit_locus_qsort);

   if (VERBOSE_DEBUG) fprintf(stdout, "** CHAIN SEEDS (total %ld)\n", hit_count);

   // Reset matchlist.
   matchlist->pos = 0;

   // Aux variables.
   long minv = 0;
   int  mini = 0;
   size_t i = 0;

   while (i < hit_count) {
      // Match start.
      int64_t span     = hits[i].depth;
      int64_t lbeg_tmp = hits[i].locus;
      int64_t lend_tmp = lbeg_tmp + span;
      int64_t rbeg_tmp = hits[i].qrypos;
      int64_t rend_tmp = rbeg_tmp + span;
      int     bulk     = hits[i].bulk;
      int64_t lbeg, lend, rbeg, rend;
      lbeg = lend = rbeg = rend = 0;
      i++;

      // Debug seeds...
      if (VERBOSE_DEBUG) {
         uint64_t g_start;
         int dir = 0;
         if (lbeg_tmp >= index->size/2) {
            g_start = index->size - lend_tmp - 2;
            dir = 1;
         } else {
            g_start = lbeg_tmp + 1;
            dir = 0;
         }
         // Search chromosome name.
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "%s,%ld%c/%ld,%ld/%ld,%ld,%d",
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 rbeg_tmp, rend_tmp-1,
                 span, span, bulk
                 );
      }

      // Extend match.
      while (i < hit_count) {
         // Compare genome and read distances.
         int64_t g_dist = (int64_t)(hits[i].locus) - (int64_t)(hits[i-1].locus);
         int64_t r_dist = (int64_t)(hits[i].qrypos) - (int64_t)(hits[i-1].qrypos);
         if (r_dist >= 0 && (g_dist <= (1+read_ref_ratio) * r_dist) && (g_dist >= (1-read_ref_ratio) * r_dist)) {
            uint64_t new_lend = hits[i].locus + hits[i].depth;
            uint64_t new_rend = hits[i].qrypos + hits[i].depth;
            if (new_rend > rend_tmp) {
               // Update hits.
               if (hits[i].qrypos <= rend_tmp) {
                  span += new_rend - rend_tmp;
               }
               else {
                  span += hits[i].depth;
                  // Check if this is the biggest stretch.
                  if (rend_tmp - rbeg_tmp > rend - rbeg) {
                     rbeg = rbeg_tmp;
                     rend = rend_tmp;
                     lbeg = lbeg_tmp;
                     lend = lend_tmp;
                  }
                  // Follow new stretch.
                  rbeg_tmp = hits[i].qrypos;
                  lbeg_tmp = hits[i].locus;
               }
               // Update ends.
               rend_tmp = new_rend;
               lend_tmp = new_lend;

            }
            bulk *= hits[i].bulk;

            // Debug seed chain...
            if (VERBOSE_DEBUG) {
               uint64_t g_start;
               int dir = 0;
               if (hits[i].locus >= index->size/2) {
                  g_start = index->size - hits[i].locus - hits[i].depth - 2;
                  dir = 1;
               } else {
                  g_start = hits[i].locus + 1;
                  dir = 0;
               }
               // Search chromosome name.
               int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
               // Print results.
               fprintf(stdout, "\t%s,%ld%c/%d,%d/%d,%ld,%d",
                       index->chr->name[chrnum],
                       g_start - index->chr->start[chrnum]+1,
                       dir ? '-' : '+',
                       hits[i].qrypos, hits[i].qrypos + hits[i].depth-1,
                       hits[i].depth, span, bulk
                       );
            }

            i++;
         } else break;
      }
      // Check if last was the biggest stretch.
      if (rend_tmp - rbeg_tmp > rend - rbeg) {
         rbeg = rbeg_tmp;
         rend = rend_tmp;
         lbeg = lbeg_tmp;
         lend = lend_tmp;
      }


      if (VERBOSE_DEBUG) fprintf(stdout, "\n");

      // Store match.
      // Non-significant streaks (bulk) will not be saved.
      if (span > minv && bulk == 0) {
         // Allocate match.
         match_t match;
         match.ref_e  = lend - 1;
         match.read_e = rend - 1;
         match.ref_s  = lbeg;
         match.read_s = rbeg;
         match.hits   = span;
         match.s_hits = 0;
         match.s_cnt  = 0;

         // Debug match
         if (VERBOSE_DEBUG) {
            uint64_t g_start;
            int dir = 0;
            if (match.ref_s >= index->size/2) {
               g_start = index->size - match.ref_e - 2;
               dir = 1;
            } else {
               g_start = match.ref_s + 1;
               dir = 0;
            }
            // Search chromosome name.
            int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
            // Print results.
            fprintf(stdout, "[%d]\t%s,%ld%c\t%ld\t%d,%d\n",
                    match.hits,
                    index->chr->name[chrnum],
                    g_start - index->chr->start[chrnum]+1,
                    dir ? '-' : '+',
                    match.ref_s,
                    match.read_s, match.read_e
                    );
         }

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
         }
         else {
            if(VERBOSE_DEBUG) {
               match_t m = matchlist->match[mini];
               fprintf(stdout, "replaces [%d]\t%ld\t%d,%d\n",m.hits,m.ref_s,m.read_s,m.read_e);
            }
            matchlist->match[mini] = match;
         }
               
         // Find the minimum that will be substituted next.
         if (matchlist->pos == matchlist->size) {
            mini = 0;
            minv = matchlist->match[mini].hits;
            for (int j = 1 ; j < matchlist->size; j++) {
               if (matchlist->match[j].hits < minv) {
                  mini  = j;
                  minv = matchlist->match[j].hits;
               }
            }
         }
      }
   }

   // Free hits.
   free(hits);

   return 0;
}

int
compar_hit_locus_qsort
(
 const void * a,
 const void * b
)
{
   hit_t * ha = (hit_t *)a;
   hit_t * hb = (hit_t *)b;
   if (ha->locus > hb->locus) return 1;
   else if (ha->locus < hb->locus) return -1;
   else return (ha->depth > hb->depth ? 1 : -1);
}

int
seed_mem
(
 uint8_t      * q,
 int            slen,
 index_t      * index,
 seedstack_t ** stack
)
{
   int pos = 0;
   while (pos < slen) {
      fmdpos_t bwpos;
      pos = extend_lr(q, pos, slen, &bwpos, index);
      bwpos_t ref_pos = {.sp = bwpos.fp, .ep = bwpos.fp + bwpos.sz - 1, .depth = bwpos.dp};
      seedstack_push((seed_t){.bulk = 0, .qry_pos = pos-bwpos.dp, .ref_pos = ref_pos},stack);
   }
   return 0;
}

int
extend_lr
(
 uint8_t  * q,
 int        p,
 int        qlen,
 fmdpos_t * bwpos,
 index_t  * index
)
{
   fmdpos_t next = *bwpos = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0};
   // Extend max to the left, starting from p.
   int i = p;
   while (i >= 0 && q[i] < 4) {
      next = extend_bw(q[i--], next, index->bwt);
      if (next.sz) *bwpos = next;
      else break;
   }
   // Extend now max to the right.
   next = *bwpos;
   i = p+1;
   while (i < qlen && q[i] < 4) {
      next = extend_fw(q[i++], next, index->bwt);
      if (next.sz) *bwpos = next;
      else { i--; break; }
   }
   // Return next R position.
   return i + (i < qlen ? (q[i] > 3) : 0);
}

int
seed_mem_bp
(
 uint8_t      * q,
 int            slen,
 int            bp,
 index_t      * index,
 seedstack_t ** stack
)
{
   int pos = 0;
   while (pos < slen) {
      fmdpos_t bwpos;
      pos = extend_lr_bp(q, pos, slen, bp, &bwpos, index);
      if (bwpos.dp) {
         bwpos_t ref_pos = {.sp = bwpos.fp, .ep = bwpos.fp + bwpos.sz - 1, .depth = bwpos.dp};
         seedstack_push((seed_t){.bulk = 0, .qry_pos = pos-bwpos.dp, .ref_pos = ref_pos},stack);
      }
   }
   return 0;
}

int
extend_lr_bp
(
 uint8_t  * q,
 int        p,
 int        qlen,
 int        bp,
 fmdpos_t * bwpos,
 index_t  * index
)
{
   fmdpos_t next = *bwpos = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0};
   // Extend max to the left, starting from p.
   int i = p;
   while (i >= 0 && q[i] < 4) {
      next = extend_bw(q[i--], next, index->bwt);
      if (next.sz == 0) break;
      if (!bwpos->dp && next.sz <= bp) *bwpos = next;
   }
   // Extend now max to the right.
   next = *bwpos;
   i = p+1;
   while (i < qlen && q[i] < 4) {
      next = extend_fw(q[i++], next, index->bwt);
      if (next.sz == 0) { i--; break; }
      if (!bwpos->dp && next.sz <= bp) *bwpos = next;
   }
   // Return next R position.
   return i + (i < qlen ? (q[i] > 3) : 0);
}

int
get_smem
(
 seedstack_t ** smem,
 seedstack_t  * mem
)
{
   if (mem->pos == 0) return 0;
   // Seeds are sorted by start position, descending.
   seed_t cur = mem->seed[mem->pos-1];
   int    len = cur.ref_pos.depth;
   int    end = cur.qry_pos + cur.ref_pos.depth;
   // Find longest non-overlapping seeds.
   for (int i = mem->pos-2; i >= 0; i--) {
      seed_t s = mem->seed[i];
      if (s.qry_pos >= end) {
         seedstack_push(cur, smem);
         cur = s;
         len = cur.ref_pos.depth;
         end = cur.qry_pos + cur.ref_pos.depth;
      } else if (s.ref_pos.depth > len) {
         cur = s;
         len = cur.ref_pos.depth;
         end = cur.qry_pos + cur.ref_pos.depth;
      }
   }
   seedstack_push(cur, smem);

   return 0;
}


seedstack_t *
naive_smem
(
 char * seq,
 seedopt_t opt,
 index_t * index
)
{
   uint32_t slen = strlen(seq);
   int32_t qry_end = slen - 1;
   int32_t last_qry_pos = slen;
   // Create stack.
   seedstack_t * stack = seedstack_new(SEEDSTACK_SIZE);

   for (int i = slen-1; i > 0; i--) {
      // Start position.
      bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
      qry_end = i;
      int j;
      uint64_t new_loci;
      for (j = 0; j <= i; j++) {
         bwpos_t newpos;
         int nt = translate[(uint8_t)seq[i-j]];
         // Extend suffix (Backward search).
         suffix_extend(nt, pos, &newpos, index->bwt);
         // Count loci.
         new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);

         if (new_loci < opt.min_loci || nt == 4) {
            // Check previous suffix.
            uint64_t loci = (pos.ep < pos.sp ? 0 : pos.ep - pos.sp + 1);
            int32_t qry_pos = qry_end + 1 - pos.depth;
            if (loci <= opt.max_loci && pos.depth >= opt.min_len && qry_pos < last_qry_pos) {
               // Seed found.
               seed_t seed = (seed_t) {.qry_pos = qry_pos, .ref_pos = pos};
               seedstack_push(seed, &stack);
               last_qry_pos = qry_pos;
            }
            if (nt == 4) {
               i = i - j;
            }
            break;
         } else if (newpos.depth == opt.max_len || j == i) {
            int32_t qry_pos = qry_end + 1 - newpos.depth;
            if (new_loci <= opt.max_loci && qry_pos < last_qry_pos) {
               // Seed found.
               seed_t seed = (seed_t) {.qry_pos = qry_pos, .ref_pos = newpos};
               seedstack_push(seed, &stack);
               last_qry_pos = qry_pos;
            }
            break;
         }
         pos = newpos;
      }
   }
   return stack;
}

seed_t
longest_mem
(
 char * seq,
 index_t * index
)
{
   uint32_t slen = strlen(seq);
   int32_t mlen = 0;
   int32_t qry_end = slen - 1;
   seed_t seed;

   for (int i = slen-1; i > 0; i--) {
      // Start position.
      bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
      qry_end = i;
      int j;
      uint64_t new_loci;
      for (j = 0; j <= i; j++) {
         bwpos_t newpos;
         int nt = translate[(uint8_t)seq[i-j]];
         // Extend suffix (Backward search).
         suffix_extend(nt, pos, &newpos, index->bwt);
         // Count loci.
         new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);
         if (!new_loci || nt == 4) {
            // Check previous suffix.
            if (pos.depth > mlen) {
               mlen = pos.depth;
               // Store longest seed.
               seed = (seed_t) {.qry_pos = qry_end + 1 - newpos.depth, .ref_pos = pos};
            }
            if (nt == 4) {
               i = i - j;
            }
            break;
         } else if (j == i) {
            if (newpos.depth > mlen) {
               mlen = newpos.depth;
               // Store longest seed.
               seed = (seed_t) {.qry_pos = qry_end + 1 - newpos.depth, .ref_pos = newpos};
            }
            break;
         }
         pos = newpos;
      }
   }
   return seed;
}


seedstack_t *
seedstack_new
(
 int size
)
{
   seedstack_t * stack = malloc(sizeof(seedstack_t) + size*sizeof(seed_t));
   if (stack == NULL) return NULL;
   stack->pos = 0;
   stack->size = size;
   return stack;
}

int
seedstack_push
(
 seed_t         seed,
 seedstack_t ** stackp
)
{
   seedstack_t * stack = *stackp;
   // Realloc if full.
   if (stack->pos >= stack->size) {
      size_t newsize = stack->size * 2;
      stack = realloc(stack, sizeof(seedstack_t) + newsize * sizeof(seed_t));
      if (stack == NULL) return -1;
      *stackp = stack;
      stack->size = newsize;
   }
   // Push seed.
   stack->seed[stack->pos++] = seed;
   return 0;
}

int
seedstack_append
(
 seedstack_t ** dest,
 seedstack_t  * src
)
{
   if (src->pos < 1) return 0;

   seedstack_t * stack = *dest;
   // Realloc if full.
   size_t newpos = stack->pos + src->pos;
   if (newpos >= stack->size) {
      size_t newsize = newpos + 1;
      stack = realloc(stack, sizeof(seedstack_t) + newsize * sizeof(seed_t));
      if (stack == NULL) return -1;
      *dest = stack;
      stack->size = newsize;
   }
   // Copy seeds.
   memcpy(stack->seed + stack->pos, src->seed, src->pos * sizeof(seed_t));
   stack->pos = newpos;
   return 0;
}

int
compar_hit_locus
(
 const void * a,
 const void * b,
 const int param
)
{
   hit_t * ha = (hit_t *)a;
   hit_t * hb = (hit_t *)b;
   if (ha->locus > hb->locus) return 1;
   else if (ha->locus < hb->locus) return -1;
   else return (ha->depth > hb->depth ? 1 : -1);
}

int
compar_seed_start
(
 const void * a,
 const void * b,
 const int param
)
{
   seed_t * sa = (seed_t *)a;
   seed_t * sb = (seed_t *)b;
   if (sa->qry_pos > sb->qry_pos) return 1;
   else if (sa->qry_pos < sb->qry_pos) return -1;
   else return (sa->ref_pos.depth > sb->ref_pos.depth ? -1 : 1);
}
