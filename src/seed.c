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
   ann_t * ann = ann_find(k,index->ann);
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
             int ann_d = ann_read(*ann, loc, NULL);
            if (ann_d > d) {
               seed_t seed = (seed_t) {.errors = 0, .qry_pos = i, .ref_pos = pos};
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
   if (beg > slen-k)
      return (hit_t) {.errors = 1};
   // Find annotation table with same k.
   ann_t * ann = ann_find(k, index->ann);
   // Query buffer.
   uint8_t * q = malloc(k);
   for (int i = beg; i <= slen - k; i++) {
      if (check_kmer(query+i, q, k)) {
         values[i] = 4;
         continue;
      }
      uint64_t key = XXH64(q, k, 0);
      int value = htable_get(sht.htable, key);
      values[i] = value;
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
            int ann_d = ann_read(*ann, loc, NULL);
            if (ann_d > d) {
               free(q);
               return (hit_t) {.locus = locus, .qrypos = i, .depth = k, .errors = 0};
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
   return (hit_t) {.errors = 1};
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
               seed_t seed = (seed_t) {.errors = 1, .qry_pos = qrypos, .ref_pos = tmp};
               seedstack_push(seed, stack);
            }
         }
      }
      pos = next;
   }
   if (pos.ep >= pos.sp) {
      seed_t seed = (seed_t) {.errors = 0, .qry_pos = qrypos, .ref_pos = pos};
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
               seed_t seed = (seed_t) {.errors = 1, .qry_pos = qrypos, .ref_pos = refpos};
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
               seed_t seed = (seed_t) {.errors = 1, .qry_pos = qrypos, .ref_pos = refpos};
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
         hits[l++] = (hit_t) {.locus=loc, .qrypos=seed.qry_pos, .depth=seed.ref_pos.depth, .errors=seed.errors};
      }
   }
   
   return hits;

}

hit_t *
chain_seeds
(
 seedstack_t * seeds,
 int           seqlen,
 index_t     * index,
 int           dist_accept,
 double        read_ref_ratio,
 size_t      * hit_cnt
)
{
   // Convert seeds to hits.
   uint64_t hit_count;
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

   // Merge-sort hits.
   qsort(hits, hit_count, sizeof(hit_t), compar_hit_locus_qsort);

   if (VERBOSE_DEBUG) fprintf(stdout, "** CHAIN SEEDS (total %ld)\n", hit_count);

   // Aux variables.
   size_t i = 0;
   size_t cnt = 0;

   while (i < hit_count) {
      // Match start.
      int64_t span     = hits[i].depth;
      int64_t lbeg_tmp = hits[i].locus;
      int64_t lend_tmp = lbeg_tmp + span;
      int64_t rbeg_tmp = hits[i].qrypos;
      int64_t rend_tmp = rbeg_tmp + span;
      int64_t lbeg = lbeg_tmp, rbeg = rbeg_tmp, rend = rend_tmp;
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
         fprintf(stdout, "%s,%ld%c/%ld,%ld/%ld,%ld",
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 rbeg_tmp, rend_tmp-1,
                 span, span
                 );
      }

      // Extend match.
      while (i < hit_count) {
         // Compare genome and read distances.
         int64_t g_dist = (int64_t)(hits[i].locus) - (int64_t)(hits[i-1].locus);
         int64_t r_dist = (int64_t)(hits[i].qrypos) - (int64_t)(hits[i-1].qrypos);
         if (r_dist >= 0 && (g_dist <= (1+read_ref_ratio) * r_dist) && (g_dist >= (1-read_ref_ratio) * r_dist)) {
            uint64_t new_rend = hits[i].qrypos + hits[i].depth;
            if (new_rend > rend_tmp) {
               // Update hits.
               if (hits[i].qrypos <= rend_tmp) {
                  span += new_rend - rend_tmp;
                  rend_tmp = new_rend;
               }
               else {
                  span += hits[i].depth;
                  rbeg_tmp = hits[i].qrypos;
                  rend_tmp = new_rend;
                  lbeg_tmp = hits[i].locus;
               }
               // Check if this is the longest seed.
               if (rend_tmp - rbeg_tmp > rend - rbeg) {
                  rbeg = rbeg_tmp;
                  rend = rend_tmp;
                  lbeg = lbeg_tmp;
               }
            }

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
               fprintf(stdout, "\t%s,%ld%c/%d,%d/%d,%ld",
                       index->chr->name[chrnum],
                       g_start - index->chr->start[chrnum]+1,
                       dir ? '-' : '+',
                       hits[i].qrypos, hits[i].qrypos + hits[i].depth-1,
                       hits[i].depth, span
                       );
            }

            i++;
         } else break;
      }

      hit_t h = (hit_t) {.locus = lbeg, .qrypos = rbeg, .depth = rend - rbeg, .errors = span};
      hits[cnt++] = h;

      if (VERBOSE_DEBUG) fprintf(stdout, "\n");
   }
   
   // Realloc hit array and set hit count.
   hits = realloc(hits, cnt*sizeof(hit_t));
   *hit_cnt = cnt;

   // Sort by hit.errors (hit count)
   qsort(hits, cnt, sizeof(hit_t), compar_hit_errors_qsort);

   // Return hits.
   return hits;
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
compar_hit_errors_qsort
(
 const void * a,
 const void * b
)
{
   hit_t * ha = (hit_t *)a;
   hit_t * hb = (hit_t *)b;
   if (ha->errors < hb->errors) return 1;
   else return -1;
}

int
mem_unique
(
 seed_t    seed,
 index_t * index,
 hit_t   * hit
)
{
   // Seed must be single loci.
   if (seed.ref_pos.sp != seed.ref_pos.ep)
      return 0;
   int k = seed.ref_pos.depth;
   uint64_t locus = get_sa(seed.ref_pos.sp, index->sar);
   uint64_t loc = locus;
   if (loc >= index->size/2) loc = index->size - 1 - loc - k;
   // Select annotation.
   ann_t * ann = ann_find(k, index->ann);
   if (ann == NULL) return 0;
   // Check neighbors for all k-mers.
   int max_d = 0;
   for (int i = 0; i <= k - ann->k; i++) {
      max_d = max(max_d, ann_read(*ann, loc + i, NULL));
   }
   if (max_d > 1) {
      *hit = (hit_t) {.locus = locus, .qrypos = seed.qry_pos, .depth = k, .errors = 0};
      return 1;
   }
   return 0;

}

int
next_mem
(
 uint8_t * q,
 int       beg,
 int       end,
 seed_t  * mem,
 index_t * index
)
{
   int pos = beg;
   fmdpos_t bwpos;
   pos = extend_lr(q, pos, end, &bwpos, index);
   bwpos_t ref_pos = {.sp = bwpos.fp, .ep = bwpos.fp + bwpos.sz - 1, .depth = bwpos.dp};
   *mem = (seed_t){.errors = 0, .qry_pos = pos-bwpos.dp, .ref_pos = ref_pos};
   return pos;
}

int
seed_mem
(
 uint8_t      * q,
 int            beg,
 int            end,
 index_t      * index,
 seedstack_t ** stack
)
{
   int pos = beg;
   while (pos <= end) {
      fmdpos_t bwpos;
      pos = extend_lr(q, pos, end, &bwpos, index);
      bwpos_t ref_pos = {.sp = bwpos.fp, .ep = bwpos.fp + bwpos.sz - 1, .depth = bwpos.dp};
      seedstack_push((seed_t){.errors = 0, .qry_pos = pos-bwpos.dp, .ref_pos = ref_pos},stack);
   }
   return 0;
}

int
extend_lr
(
 uint8_t  * q,
 int        beg,
 int        end,
 fmdpos_t * bwpos,
 index_t  * index
)
{
   fmdpos_t next = *bwpos = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0};
   // Extend max to the left, starting from p.
   int i = beg;
   while (i >= 0 && q[i] < 4) {
      next = extend_bw(q[i--], next, index->bwt);
      if (next.sz) *bwpos = next;
      else break;
   }
   // Extend now max to the right.
   next = *bwpos;
   i = beg+1;
   while (i <= end && q[i] < 4) {
      next = extend_fw(q[i++], next, index->bwt);
      if (next.sz) *bwpos = next;
      else { i--; break; }
   }
   // Return next R position.
   return i + (i <= end ? (q[i] > 3) : 0);
}

/*
** This version of the function uses suffix_extend attempting
** to do forward search. Either use reverse complement and
** then correct the locus (which requires querying the SA),
** or start from the end (results may differ from BWA).
int
seed_interv
(
 uint8_t      * q,
 int            beg,
 int            end,
 int            min_len,
 int            max_loci,
 index_t      * index,
 seedstack_t ** stack
)
{
   if (min_len < 1) min_len = 1;
   while (beg <= end - min_len) {
      bwpos_t pos = index->bwt->bwt_base;
      while (pos.ep - pos.sp >= max_loci || pos.depth < min_len) {
         if (beg >= end) return 0;
         int c = q[beg++];
         if (c > 3) {
            if (beg > end - min_len) return 0;
            pos = index->bwt->bwt_base;
            continue;
         }
         suffix_extend(c,pos,&pos,index->bwt);
         fprintf(stdout, "pos[%d]: depth=%d, loci=%ld\n", beg-1, pos.depth, pos.ep-pos.sp+1);
      }
      if (pos.ep >= pos.sp) 
         seedstack_push((seed_t){.errors = 0, .qry_pos = beg-pos.depth, .ref_pos = pos}, stack);
   }
   return 0;
} 
*/ 

int
seed_interv
(
 uint8_t      * q,
 int            beg,
 int            end,
 int            min_len,
 int            max_loci,
 index_t      * index,
 seedstack_t ** stack
)
{
   if (min_len < 2) min_len = 2;
   while (beg <= end - min_len) {
      fmdpos_t pos = index->bwt->fmd_base;
      while (pos.sz > max_loci || pos.dp < min_len) {
         if (beg >= end) return 0;
         int c = q[beg++];
         if (c > 3) {
            pos = index->bwt->fmd_base;
            continue;
         }
         pos = extend_fw(c,pos,index->bwt);
      }
      if (pos.sz) 
         seedstack_push((seed_t){.errors = 0, .qry_pos = beg-pos.dp, .ref_pos = 
                  (bwpos_t){.sp=pos.fp, .ep=pos.fp+pos.sz-1, .depth=pos.dp}}, stack);
   }
   if (beg <= end) {
      fmdpos_t pos = index->bwt->fmd_base;
      // Extend min_len nt from the end.
      int i = end;
      while (pos.sz > max_loci || pos.dp < min_len) {
         pos = extend_bw(q[i--], pos, index->bwt);
      }
      if (pos.sz) {
         seedstack_push((seed_t){.errors = 0, .qry_pos = end-pos.dp+1, .ref_pos = 
                  (bwpos_t){.sp=pos.fp, .ep=pos.fp+pos.sz-1, .depth=pos.dp}}, stack);
      }
   }
   return 0;
} 



int
reseed_mem
(
 uint8_t      * q,
 seed_t         mem,
 int            slen,
 int            min_len,
 index_t      * index,
 seedstack_t ** stack
)
{
   /*
   int beg = mem.qry_pos;
   int end = mem.qry_pos + mem.ref_pos.depth;
   int c = (beg + end) / 2;
   */
   int c = mem.qry_pos + mem.ref_pos.depth/2;
   int mem_loc = mem.ref_pos.ep - mem.ref_pos.sp + 1;
   int cnt = 0;
   fmdpos_t * tmp = malloc((slen + 1)*sizeof(fmdpos_t));
   fmdpos_t next, cur = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0};
   // Query 1st base.
   cur = extend_fw(q[c],cur,index->bwt);
   // Extend forward.
   for (int i = c + 1; i < slen && cur.sz > mem_loc; i++) {
      next = extend_fw(q[i],cur,index->bwt);
      if (cur.sz > next.sz) {
         tmp[cnt++] = cur;
      }
      cur = next;
   }
   if (cur.sz > mem_loc) tmp[cnt++] = cur;
   // Resume and extend backward.
   int last = c;
   for (int t = cnt-1; t >= 0; t--) {
      fmdpos_t cur = tmp[t];
      if (cur.dp + c < min_len) break;
      int i = c-1;
      for (; i >= 0; i--) {
         next = extend_bw(q[i], cur, index->bwt);
         if (next.sz <= mem_loc) {
            if (cur.dp >= min_len && i + 1 < last) {
               last = i + 1;
               seedstack_push((seed_t){.errors = 0, .qry_pos = i+1, .ref_pos = (bwpos_t){.sp = cur.fp, .ep = cur.fp + cur.sz - 1, .depth = cur.dp}}, stack);
            }
            break;
         } else if (i == 0 && next.dp >= min_len && i < last) {
            last = i;
            seedstack_push((seed_t){.errors = 0, .qry_pos = i, .ref_pos = (bwpos_t){.sp = next.fp, .ep = next.fp + next.sz - 1, .depth = next.dp}}, stack);
            break;
         }
         cur = next;
      }
   }
   free(tmp);
   return 0;
}

int
seed_thr
(
 uint8_t      * q,
 int            beg,
 int            end,
 int            thr,
 index_t      * index,
 seedstack_t ** stack
)
{
   int pos = beg;
   while (pos <= end) {
      fmdpos_t bwpos;
      if (extend_lr_thr(q, &pos, end, thr, &bwpos, index)) {
         bwpos_t ref_pos = {.sp = bwpos.fp, .ep = bwpos.fp + bwpos.sz - 1, .depth = bwpos.dp};
         seedstack_push((seed_t){.errors = 0, .qry_pos = pos-bwpos.dp, .ref_pos = ref_pos},stack);
      }
   }
   return 0;
}

int
extend_lr_thr
(
 uint8_t  * q,
 int      * p,
 int        qlen,
 int        thr,
 fmdpos_t * bwpos,
 index_t  * index
)
{
   fmdpos_t cur = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0};
   fmdpos_t next;
   int unset = 1;
   // Extend max to the left, starting from p.
   int i = *p;
   while (i >= 0 && q[i] < 4) {
      next = extend_bw(q[i--], cur, index->bwt);
      if (unset && cur.sz <= thr) {
         unset = 0;
         *bwpos = cur;
      }
      if (next.sz == 0) break;
      cur = next;
   }
   // Extend now max to the right.
   i = *p + 1;
   while (i < qlen && q[i] < 4) {
      next = extend_fw(q[i++], cur, index->bwt);
      if (unset && cur.sz <= thr) {
         unset = 0;
         *bwpos = cur;
      }
      if (next.sz == 0) {i--; break;}
      cur = next;
   }
   *p = i + (i < qlen ? (q[i] > 3) : 0);
   // Return next R position.
   return !unset;
}   

int
seed_mem_bp
(
 uint8_t      * q,
 int            beg,
 int            end,
 int            bp,
 index_t      * index,
 seedstack_t ** stack
)
{
   int pos = beg;
   while (pos <= end) {
      fmdpos_t bwpos;
      pos = extend_lr_bp(q, pos, end, bp, &bwpos, index);
      if (bwpos.dp) {
         bwpos_t ref_pos = {.sp = bwpos.fp, .ep = bwpos.fp + bwpos.sz - 1, .depth = bwpos.dp};
         seedstack_push((seed_t){.errors = 0, .qry_pos = pos-bwpos.dp, .ref_pos = ref_pos},stack);
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
      if (next.sz >= bp) *bwpos = next;
      else break;
   }
   // Extend now max to the right.
   next = *bwpos;
   i = p+1;
   while (i < qlen && q[i] < 4) {
      next = extend_fw(q[i++], next, index->bwt);
      if (next.sz >= bp) *bwpos = next;
      else {i--; break;}
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
   // Seeds are sorted by start position, ascending.
   seed_t cur = mem->seed[0];
   int    len = cur.ref_pos.depth;
   int    end = cur.qry_pos + cur.ref_pos.depth;
   // Find longest non-overlapping seeds.
   for (int i = 1; i < mem->pos; i++) {
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
