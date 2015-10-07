#include "seed.h"

char *
rev_comp
(
 char * seq
)
{
   int64_t slen = strlen(seq);
   char * rseq = malloc(slen+1);
   for (int64_t i = 0; i < slen; i++)
      rseq[slen-1-i] = revcomp[(uint8_t)seq[i]];
   rseq[slen] = 0;
   return rseq;
}

int
merge_seeds
(
 seedstack_t  * fw,
 seedstack_t  * rv,
 seedstack_t ** m
)
{
   int64_t i = 0, j = rv->pos-1;
   while (i < fw->pos && j >= 0) {
      // A > B. (push A)
      if (rv->seed[j].qry_pos < fw->seed[i].qry_pos) {
         seedstack_push(fw->seed[i++], m);
         continue;
      }
      // A = B. (ignore seed)
      if (rv->seed[j].qry_pos == fw->seed[i].qry_pos) {
         seedstack_push(fw->seed[i], m);
         if (rv->seed[j].ref_pos.depth != -fw->seed[i].ref_pos.depth)
            seedstack_push(rv->seed[j], m);
         i++; j--;
         continue;
      }
      // A < B. (push B)
      while (rv->seed[j].qry_pos > fw->seed[i].qry_pos) {
         seedstack_push(rv->seed[j--], m);
         if (j < 0) break;
      }
   }
   // Push remaining seeds.
   while (i < fw->pos) seedstack_push(fw->seed[i++], m);
   while (j >= 0) seedstack_push(rv->seed[j--], m);

   return 0;
}

int
seed_thr
(
 char         * seq,
 size_t         slen,
 seedstack_t ** stack,
 seedopt_t      opt,
 index_t      * index
)
{
   int32_t i = slen-1, qry_end = slen-1;
   
   // Start position.
   bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
   while (i >= 0) {
      int nt = translate[(uint8_t)seq[i]];
      if (nt > 3) {
         pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
         i--;
         qry_end = i;
         continue;
      }
      bwpos_t newpos;
      // Extend suffix (Backward search).
      suffix_extend(nt, pos, &newpos, index);
      // Count loci.
      uint64_t new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);
      // If seed is broken, shrink previous seed.
      if (new_loci == 0) {
         int depth = pos.depth;
         suffix_shrink(pos, &newpos, index);
         qry_end -= depth - pos.depth;
      // Seed found, save and shrink suffix.
      } else if (new_loci <= opt.thr_seed && newpos.depth >= opt.min_len) {
         // Save seed.
         int32_t qry_pos = qry_end + 1 - newpos.depth;
         seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qry_pos, .ref_pos = newpos};
         seedstack_push(seed, stack);
         // Shrink suffix.
         int depth = newpos.depth;
         suffix_shrink(newpos, &newpos, index);
         qry_end -= depth - newpos.depth;

         i--;
      } else {
         i--;
      }
      pos = newpos;
   }

   return 0;
}

int
reseed_mem
(
 char         * seq,
 seedstack_t ** seeds,
 seedopt_t      opt,
 index_t      * index
)
{
   int n_seeds = (*seeds)->pos;
   for (size_t i = 0; i < n_seeds; i++) {
      // Do not reseed short seeds.
      seed_t seed = (*seeds)->seed[i];
      if (seed.ref_pos.depth < opt.reseed_len) continue;
      // Reseed, min length = seedlen/2, min loci = seedloci + 1.
      seedopt_t rseed_opt = opt;
      rseed_opt.min_len = seed.ref_pos.depth/2;
      rseed_opt.min_loci = seed.ref_pos.ep - seed.ref_pos.sp + 2;
      // Seed forward.
      int beg = seed.qry_pos;
      int end = beg + seed.ref_pos.depth;
      seed_mem(seq, beg, end, seeds, rseed_opt, index, 0);
   }
   return 0;
}

int
seed_mem
(
 char         * seq,
 size_t         beg,
 size_t         end,
 seedstack_t ** stack,
 seedopt_t      opt,
 index_t      * index,
 int            reverse
)
{
   int32_t i = end - 1, qry_end = end - 1, last_qry_pos = end;
   uint64_t new_loci = index->size, last_loci;

   // Start position.
   bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
   while (i >= beg) {
      bwpos_t newpos;
      int nt = translate[(uint8_t)seq[i]];
      // Extend suffix (Backward search).
      suffix_extend(nt, pos, &newpos, index);
      // Count loci.
      last_loci = new_loci;
      new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);

      // No hits or depth exceeded.
      if (new_loci < opt.min_loci || nt == 4) {
         // Check previous suffix.
         int32_t qry_pos = qry_end + 1 - pos.depth;
         if (last_loci <= opt.aux_loci && pos.depth >= opt.min_len && qry_pos < last_qry_pos) {
            // Seed found.
            seed_t seed = (seed_t) {.bulk = (last_loci > opt.max_loci), .qry_pos = qry_pos, .ref_pos = pos};
            if (reverse) {
               seed.qry_pos = end - 1 - qry_end;
               seed.ref_pos.depth = -seed.ref_pos.depth;
            }
            // Save seed.
            seedstack_push(seed, stack);
            last_qry_pos = qry_pos;
         }
         // Shrink suffix.
         if (nt == 4) {
            i--;
            qry_end = i;
            pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
         } else {
            int depth = pos.depth;
            bwpos_t tmp;
            if (suffix_shrink(pos, &tmp, index) == 1) {
               suffix_string(seq + i + 1, qry_end - i, last_loci + 1, &tmp, index);
            }
            pos = tmp;
            qry_end -= depth - pos.depth;
         }
      } else if (newpos.depth == opt.max_len || i == beg) {
         // Update suffix.
         pos = newpos;
         // Check options.
         int32_t qry_pos = qry_end + 1 - pos.depth;
         if (new_loci <= opt.aux_loci && qry_pos < last_qry_pos && pos.depth >= opt.min_len) {
            // Seed found.
            seed_t seed = (seed_t) {.bulk = (new_loci > opt.max_loci), .qry_pos = qry_pos, .ref_pos = pos};
            if (reverse) {
               seed.qry_pos = end - 1 - qry_end;
               seed.ref_pos.depth = -seed.ref_pos.depth;
            }
            // Save seed.
            seedstack_push(seed, stack);
            last_qry_pos = qry_pos;
         }
         if (i == beg) break;

         // Shrink suffix.
         int depth = pos.depth;
         bwpos_t tmp;
         if (suffix_shrink(pos, &tmp, index) == 1) {
            suffix_string(seq + i, qry_end - i - 1, new_loci + 1, &tmp, index);
         }
         pos = tmp;
         qry_end -= depth - pos.depth;
         i--;
         /*
      } else if (new_loci <= opt.thr_seed && pos.depth >= opt.min_len) {
         // Update suffix.
         pos = newpos;
         // Check options.
         int32_t qry_pos = qry_end + 1 - pos.depth;
         if (new_loci < last_loci) {
            // Seed found.
            seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qry_pos, .ref_pos = pos};
            seedstack_push(seed, stack);
         }
         i--;
         */
      } else {
         // Update suffix.
         pos = newpos;
         i--;
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

   if (VERBOSE_DEBUG) fprintf(stdout, "seed list:\n");

   // Find loci in Suffix Array.
   hit_t * hits = malloc(cnt * sizeof(hit_t));
   size_t l = 0;
   for (size_t i = 0; i < seeds->pos; i++) {
      seed_t seed = seeds->seed[i];
      for (size_t j = seed.ref_pos.sp; j <= seed.ref_pos.ep; j++) {
         int64_t loc = get_sa(j,index->sa,index->sa_bits);
         if (seed.ref_pos.depth > 0) {
            hits[l++] = (hit_t) {.locus = loc, .qrypos = seed.qry_pos, .depth = seed.ref_pos.depth, .bulk = seed.bulk};
         }
         else {
            loc = index->size - 1 - loc - seed.ref_pos.depth;
            hits[l++] = (hit_t) {.locus = loc, .qrypos = seed.qry_pos, .depth = -seed.ref_pos.depth, .bulk = seed.bulk};
         }
      }
   }
   
   return hits;
}


int
match_seeds
(
 seedstack_t * seeds,
 int           seqlen,
 matchlist_t * matchlist,
 index_t     * index,
 int           dist_accept,
 double        read_ref_ratio
)
{
   if (matchlist->size < 1) return -1;

   // Convert seeds to hits.
   uint64_t hit_count;
   hit_t * hits = compute_hits(seeds, index, &hit_count);

   // Merge-sort hits.
   mergesort_mt(hits, hit_count, sizeof(hit_t), 0, 1, compar_hit_locus);

   // Reset matchlist.
   matchlist->pos = 0;

   // Aux variables.
   long minv = 0;
   int  min  = 0;
   size_t i = 0;

   while (i < hit_count) {
      // Match start.
      int64_t span   = hits[i].depth;
      int64_t lstart = hits[i].locus;
      int64_t lend   = lstart + span;
      int64_t rstart = hits[i].qrypos;
      int64_t rend   = rstart + span;
      int     bulk   = hits[i].bulk;
      i++;

      // Debug seeds...
      if (VERBOSE_DEBUG) {
         uint64_t g_start;
         int dir = 0;
         if (lstart >= index->size/2) {
            g_start = index->size - lend - 2;
            dir = 1;
         } else {
            g_start = lstart + 1;
            dir = 0;
         }
         // Search chromosome name.
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "%s,%ld%c/%ld,%ld/%ld,%ld,%d",
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 rstart, rend-1,
                 span, span, bulk
                 );
      }

      // Extend match.
      while (i < hit_count) {
         // Compare genome and read distances.
         int64_t g_dist = (int64_t)(hits[i].locus) - (int64_t)(hits[i-1].locus);
         int64_t r_dist = (int64_t)(hits[i].qrypos) - (int64_t)(hits[i-1].qrypos);
         if (VERBOSE_DEBUG) {
            fprintf(stdout, "\t[r:%ld,%ld(%ld) q:%d,%d(%ld)]", hits[i-1].locus, hits[i].locus, g_dist, hits[i-1].qrypos, hits[i].qrypos,r_dist);
         }
         if (r_dist >= 0 && (g_dist <= (read_ref_ratio * r_dist) || (g_dist < dist_accept && r_dist < dist_accept))) {
            //            int hitdiff = hits[i].depth - (lend - (int64_t)hits[i].locus > 0 ? lend - (int64_t)hits[i].locus : 0);
            //            span += (hitdiff > 0 ? hitdiff : 0);
            if (hits[i].locus < lend) span += hits[i].locus + hits[i].depth - lend;
            else span += hits[i].depth;
            lend = hits[i].locus + hits[i].depth;
            rend = hits[i].qrypos + hits[i].depth;
            bulk *= hits[i].bulk;
            // Debug seeds...
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

      if (VERBOSE_DEBUG) fprintf(stdout, "\n");

      // Store match.
      // Non-significant streaks (bulk) will not be saved.
      if (span > minv && bulk == 0) {
         // Allocate match.
         match_t match;
         match.ref_e  = lend;
         match.read_e = rend;
         match.ref_s  = lstart;
         match.read_s = rstart;
         match.hits   = span;
         match.flags  = 0;
         match.score  = -1;

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
         }
         else {
            matchlist->match[min] = match;
         }
               
         // Find the minimum that will be substituted next.
         if (matchlist->pos == matchlist->size) {
            minv = matchlist->match[0].hits;
            for (int j = 1 ; j < matchlist->pos; j++) {
               if (minv > matchlist->match[j].hits) {
                  min  = j;
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
         suffix_extend(nt, pos, &newpos, index);
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
