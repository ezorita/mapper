#include "seed.h"

seedstack_t *
seed
(
 char      * seq,
 seedopt_t   opt,
 index_t   * index
)
{
   uint32_t slen = strlen(seq);
   int32_t i = slen-1, qry_end = slen-1, last_qry_pos = slen;
   uint64_t new_loci;

   // Create stack.
   seedstack_t * stack = seedstack_new(SEEDSTACK_SIZE);

   // Start position.
   bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
   while (i >= 0) {
      bwpos_t newpos;
      int nt = translate[(uint8_t)seq[i]];
      // Extend suffix (Backward search).
      suffix_extend(nt, pos, &newpos, index);
      // Count loci.
      new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);

      // No hits or depth exceeded.
      if (new_loci < opt.min_loci || nt == 4) {
         // Check previous suffix.
         uint64_t loci = (pos.ep < pos.sp ? 0 : pos.ep - pos.sp + 1);
         int32_t qry_pos = qry_end + 1 - pos.depth;
         if (loci <= opt.aux_loci && pos.depth >= opt.min_len && qry_pos < last_qry_pos) {
            // Seed found.
            seed_t seed = (seed_t) {.bulk = (loci > opt.max_loci), .qry_pos = qry_pos, .ref_pos = pos};
            seedstack_push(seed, &stack);
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
               suffix_string(seq + i + 1, qry_end - i, loci + 1, &tmp, index);
            }
            pos = tmp;
            qry_end -= depth - pos.depth;
         }
      } else if (newpos.depth == opt.max_len || i == 0) {
         // Update suffix.
         pos = newpos;
         // Check options.
         int32_t qry_pos = qry_end + 1 - pos.depth;
         if (new_loci <= opt.aux_loci && qry_pos < last_qry_pos) {
            // Seed found.
            seed_t seed = (seed_t) {.bulk = (new_loci > opt.max_loci), .qry_pos = qry_pos, .ref_pos = pos};
            seedstack_push(seed, &stack);
            last_qry_pos = qry_pos;
            if (i == 0) break;
         }
         if (i == 0) break;

         // Shrink suffix.
         int depth = pos.depth;
         bwpos_t tmp;
         if (suffix_shrink(pos, &tmp, index) == 1) {
            suffix_string(seq + i, qry_end - i - 1, new_loci + 1, &tmp, index);
         }
         pos = tmp;
         qry_end -= depth - pos.depth;
         i--;
      } else {
         // Update suffix.
         pos = newpos;
         i--;
      }
   } while (i > 0);

   return stack;
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
         hits[l++] = (hit_t) {.locus = get_sa(j,index->sa,index->sa_bits), .qrypos = seed.qry_pos, .depth = seed.ref_pos.depth, .bulk = seed.bulk};
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
      // Extend match.
      while (i < hit_count) {
         // Compare genome and read distances.
         int64_t g_dist = (int64_t)hits[i].locus - (int64_t)hits[i-1].locus;
         int64_t r_dist = (int64_t)hits[i].qrypos - (int64_t)hits[i-1].qrypos;
         if (r_dist > 0 && (g_dist < (read_ref_ratio * r_dist) || (g_dist < dist_accept && r_dist < dist_accept))) {
            span += hits[i].depth - (lend - (int64_t)hits[i].locus > 0 ? lend - (int64_t)hits[i].locus : 0);
            lend = hits[i].locus + hits[i].depth;
            rend = hits[i].qrypos + hits[i].depth;
            bulk *= hits[i].bulk;
            i++;
         } else break;
      }
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
   else return (ha->qrypos > hb->qrypos ? 1 : -1);
}
