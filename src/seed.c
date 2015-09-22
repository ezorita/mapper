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
