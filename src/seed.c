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
   uint32_t i = slen-1, qry_pos = 0;

   // Create stack.
   seedstack_t * stack = seedstack_new(SEEDSTACK_SIZE);

   // Start position.
   bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
   while (i >= 0) {
      bwpos_t newpos;
      // Extend suffix (Backward search).
      suffix_extend(translate[(uint8_t)seq[i]], pos, &newpos, index);
      // Count loci.
      uint64_t new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);
      // No hits or depth exceeded.
      if (new_loci < opt.min_loci) {
         // Check previous suffix.
         uint64_t loci = (pos.ep < pos.sp ? 0 : pos.ep - pos.sp + 1);
         if (loci <= opt.max_loci && pos.depth >= opt.min_len) {
            // Seed found.
            seed_t seed = (seed_t) {.qry_pos = qry_pos, .ref_pos = pos};
            seedstack_push(seed, &stack);
         }
         // Shrink suffix.
         int depth = pos.depth;
         suffix_shrink(pos, &pos, index);
         qry_pos += depth - pos.depth;
      } else if (newpos.depth == opt.max_len) {
         // Update suffix.
         pos = newpos;
         i--;
         // Check options.
         if (new_loci <= opt.max_loci) {
            // Seed found.
            seed_t seed = (seed_t) {.qry_pos = qry_pos, .ref_pos = pos};
            seedstack_push(seed, &stack);
         }
         // Shrink suffix.
         int depth = pos.depth;
         suffix_shrink(pos, &pos, index);
         qry_pos += depth - pos.depth;
      } else {
         // Update suffix.
         pos = newpos;
         i--;
      }
   } while (i > 0);

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
