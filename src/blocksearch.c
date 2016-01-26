#include "blocksearch.h"

pathstack_t *
pathstack_new
(
 int n_elem
)
{
   // Set size to positive number.
   if (n_elem < 1) n_elem = 1;
   // Alloc structure.
   pathstack_t * stack = malloc(sizeof(pathstack_t) + n_elem * sizeof(spath_t));
   if (stack == NULL)
      return NULL;
   // Empty stack.
   stack->pos = 0;
   stack->size = n_elem;
   return stack;
}

int
path_push
(
 spath_t        p,
 pathstack_t ** pstack
)
{
   pathstack_t * stack = *pstack;
   // If full, realloc with double size.
   if (stack->pos >= stack->size) {
      size_t newsize = stack->size * 2;
      *pstack = stack = realloc(stack, sizeof(pathstack_t) + newsize * sizeof(spath_t));
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   // Store element.
   stack->path[stack->pos++] = p;
   return 0;
}

void
blocksearch
(
 uint8_t      * query,
 int            slen,
 int            tau,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Divide query in tau+1 blocks and perform recursive search.
   return blocksearch_rec(query, slen, tau+1, index, hits);
}

void
blocksearch_rec
(
 uint8_t      * query,
 int            slen,
 int            blocks,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Compute single block and return.
   if (blocks == 1) {
      spath_t empty = {.pos = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0}, .score = 0};
      seqsearch_bw(empty, query, 0, slen, 0, 0, 0, index, hits);
      return;
   }
   // Split block.
   int blk_l = (blocks >> 1) + (blocks & 1);
   int blk_r = (blocks >> 1);
   int beg_r = (slen >> 1) + (slen & 1);
   
   // Alloc temp hit buffer.
   pathstack_t * tmp_hits = pathstack_new(PATHSTACK_DEF_SIZE);
   // Recursive call on left block.
   blocksearch_rec(query, beg_r, blk_l, index, &tmp_hits);
   // Extend left block to the right.
   for (int i = 0; i < tmp_hits->pos; i++) {
      seqsearch_fw(tmp_hits->path[i], query+beg_r, 0, slen-beg_r, blocks-1, tmp_hits->path[i].score, 0, index, hits);
   }

   // Reset temp hits.
   tmp_hits->pos = 0;

   // Recursive call on right block.
   blocksearch_rec(query + beg_r, slen-beg_r, blk_r, index, &tmp_hits);
   // Extend right block to the left.
   for (int i = 0; i < tmp_hits->pos; i++) {
      seqsearch_bw(tmp_hits->path[i], query, 0, beg_r, blocks-1, tmp_hits->path[i].score, blk_l, index, hits);
   }

   // Free memory.
   free(tmp_hits);

   return;
}

int
seqsearch_bw
(
 spath_t        path,
 uint8_t      * query,
 int            depth,
 int            slen,
 int            tau,
 int            score_ref,
 int            score_diff,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Extend sequence.
   fmdpos_t tmp[NUM_BASES];
   extend_bw_all(path.pos, tmp, index);

   // Iterate over nt.
   for (int nt = 0; nt < NUM_BASES; nt++) {
      // Check whether prefix exists.
      if (tmp[nt].sz < 1)
         continue;
      // Update score.
      int s = path.score + (nt != query[slen-1-depth]);
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      int d = depth;
      spath_t p = {.pos = tmp[nt], .score = s};
      // Dash if max tau is reached.
      if (s == tau) {
         if (seqdash_bw(&p, query, d+1, slen, index))
            continue;
         d = slen-1;
      }
      // If depth is reached.
      if (d == slen-1) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         seqsearch_bw(p, query, d+1, slen, tau, score_ref, score_diff, index, hits);
      }
   }

   return 0;
}

int
seqdash_bw
(
 spath_t * path,
 uint8_t * query,
 int       depth,
 int       slen,
 index_t * index
)
{
   fmdpos_t p = path->pos;
   // Dash until the end of the search region.
   for (int i = slen-1-depth; i >= 0; i--) {
      // Do not allow any mismatch. If sequence does not exits, return.
      p = extend_bw(query[i], p, index);
      if (p.sz < 1) return -1;
   }
   // Sequence found, return path in pointer.
   path->pos = p;
   return 0;
}


int
seqsearch_fw
(
 spath_t        path,
 uint8_t      * query,
 int            depth,
 int            slen,
 int            tau,
 int            score_ref,
 int            score_diff,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Extend sequence.
   fmdpos_t tmp[NUM_BASES];
   extend_fw_all(path.pos, tmp, index);

   // Iterate over nt.
   for (int nt = 0; nt < NUM_BASES; nt++) {
      // Check whether prefix exists.
      if (tmp[nt].sz < 1)
         continue;
      // Update score.
      int s = path.score + (nt != query[depth]);
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      int d = depth;
      spath_t p = {.pos = tmp[nt], .score = s};
      // Dash if max tau is reached.
      if (s == tau) {
         if (seqdash_fw(&p, query, d+1, slen, index))
            continue;
         d = slen-1;
      }
      // If depth is reached.
      if (d == slen-1) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         seqsearch_fw(p, query, d+1, slen, tau, score_ref, score_diff, index, hits);
      }
   }

   return 0;
}

int
seqdash_fw
(
 spath_t * path,
 uint8_t * query,
 int       depth,
 int       slen,
 index_t * index
)
{
   fmdpos_t p = path->pos;
   // Dash until the end of the search region.
   for (int i = depth; i < slen; i++) {
      // Do not allow any mismatch. If sequence does not exits, return.
      p = extend_fw(query[i], p, index);
      if (p.sz < 1) return -1;
   }
   // Sequence found, return path in pointer.
   path->pos = p;
   return 0;
}
