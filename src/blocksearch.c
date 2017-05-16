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
blocksc_trail
(

 uint8_t   * query,
 fmdpos_t  * qarray,
 int         slen,
 int         tau,
 int         trail,
 index_t   * index,
 pstree_t  * tree
)
// Seeq & Construct version of blocksearch algorithm over a Forward/Reverse index.
// This version shall only be used for exhaustive genome neighbors computaion.

// Out of the tau mismatches, make first division and assign mismatches as follows:
// 
//  - Sequence is lexicographically smaller compared to its revcomp:
//     Left block:
//       length = L/2 
//       tau    = floor(tau_max/2) - (1 - tau_max % 2)
//     Then extend right up to tau_max.
//
//  - Sequence is lexicographically bigger or equal compared to its revcomp:
//     Left block:
//       length = L/2 + L%2
//       tau    = floor(tau_max/2)
//     Then extend right up to tau_max.
{
   if (trail >= slen) return;

   // Reset hits.
   tree->stack->pos = 0;

   // Which strand.
   int rc_last = qarray[slen].fp >= qarray[slen].rp;

   // If the query contains N in either block (left or right), reduce the value of
   // tau by num(N) mismatches.
   int n_cnt = 0;
   for (int i = 0; i < slen && n_cnt <= tau; i++) n_cnt += query[i] == UNKNOWN_BASE;

   // Return if num(N) > tau or tau == 0 and not in last strand.
   tau -= n_cnt;
   if (tau < 0 || (tau == 0 && !rc_last)) return;

   // Split block.
   int pos_r = slen/2 + (rc_last ? slen%2 : 0);
   int tau_l = tau/2 - (rc_last ? 0 : (1 - tau % 2));

   // Recursive call on left block, don't compute if current data is valid (trail).
   if (trail < pos_r) {
      blocksearch_trail_rec(query, 0, pos_r-1, tau_l+1, trail, index, tree->next_l);
      // Remove out-of-boundary hits (lexicographically bigger than query).
      int64_t max_sa_pos = qarray[pos_r].fp + qarray[pos_r].sz;
      int i = 0;
      while (i < tree->next_l->stack->pos) {
         if (tree->next_l->stack->path[i].pos.fp < max_sa_pos)
            i++;
         else
            tree->next_l->stack->path[i] = tree->next_l->stack->path[--tree->next_l->stack->pos];
      }
   }

   // Extend left block to the right.
   for (int i = 0; i < tree->next_l->stack->pos; i++) {
      spath_t p = tree->next_l->stack->path[i];
      scsearch_fw(p, query, pos_r, slen-1, tau, p.score, 0, 1, index, &(tree->stack));
   }

   // Add num(N) mismatches to hits.
   if (n_cnt)
      for (int i = 0; i < tree->stack->pos; i++)
         tree->stack->path[i].score += n_cnt;

   return;
}


void
blocksearch_trail_rec
(
 uint8_t   * query,
 int         pos,
 int         end,
 int         blocks,
 int         trail,
 index_t   * index,
 pstree_t  * tree
)
{
   // Reset hits.
   tree->stack->pos = 0;

   // Compute single block and return.
   if (blocks == 1) {
      spath_t empty = {.pos = index->bwt->fmd_base, .score = 0};
      seqsearch_bw(empty, query, end, pos, 0, 0, 0, index, &(tree->stack));
      return;
   }

   // Split block.
   int blk_l = (blocks >> 1) + (blocks & 1);
   int blk_r = (blocks >> 1);

   // Read params.
   int slen = end - pos + 1;
   int pos_r = pos + (slen >> 1) + (slen & 1);
   int end_l = pos_r - 1;

   // Recursive call on left block, don't compute if current data is valid (trail).
   if (trail < pos_r)
      blocksearch_trail_rec(query, pos, end_l, blk_l, trail, index, tree->next_l);
   // Extend left block to the right.
   for (int i = 0; i < tree->next_l->stack->pos; i++) {
      spath_t p = tree->next_l->stack->path[i];
      seqsearch_fw(p, query, pos_r, end, blocks-1, p.score, 0, index, &(tree->stack));
   }

   // Recursive call on right block.
   blocksearch_trail_rec(query, pos_r, end, blk_r, trail, index, tree->next_r);
   // Extend right block to the left.
   for (int i = 0; i < tree->next_r->stack->pos; i++) {
      spath_t p = tree->next_r->stack->path[i];
      seqsearch_bw(p, query, end_l, pos, blocks-1, p.score, blk_l, index, &(tree->stack));
   }

   return;
}

pstree_t *
alloc_stack_tree
(
 int tau
)
{
   return alloc_stack_tree_rec(tau+1);
}

pstree_t *
alloc_stack_tree_rec
(
 int        block
)
{
   // End of tree, return null.
   if (block == 0) return NULL;
   // Alloc stack for this node.
   pstree_t * node = calloc(1,sizeof(pstree_t));
   node->stack = pathstack_new(PATHSTACK_DEF_SIZE);
   // Build tree with recursive calls.
   if (block > 1) {
      node->next_l = alloc_stack_tree_rec((block >> 1) + (block & 1));
      node->next_r = alloc_stack_tree_rec(block >> 1);
   } else {
      node->next_l = node->next_r = NULL;
   }
   return node;
}

void
free_stack_tree
(
 pstree_t * tree
)
{
   free(tree->stack);
   if (tree->next_l != NULL) free_stack_tree(tree->next_l);
   if (tree->next_r != NULL) free_stack_tree(tree->next_r);
   free(tree);
}

void
mpos_set
(
 spath_t * path,
 int       pos
)
{
   int w = pos/ALIGN_WORD_SIZE;
   int b = pos%ALIGN_WORD_SIZE;
   path->align[w] |= (uint64_t)1 << b;
   return;
}

int
seqsearch_bw
(
 spath_t        path,
 uint8_t      * query,
 int            pos,
 int            end,
 int            tau,
 int            score_ref,
 int            score_diff,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Extend sequence.
   fmdpos_t tmp[NUM_BASES];
   extend_bw_all(path.pos, tmp, index->bwt);

   // Iterate over nt.
   for (int nt = 0; nt < NUM_BASES; nt++) {
      // Check whether prefix exists.
      if (tmp[nt].sz < 1)
         continue;
      // Update score.
      int ds = (nt != query[pos] && query[pos] != UNKNOWN_BASE);
      int  s = path.score + ds;
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      spath_t p = {.pos = tmp[nt], .score = s};
      memcpy(&(p.align), &(path.align), ALIGN_WORDS*sizeof(uint64_t));
      // Update mismatched positions.
      if (ds || query[pos] == UNKNOWN_BASE) mpos_set(&p, pos);      
      // Dash if max tau is reached.
      if (s == tau) {
         if (s - score_ref >= score_diff)
            seqdash_bw(p, query, pos-1, end, index, hits);
         continue;
      }
      // If depth is reached.
      if (pos == end) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         seqsearch_bw(p, query, pos-1, end, tau, score_ref, score_diff, index, hits);
      }
   }

   return 0;
}

int
seqsearch_fw
(
 spath_t        path,
 uint8_t      * query,
 int            pos,
 int            end,
 int            tau,
 int            score_ref,
 int            score_diff,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Extend sequence.
   fmdpos_t tmp[NUM_BASES];
   extend_fw_all(path.pos, tmp, index->bwt);

   // Iterate over nt.
   for (int nt = 0; nt < NUM_BASES; nt++) {
      // Check whether prefix exists.
      if (tmp[nt].sz < 1)
         continue;
      // Update score.
      int ds = (nt != query[pos] && query[pos] != UNKNOWN_BASE);
      int  s = path.score + ds;
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      spath_t p = {.pos = tmp[nt], .score = s};
      memcpy(&(p.align), &(path.align), ALIGN_WORDS*sizeof(uint64_t));
      // Update mismatched positions.
      if (ds || query[pos] == UNKNOWN_BASE) mpos_set(&p, pos);
      // Dash if max tau is reached.
      if (s == tau) {
         if (s - score_ref >= score_diff)
            seqdash_fw(p, query, pos+1, end, index, hits);
         continue;
      }
      // If depth is reached.
      if (pos == end) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         seqsearch_fw(p, query, pos+1, end, tau, score_ref, score_diff, index, hits);
      }
   }

   return 0;
}

int
scsearch_fw
(
 spath_t        path,
 uint8_t      * query,
 int            pos,
 int            end,
 int            tau,
 int            score_ref,
 int            score_diff,
 int            boundary,
 index_t      * index,
 pathstack_t ** hits
)
{
   // Extend sequence.
   fmdpos_t tmp[NUM_BASES];
   extend_fw_all(path.pos, tmp, index->bwt);

   // Iterate over nt.
   for (int nt = 0; nt < (boundary ? query[pos] + 1 : NUM_BASES); nt++) {
      // Check whether prefix exists.
      if (tmp[nt].sz < 1)
         continue;
      // Update score.
      // 'N' is treated as wildcard -- The maximum distance should have been lowered
      // by the number of 'N' present in the query.
      int ds = (nt != query[pos] && query[pos] != UNKNOWN_BASE);
      int  s = path.score + ds;
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      spath_t p = {.pos = tmp[nt], .score = s};
      memcpy(&(p.align), &(path.align), ALIGN_WORDS*sizeof(uint64_t));
      if (ds || query[pos] == UNKNOWN_BASE) mpos_set(&p, pos);
      // Dash if max tau is reached.
      if (s == tau) {
         if (s - score_ref >= score_diff)
            seqdash_fw(p, query, pos+1, end, index, hits);
         continue;
      }
      // If depth is reached.
      if (pos == end) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         int bnd = boundary && (nt == query[pos]);
         scsearch_fw(p, query, pos+1, end, tau, score_ref, score_diff, bnd, index, hits);
      }
   }

   return 0;
}

int
seqdash_fw
(
 spath_t        path,
 uint8_t      * query,
 int            pos,
 int            end,
 index_t      * index,
 pathstack_t ** hits
)
{
   fmdpos_t p = path.pos;
   // Dash until the end of the search region.
   for (int d = pos; d <= end; d++) {
      // Do not allow any mismatch. If sequence does not exist, return.
      if (__builtin_expect(query[d] != UNKNOWN_BASE,1)) {
         p = extend_fw(query[d], p, index->bwt);
         if (p.sz < 1) return 0;
      } 
      // When querying UNKNOWN_BASE, follow all branches.
      else {
         // Set mismatch position.
         mpos_set(&path, d);
         // Extend all branches.
         fmdpos_t newpos[NUM_BASES];
         extend_fw_all(p, newpos, index->bwt);
         // Follow branches.
         for (int i = 0; i < NUM_BASES; i++) {
            if (newpos[i].sz < 1) continue;
            path.pos = newpos[i];
            seqdash_fw(path, query, d+1, end, index, hits);
         }
         return 0;
      }
   }
   // Sequence found, push to hit stack.
   path.pos = p;
   path_push(path, hits);
   return 0;
}

int
seqdash_bw
(
 spath_t        path,
 uint8_t      * query,
 int            pos,
 int            end,
 index_t      * index,
 pathstack_t ** hits
)
{
   fmdpos_t p = path.pos;
   // Dash until the end of the search region.
   for (int d = pos; d >= end; d--) {
      // Do not allow any mismatch. If sequence does not exist, return.
      if (__builtin_expect(query[d] != UNKNOWN_BASE,1)) {
         p = extend_bw(query[d], p, index->bwt);
         if (p.sz < 1) return -1;
      }
      // When querying UNKNOWN_BASE, follow all branches.
      else {
         // Set mismatch position.
         mpos_set(&path, d);
         // Extend all branches.
         fmdpos_t newpos[NUM_BASES];
         extend_bw_all(p, newpos, index->bwt);
         // Follow branches.
         for (int i = 0; i < NUM_BASES; i++) {
            if (newpos[i].sz < 1) continue;
            path.pos = newpos[i];
            // Dash all branches.
            seqdash_bw(path, query, d-1, end, index, hits);
         }
         return 0;
      }
   }
   // Sequence found, push to hit stack.
   path.pos = p;
   path_push(path, hits);
   return 0;
}
