#include "blocksearch.h"

int
which_strand
(
 uint8_t * query,
 int       kmer
 )
/*
** This function takes a seq and compares it with its reverse
** complement. Returns 1 if the reverse complement is 
** lexicographically smaller/equal than the original sequence.
** The aim of this is to avoid computing the same thing twice.
** Since we have an index with forward and reverse strands,
** the forward and reverse complement of each kmer must yield
** exactly the same number of hits. Therefore we will only
** compute half of the job for the lexicographically smaller
** and the other half (maybe +1 nucleotide if k is odd) for
** the legicographically bigger.
*/
{
   for (int i = 0, j = kmer-1; i < kmer; i++, j--) {
      if (3-query[j] == query[i]) continue;
      return 3-query[j] < query[i];
   }
   return 1;
}

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
blocksearch_trail_sc
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

// If the query contains N in either block (left or right), reduce the value of
// tau by num(N) mismatches.
{
   if (trail >= slen) return;

   // Reset hits.
   tree->stack->pos = 0;

   // Which strand.
   int rc_last = which_strand(query, slen);

   // Split block.
   int len = slen/2 + (rc_last ? slen%2 : 0);
   int tau = tau_max/2 - (rc_last ? 0 : (1 - tau_max % 2));

   // Recursive call on left block, don't compute if current data is valid (trail).
   if (trail < len) {
      blocksearch_trail_rec(query, pos, beg_r, blk_l, trail, index, tree->next_l);
      // Remove out-of-boundary hits (lexicographically bigger than query).
      int64_t max_sa_pos = qarray[len].fp + qarray[len].sz;
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
      seqsearch_fw_sc(p, query+len, 0, slen-len, tau, p.score, 0, index, &(tree->stack));
   }

   return;
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

void
blocksearch_trail
(
 uint8_t   * query,
 int         slen,
 int         tau,
 int         trail,
 index_t   * index,
 pstree_t  * tree
)
{
   if (trail >= slen) return;
   if (trail < 0) trail = 0;
   return blocksearch_trail_rec(query,0,slen,tau+1,trail,index,tree);
}

void
blocksearch_trail_rec
(
 uint8_t   * query,
 int         pos,
 int         slen,
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
      spath_t empty = {.pos = (fmdpos_t){.fp = 0, .rp = 0, .sz = index->size, .dp = 0}, .score = 0};
      seqsearch_bw(empty, query+pos, 0, slen, 0, 0, 0, index, &(tree->stack));
      return;
   }
   // Split block.
   int blk_l = (blocks >> 1) + (blocks & 1);
   int blk_r = (blocks >> 1);
   int beg_r = (slen >> 1) + (slen & 1);

   // Recursive call on left block, don't compute if current data is valid (trail).
   if (trail < pos+beg_r)
      blocksearch_trail_rec(query, pos, beg_r, blk_l, trail, index, tree->next_l);
   // Extend left block to the right.
   for (int i = 0; i < tree->next_l->stack->pos; i++) {
      spath_t p = tree->next_l->stack->path[i];
      seqsearch_fw(p, query+pos+beg_r, 0, slen-beg_r, blocks-1, p.score, 0, index, &(tree->stack));
   }

   // Recursive call on right block.
   blocksearch_trail_rec(query, pos + beg_r, slen-beg_r, blk_r, trail, index, tree->next_r);
   // Extend right block to the left.
   for (int i = 0; i < tree->next_r->stack->pos; i++) {
      spath_t p = tree->next_r->stack->path[i];
      seqsearch_bw(p, query+pos, 0, beg_r, blocks-1, p.score, blk_l, index, &(tree->stack));
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
   extend_bw_all(path.pos, tmp, index->bwt);

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
      p = extend_bw(query[i], p, index->bwt);
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
   extend_fw_all(path.pos, tmp, index->bwt);

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
seqsearch_fw_sc
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
   extend_fw_all(path.pos, tmp, index->bwt);

   // Iterate over nt.
   for (int nt = 0; nt < (score_ref ? NUM_BASES : query[depth] + 1); nt++) {
      // Check whether prefix exists.
      if (tmp[nt].sz < 1)
         continue;
      // Update score.
      // 'N' is treated as wildcard -- The maximum distance should have been lowered
      // by the number of 'N' present in the query.
      int s = path.score + (nt != query[depth] && query[depth] < 4);
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
         seqsearch_fw_sc(p, query, d+1, slen, tau, score_ref, score_diff, index, hits);
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
      p = extend_fw(query[i], p, index->bwt);
      if (p.sz < 1) return -1;
   }
   // Sequence found, return path in pointer.
   path->pos = p;
   return 0;
}

int
compar_path_score
(
 const void * pa,
 const void * pb
)
{
   spath_t * a = (spath_t *)pa;
   spath_t * b = (spath_t *)pb;

   if (a->score < b->score) return -1;
   else if (a->score > b->score) return 1;
   else return (a->pos.sz < b->pos.sz ? -1 : 1);
}
