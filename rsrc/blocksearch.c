#include "blocksearch.h"

/*
** Block-search functions
*/

void
blocksc_trail
(
 uint8_t     * query,
 bwtquery_t ** qarray,
 int           slen,
 int           tau,
 int           trail,
 pstree_t    * tree
)
// Seeq & Construct version of blocksearch algorithm over a Forward/Reverse index.
// This version shall only be used for exhaustive genome neighbors computaion.

// Out of tau mismatches, make first division and assign mismatches as follows:
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

   // Index info.
   bwt_t * bwt = bwt_get_bwt(qarray[slen]);

   // Reset hits. (Free bwtquery first)
   for (int i = 0; i < tree->stack->pos; i++) {
      free(tree->stack->path[i].bwtq);
   }
   tree->stack->pos = 0;

   // Which strand.
   int last_fragment = bwt_start(qarray[slen]) >= bwt_rcstart(qarray[slen]);

   // If the query contains N in either block (left or right), reduce the value of
   // tau by num(N) mismatches.
   int n_cnt = 0;
   for (int i = 0; i < slen && n_cnt <= tau; i++) n_cnt += query[i] == UNKNOWN_BASE;

   // Return if num(N) > tau or tau == 0 and not in last strand.
   tau -= n_cnt;
   if (tau < 0 || (tau == 0 && !last_fragment)) return;

   // Split block.
   int pos_r = slen/2 + (last_fragment ? slen%2 : 0);
   int tau_l =  tau/2 - (last_fragment ? 0 : (1 - tau % 2));

   // Recursive call on left block, don't compute if current data is valid (trail).
   if (trail < pos_r) {
      blocksearch_trail_rec(query, 0, pos_r-1, tau_l+1, trail, bwt, tree->next_l);
      // Remove out-of-boundary hits (lexicographically bigger than query).
      int64_t max_sa_pos = bwt_start(qarray[pos_r]) + bwt_size(qarray[pos_r]);
      int i = 0;
      while (i < tree->next_l->stack->pos) {
         if (bwt_start(tree->next_l->stack->path[i].bwtq) < max_sa_pos) {
            i++;
         } else {
            free(tree->next_l->stack->path[i].bwtq);
            tree->next_l->stack->path[i] = tree->next_l->stack->path[--tree->next_l->stack->pos];
         }
      }
   }

   // Extend left block to the right.
   for (int i = 0; i < tree->next_l->stack->pos; i++) {
      spath_t p = tree->next_l->stack->path[i];
      scsearch_fw(p, query, pos_r, slen-1, tau, p.score, 0, 1, &(tree->stack));
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
 bwt_t     * bwt,
 pstree_t  * tree
)
{
   // Reset hits. (Free bwtquery first)
   for (int i = 0; i < tree->stack->pos; i++) {
      free(tree->stack->path[i].bwtq);
   }
   tree->stack->pos = 0;

   // Compute single block and return.
   if (blocks == 1) {
      spath_t empty = {.bwtq = bwt_new_query(bwt), .score = 0};
      seqsearch_bw(empty, query, end, pos, 0, 0, 0, &(tree->stack));
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
      blocksearch_trail_rec(query, pos, end_l, blk_l, trail, bwt, tree->next_l);
   // Extend left block to the right.
   for (int i = 0; i < tree->next_l->stack->pos; i++) {
      spath_t p = tree->next_l->stack->path[i];
      seqsearch_fw(p, query, pos_r, end, blocks-1, p.score, 0, &(tree->stack));
   }

   // Recursive call on right block.
   blocksearch_trail_rec(query, pos_r, end, blk_r, trail, bwt, tree->next_r);
   // Extend right block to the left.
   for (int i = 0; i < tree->next_r->stack->pos; i++) {
      spath_t p = tree->next_r->stack->path[i];
      seqsearch_bw(p, query, end_l, pos, blocks-1, p.score, blk_l, &(tree->stack));
   }

   return;
}


/*
** Index query functions
*/

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
 pathstack_t ** hits
)
{
   bwtquery_t  * q   = path.bwtq;
   bwt_t       * bwt = bwt_get_bwt(q);
   bwtquery_t ** qv  = bwt_new_vec(bwt);

   // Extend sequence.
   bwt_query_all(BWT_QUERY_PREFIX, q, qv);

   // Iterate over nt.
   int num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   for (int nt = 0; nt < num_symb; nt++) {
      // Check whether prefix exists.
      if (bwt_size(qv[nt]) < 1)
         continue;
      // Update score.
      int ds = (nt != query[pos] && query[pos] != UNKNOWN_BASE);
      int  s = path.score + ds;
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      spath_t p = {.bwtq = qv[nt], .score = s};
      memcpy(&(p.align), &(path.align), ALIGN_WORDS*sizeof(uint64_t));
      // Update mismatched positions.
      if (ds || query[pos] == UNKNOWN_BASE) aln_bit_set(&p, pos);
      // Dash if max tau is reached.
      if (s == tau) {
         if (s - score_ref >= score_diff) {
            seqdash_bw(p, query, pos-1, end, hits);
         }
      }
      // If depth is reached.
      else if (pos == end) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff) {
            path_push(p, hits);
         }
      } else {
         // Recursive call.
         seqsearch_bw(p, query, pos-1, end, tau, score_ref, score_diff, hits);
      }
   }

   bwt_free_vec(qv);

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
 pathstack_t ** hits
)
{
   bwtquery_t  * q   = path.bwtq;
   bwt_t       * bwt = bwt_get_bwt(q);
   bwtquery_t ** qv  = bwt_new_vec(bwt);

   // Extend sequence.
   bwt_query_all(BWT_QUERY_SUFFIX, q, qv);

   // Iterate over nt.
   int num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   for (int nt = 0; nt < num_symb; nt++) {
      // Check whether prefix exists.
      if (bwt_size(qv[nt]) < 1)
         continue;
      // Update score.
      int ds = (nt != query[pos] && query[pos] != UNKNOWN_BASE);
      int  s = path.score + ds;
      // Check score boundary.
      if (s > tau)
         continue;
      // Update path.
      spath_t p = {.bwtq = qv[nt], .score = s};
      memcpy(&(p.align), &(path.align), ALIGN_WORDS*sizeof(uint64_t));
      // Update mismatched positions.
      if (ds || query[pos] == UNKNOWN_BASE) aln_bit_set(&p, pos);
      // Dash if max tau is reached.
      if (s == tau) {
         if (s - score_ref >= score_diff)
            seqdash_fw(p, query, pos+1, end, hits);
      }
      // If depth is reached.
      else if (pos == end) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         seqsearch_fw(p, query, pos+1, end, tau, score_ref, score_diff, hits);
      }
   }

   bwt_free_vec(qv);

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
 pathstack_t ** hits
)
{
   bwtquery_t  * q   = path.bwtq;
   bwt_t       * bwt = bwt_get_bwt(q);
   bwtquery_t ** qv  = bwt_new_vec(bwt);

   // Extend sequence.
   bwt_query_all(BWT_QUERY_SUFFIX, q, qv);

   // Iterate over nt.
   int num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   for (int nt = 0; nt < (boundary ? query[pos] + 1 : num_symb); nt++) {
      // Check whether prefix exists.
      if (bwt_size(qv[nt]) < 1)
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
      spath_t p = {.bwtq = qv[nt], .score = s};
      memcpy(&(p.align), &(path.align), ALIGN_WORDS*sizeof(uint64_t));
      if (ds || query[pos] == UNKNOWN_BASE) aln_bit_set(&p, pos);
      // Dash if max tau is reached.
      if (s == tau) {
         if (s - score_ref >= score_diff)
            seqdash_fw(p, query, pos+1, end, hits);
      }
      // If depth is reached.
      else if (pos == end) {
         // Check score difference and store hit.
         if (s - score_ref >= score_diff)
            path_push(p, hits);
      } else {
         // Recursive call.
         int bnd = boundary && (nt == query[pos]);
         scsearch_fw(p, query, pos+1, end, tau, score_ref, score_diff, bnd, hits);
      }
   }

   bwt_free_vec(qv);

   return 0;
}


int
seqdash_fw
(
 spath_t        path,
 uint8_t      * query,
 int            pos,
 int            end,
 pathstack_t ** hits
)
{
   bwtquery_t  * q   = path.bwtq;
   bwt_t       * bwt = bwt_get_bwt(q);

   // Dash until the end of the search region.
   for (int d = pos; d <= end; d++) {
      // Do not allow any mismatch. If sequence does not exist, return.
      if (__builtin_expect(query[d] != UNKNOWN_BASE,1)) {
         bwt_query(query[d], BWT_QUERY_SUFFIX, q, q);
         if (bwt_size(q) < 1)
            return 0;
      } 
      // When querying UNKNOWN_BASE, follow all branches.
      else {
         // Set mismatch position.
         aln_bit_set(&path, d);

         // Extend all branches.
         bwtquery_t ** qv  = bwt_new_vec(bwt);
         bwt_query_all(BWT_QUERY_SUFFIX, q, qv);

         // Follow branches.
         int num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));
         for (int i = 0; i < num_symb; i++) {
            if (bwt_size(qv[i]) < 1)
               continue;
            path.bwtq = qv[i];
            seqdash_fw(path, query, d+1, end, hits);
         }
         
         // Free query vector.
         bwt_free_vec(qv);
         return 0;
      }
   }
   // Sequence found, push to hit stack.
   path.bwtq = q;
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
 pathstack_t ** hits
)
{
   bwtquery_t  * q   = path.bwtq;
   bwt_t       * bwt = bwt_get_bwt(q);

   // Dash until the end of the search region.
   for (int d = pos; d >= end; d--) {
      // Do not allow any mismatch. If sequence does not exist, return.
      if (__builtin_expect(query[d] != UNKNOWN_BASE,1)) {
         bwt_query(query[d], BWT_QUERY_PREFIX, q, q);
         if (bwt_size(q) < 1)
            return 0;
      }
      // When querying UNKNOWN_BASE, follow all branches.
      else {
         // Set mismatch position.
         aln_bit_set(&path, d);

         // Extend all branches.
         bwtquery_t ** qv  = bwt_new_vec(bwt);
         bwt_query_all(BWT_QUERY_PREFIX, q, qv);

         // Follow branches.
         int num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));
         for (int i = 0; i < num_symb; i++) {
            if (bwt_size(qv[i]) < 1)
               continue;
            path.bwtq = qv[i];
            seqdash_bw(path, query, d-1, end, hits);
         }
         
         // Free query vector.
         bwt_free_vec(qv);
         return 0;
      }
   }
   // Sequence found, push to hit stack.
   path.bwtq = q;
   path_push(path, hits);
   return 0;
}


/*
** Stack tree functions
*/

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
   for (int i = 0; i < tree->stack->pos; i++) {
      free(tree->stack->path[i].bwtq);
   }
   free(tree->stack);
   if (tree->next_l != NULL) free_stack_tree(tree->next_l);
   if (tree->next_r != NULL) free_stack_tree(tree->next_r);
   free(tree);
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
   stack->path[stack->pos] = p;
   stack->path[stack->pos].bwtq = bwt_dup_query(p.bwtq);
   stack->pos++;
   return 0;
}


/*
** Misc functions.
*/

void
aln_bit_set
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

