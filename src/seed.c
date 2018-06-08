#include "seed.h"

// Define data structures.

struct seed_t {
   int64_t      beg;
   int64_t      end;
   bwtquery_t * bwtq;
   seqread_t  * read;
};

// Interface functions.

seed_t *
seed_new
(
 int64_t      beg,
 int64_t      end,
 bwtquery_t * bwtq,
 seqread_t  * read
)
{
   // Declare variables.
   seed_t * seed = NULL;

   // Check arguments.
   error_test_msg(beg < 0, "argument 'beg' must be positive.");
   error_test_msg(end < 0, "argument 'end' must be positive.");
   error_test_msg(end < beg, "argument 'end' must greater than 'beg'.");
   error_test_msg(bwtq == NULL, "argument 'bwtq' is NULL.");
   error_test_msg(read == NULL, "argument 'read' is NULL.");

   // Alloc memory.
   seed = malloc(sizeof(seed_t));
   error_test_mem(seed);

   seed->beg = beg;
   seed->end = end;
   seed->bwtq = bwtq;
   seed->read = read;

   // Return seed.
   return seed;

 failure_return:
   free(seed);
   return NULL;
}


void
seed_free
(
  void * seedptr
)
{
   seed_t * seed = (seed_t *) seedptr;
   if (seed != NULL) {
      free(seed->bwtq);
      free(seed);
   }

   return;
}


gstack_t *
seed_stack
(
  size_t max_elm
)
{
   gstack_t * stack = gstack_new(max_elm, seed_free);
   error_test(stack == NULL);
   
   return stack;

 failure_return:
   return NULL;
}


seed_t *
seed_pop
(
  gstack_t * stack
)
{
   seed_t * seed = (seed_t *) gstack_pop(stack);
   return seed;
}


seed_t *
seed_get
(
  size_t     index,
  gstack_t * stack
)
{
   seed_t * seed = (seed_t *) gstack_get(index, stack);
   error_test(seed == NULL);

   return seed;
   
 failure_return:
   return NULL;
}


seed_t *
seed_next_mem
(
 seed_t    * last_mem,
 seqread_t * read,
 index_t   * index
)
{
   // Declare variables.
   seed_t     * seed = NULL;
   bwtquery_t * bwtq = NULL;
   bwtquery_t * bwtq_next = NULL;
   uint8_t    * syms = NULL;
   
   // Check arguments.
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   error_test_msg(index == NULL, "argument 'index' is NULL.");
   
   // Determine MEM center position.
   int center = 0;
   if (last_mem != NULL) {
      center = seed_end(last_mem);
   }

   // Return NULL at end of sequence.
   if (center >= seqread_len(read)) {
      errno = SEED_ERRNO_END;
      return NULL;
   }
   
   // Get symbols from read.
   syms = sym_str_index(seqread_seq(read), index->sym);
   error_test(syms == NULL);

   // Begin index query.
   bwtq = bwt_new_query(index->bwt);
   error_test(bwtq == NULL);

   bwtq_next = bwt_new_query(index->bwt);
   error_test(bwtq_next == NULL);
   
   // Extend max to the left, starting from center.
   int beg = center;
   while (beg >= 0 && beg < seqread_len(read)) {
      bwt_query(syms[beg], BWT_QUERY_PREFIX, bwtq, bwtq_next);
      // If query succeeded, save to 'bwtq'.
      if (bwt_size(bwtq_next) > 0) {
	 bwtquery_t * tmp = bwtq;
	 bwtq = bwtq_next;
	 bwtq_next = tmp;
	 beg--;
      }
      // Otherwise break search.
      else {
	 // Move center, this only happens when we query a symbol
	 // that is not present in the index.
	 if (beg == center && beg < seqread_len(read)) {
	    center++;
	    beg++;
	 } else break;
      }
   }
   // Make beg point to the leftmost queried symbol.
   beg++;
   
   // Extend max to the right.
   int end = center + 1;
   while (end < seqread_len(read)) {
      bwt_query(syms[end], BWT_QUERY_SUFFIX, bwtq, bwtq_next);
      if (bwt_size(bwtq_next) > 0) {
	 bwtquery_t * tmp = bwtq;
	 bwtq = bwtq_next;
	 bwtq_next = tmp;
	 end++;
      }
      // Otherwise break search.
      else break;
   }


   // Return valid seed or NULL if nothing was found.
   if (bwt_depth(bwtq) > 0) {
      seed = seed_new(beg, end, bwtq, read);
      error_test(seed == NULL);
   } else {
      errno = SEED_ERRNO_NOT_FOUND;
      free(bwtq);
   }

   // Free memory.
   free(syms);
   free(bwtq_next);

   return seed;
   
 failure_return:
   seed_free(seed);
   free(bwtq);
   free(bwtq_next);   
   free(syms);
   return NULL;
}


gstack_t *
seed_mems
(
 seqread_t * read,
 index_t   * index
)
{
   // Declare variables.
   gstack_t * stack = NULL;

   // Check arguments.
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   error_test_msg(index == NULL, "argument 'index' is NULL.");

   // Alloc stack.
   stack = seed_stack(SEEDSTACK_DEFAULT_SIZE);
   error_test(stack == NULL);

   // Iterate over MEMs.
   seed_t * seed = NULL;
   while ((seed = seed_next_mem(seed, read, index)) != NULL) {
      error_test(gstack_push(seed, stack) == -1);
   }

   // Check cause of loop break.
   if (errno != SEED_ERRNO_END && errno != SEED_ERRNO_NOT_FOUND) {
      error_throw();
   }

   return stack;

 failure_return:
   gstack_free(stack);
   return NULL;
}


gstack_t **
seed_ann_dist
(
  gstack_t  * seeds,
  index_t   * index
)
/*
** This function queries the annotation index and classifies the seeds
** inside an array of stacks by their distance to the closest neighbor.
**
** Issue: This function copies the content of the stack to another stack,
** so using gstack_free on both will produce a double free error.
**
** Solution: Transfer the seeds, read the first stack and copy all of its
** content to the new stacks, then set num_elm to 0 and max_elm to 1, and
** realloc the old stack.
*/
{
   // Declare variables.
   gstack_t ** seed_dist = NULL;
   seed_t    * s         = NULL;
   int         max_d     = 0;
   
   // Check arguments.
   error_test_msg(seeds == NULL, "argument 'seeds' is NULL.");
   error_test_msg(index == NULL, "argument 'index' is NULL.");

   // Return the same stack when no info is available.
   if (gstack_num_elm(seeds) == 0 || index->ann_cnt == 0) {
      // Create stack array.
      seed_dist = calloc(sizeof(gstack_t *), 2);
      error_test_mem(seed_dist);
      
      // Initialize stack array.
      seed_dist[0] = seed_stack(1);
      error_test(seed_dist[0]);

      // Transfer contents.
      gstack_transfer_all(seed_dist[0], seeds);
      
      return seed_dist;
   }

   // Find max distance.
   for (int i = 0; i < index->ann_cnt; i++) {
      int d = ann_get_dist(index->ann[i]) + 1;
      if (d > max_d) max_d = d;
   }

   // Alloc stack arrays.
   seed_dist = calloc(sizeof(gstack_t *), max_d + 2);
   error_test_mem(seed_dist);
      
   // Initialize stack array.
   for (int i = 0; i <= max_d; i++) {
      seed_dist[i] = seed_stack(1);
      error_test(seed_dist[i] == NULL);
   }

   // Find fragments with far-neighbors in seeds.
   while((s = seed_pop(seeds)) != NULL) {
      // Get bwtq
      bwtquery_t * bwtq = seed_bwtq(s);
     
      // Exact neighbors.
      if (bwt_size(bwtq) > 1) {
	 error_test(gstack_push(s, seed_dist[0]) == -1);
	 continue;
      }

      // Otherwise, check annotations.
      int64_t s_dst = 0;
      int64_t s_len = seed_len(s);
      // Find locus.
      int64_t s_loc = sar_get(bwt_start(bwtq), index->sar);
      // Iterate over all annotations.
      for (int a = 0; a < index->ann_cnt; a++) {
	 ann_t * ann      = index->ann[a];
	 int     ann_kmer = ann_get_kmer(ann);
	 // Only annotations with k <= seed_length.
	 if (ann_kmer <= s_len) {
	    // Check all annotations on overlaping positions.
	    for (int i = 0; i <= s_len-ann_kmer; i++) {
	       locinfo_t * loc = ann_query(s_loc+i, ann);
	       error_test_mem(loc);
	       s_dst = (loc->dist > s_dst ? loc->dist : s_dst);
	       free(loc);
	    }
	 }
      }

      // Store seed in stack based on max neighbor distance.
      error_test(gstack_push(s, seed_dist[s_dst]) == -1);
      // Delete reference to stored seed (avoids double free in case of error).
      s = NULL;
   }

   return seed_dist;

 failure_return:
   for (int i = 0; i <= max_d; i++) {
      gstack_free(seed_dist[i]);
   }
   free(seed_dist);
   seed_free(s);
   return NULL;
}


// Helper function sources.

int64_t
seed_beg
(
  seed_t * seed
)
{
   // Check arguments.
   error_test_msg(seed == NULL, "argument 'seed' is NULL.");

   return seed->beg;

 failure_return:
   return -1;
}

int64_t
seed_end
(
  seed_t * seed
)
{
   // Check arguments.
   error_test_msg(seed == NULL, "argument 'seed' is NULL.");

   return seed->end;

 failure_return:
   return -1;
}

int64_t
seed_len
(
  seed_t * seed
)
{
   // Check arguments.
   error_test_msg(seed == NULL, "argument 'seed' is NULL.");

   return seed->end - seed->beg;

 failure_return:
   return -1;
}


bwtquery_t *
seed_bwtq
(
  seed_t * seed
)
{
   // Check arguments.
   error_test_msg(seed == NULL, "argument 'seed' is NULL.");

   return seed->bwtq;

 failure_return:
   return NULL;
}

seqread_t *
seed_read
(
  seed_t * seed
)
{
   // Check arguments.
   error_test_msg(seed == NULL, "argument 'seed' is NULL.");

   return seed->read;

 failure_return:
   return NULL;
}
