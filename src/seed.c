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
   // Check arguments.
   error_test_msg(beg < 0, "argument 'beg' must be positive.");
   error_test_msg(end < 0, "argument 'end' must be positive.");
   error_test_msg(end < beg, "argument 'end' must greater than 'beg'.");
   error_test_msg(bwtq == NULL, "argument 'bwtq' is NULL.");
   error_test_msg(read == NULL, "argument 'read' is NULL.");

   // Declare variables.
   seed_t * seed = NULL;

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
   error_test_msg(last_mem == NULL, "argument 'last_mem' is NULL.");
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   error_test_msg(index == NULL, "argument 'index' is NULL.");
   
   // Determine MEM center position.
   int center = 0;
   if (last_mem != NULL) {
      center = seed_end(last_mem);
   }

   // Return NULL at end of sequence.
   if (center >= seqread_len(read)) {
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
   while (beg >= 0) {
      bwt_query(syms[beg], BWT_QUERY_PREFIX, bwtq, bwtq_next);
      // If query succeeded, save to 'bwtq'.
      if (bwt_size(bwtq_next) > 0) {
	 bwtquery_t * tmp = bwtq;
	 bwtq = bwtq_next;
	 bwtq_next = tmp;
	 beg--;
      }
      // Otherwise break search.
      else break;
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

   // Free memory.
   free(syms);
   free(bwtq_next);

   // Return valid seed or NULL if nothing was found.
   if (beg < end) {
      seed = seed_new(beg, end, bwtq, read);
      error_test(seed == NULL);
   }
   
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

   return stack;

 failure_return:
   gstack_free(stack);
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
