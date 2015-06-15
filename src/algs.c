#include "algs.h"

/*********************/
/** seq_t functions **/
/*********************/

long
bisect_search
(
 long   start,
 long   end,
 long * set,
 long   value
)
{
   if (start >= end - 1) {
      if (value < set[start]) return start;
      return (value >= set[end] ? end : start) + 1;
   }
   long middle = (start + end) / 2;
   if (set[middle] >= value) return bisect_search(start, middle, set, value);
   else return bisect_search(middle, end, set, value);
}

/*********************/
/** seq_t functions **/
/*********************/

int
seq_push
(
 seqstack_t ** stackp,
 const char  * tag,
 const char  * seq,
 const char  * q,
 const int     reverse
)
{
   char rcode[256] = {[0 ... 255] = 0,
                   ['a'] = 'T', ['c'] = 'G', ['g'] = 'C', ['t'] = 'A', ['u'] = 'A', ['n'] = 'N',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N' };

   seqstack_t * stack = *stackp;

   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(seqstack_t) + newsize*sizeof(seq_t));
      if (stack == NULL) {
         return -1;
      }
      stack->size = newsize;
   }

   seq_t * seqt = &(stack->seq[stack->pos++]);
   
   // Copy tag
   seqt->tag = strdup(tag);
   seqt->seq = strdup(seq);
   if (q) seqt->q = strdup(q);
   else seqt->q = calloc(strlen(seqt->seq),sizeof(char));

   // Copy sequence (or reverse-complement)
   if (reverse) {
      int len = strlen(seq);
      char * rseq = malloc(len+1);
      for (int i = 0; i < len; i++)
         rseq[len-1-i] = rcode[(int)seq[i]];
      rseq[len] = 0;
      seqt->rseq = rseq;
   } else {
      seqt->rseq = NULL;
   }
   return 0;
}

seqstack_t *
new_seqstack
(
 int size
 )
{
   if (size < 1) size = 1;
   seqstack_t * stack = malloc(sizeof(seqstack_t) + size * sizeof(seq_t));
   if (stack == NULL)
      return NULL;
   
   stack->pos = 0;
   stack->size = size;
   
   return stack;
}



/*********************/
/** stack functions **/
/*********************/


vstack_t *
new_stack
(
 long size
)
{
   if (size < 1) size = 1;
   vstack_t * stack = malloc(sizeof(vstack_t) + size * sizeof(long));
   if (stack == NULL) return NULL;
   stack->pos  = 0;
   stack->size = size;
   
   return stack;
}


int
push
(
 vstack_t ** stackp,
 long      value
)
{
   vstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(vstack_t) + newsize * sizeof(long));
      if (stack == NULL) {
         fprintf(stderr, "error in 'push' (realloc): %s\n", strerror(errno));
         return -1;
      }
      stack->size = newsize;
   }

   stack->val[stack->pos++] = value;
   return 0;
}

int
pushvec
(
 vstack_t ** stackp,
 long      * vector,
 int         vecsize
)
{
   vstack_t * stack = *stackp;
   if (stack->pos + vecsize >= stack->size) {
      long newsize = stack->size * 2;
      while(newsize <= stack->pos + vecsize) newsize *= 2;
      *stackp = stack = realloc(stack, sizeof(vstack_t) + newsize * sizeof(long));
      if (stack == NULL) {
         fprintf(stderr, "error in 'push' (realloc): %s\n", strerror(errno));
         return -1;
      }
      stack->size = newsize;
   }
   memcpy(stack->val + stack->pos, vector, vecsize * sizeof(long));
   stack->pos += vecsize;
   return 0;
}



pstack_t *
new_pstack
(
 long size
 )
{
   if (size < 1) size = 1;
   pstack_t * stack = malloc(sizeof(pstack_t) + size * sizeof(pebble_t));
   if (stack == NULL) return NULL;
   stack->pos  = 0;
   stack->size = size;

   return stack;
}

int
ppush
(
 pstack_t ** stackp,
 pebble_t pebble
 )
{
   pstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(pstack_t) + newsize * sizeof(pebble_t));
      if (stack == NULL) {
         fprintf(stderr, "error in 'ppush' (realloc): %s\n", strerror(errno));
         return -1;
      }
      stack->size = newsize;
   }

   stack->pebble[stack->pos++] = pebble;
   return 0;
}

/*********************/
/** trie  functions **/
/*********************/

trie_t *
trie_new
(
 int initial_size
)
// TODO: UPDATE.
// SYNOPSIS:                                                              
//   Creates and initializes a new NW-row trie preallocated with the specified
//   number of nodes.
//                                                                        
// PARAMETERS:                                                            
//   initial_size : the number of preallocated nodes.
//   height       : height of the trie. It must be equal to the number of keys
//                  as returned by parse.
//                                                                        
// RETURN:                                                                
//   On success, the function returns a pointer to the new trie_t structure.
//   A NULL pointer is returned in case of error.
//
// SIDE EFFECTS:
//   The returned trie_t struct is allocated using malloc and must be manually freed.
{

   // Allocate at least one node.
   if (initial_size < 1) initial_size = 1;

   trie_t * trie = malloc(sizeof(trie_t) + initial_size*sizeof(node_t));
   if (trie == NULL) {
      fprintf(stderr, "error in 'trie_new' (malloc) trie_t: %s\n", strerror(errno));
      return NULL;
   }

   // Initialize root node.
   memset(&(trie->nodes[0]), 0, initial_size*sizeof(node_t));

   // Initialize trie struct.
   trie->pos = 1;
   trie->size = initial_size;

   return trie;
}


int
trie_getrow
(
 trie_t * trie,
 uint     nodeid,
 int      refval,
 int    * wingsz,
 uint   * nwrow
)
// TODO: UPDATE.
// SYNOPSIS:                                                              
//   Recomputes the NW row that terminates at nodeid.
//                                                                        
// PARAMETERS:                                                            
//   trie   : Pointer to the trie.
//   nodeid : Id of the leaf at which the NW row terminates.
//   refval : The initial condition of the alignment. The score of
//            the rightmost value (top of the inverted L).
//                                                                        
// RETURN:                                                                
//   trie_getrow returns a pointer to the start of the NW row. If an
//   error occurred during the row computation or nodeid did not point
//   to a leaf node, a NULL pointer is returned.
//
// SIDE EFFECTS:
//   An array containing the NW row is allocated using malloc and must be
//   manually freed.
{
   // Return if path starts at root node.
   if (nodeid == 0) {
      *wingsz = 0;
      *nwrow  = refval;
      return 0;
   }

   node_t * nodes  = &(trie->nodes[0]);
   uint     parent = nodeid;
   int      height = 1;   
   while ((parent = nodes[parent].parent) != 0) height++;
   *wingsz = height/2;

   // Match value.
   int i = *wingsz;
   nwrow[i] = refval;
   uint id = nodeid;
   while (id != 0 && i > -(*wingsz)) {
      uint next_id = nodes[id].parent;
      nwrow[i-1] = nwrow[i] + (nodes[next_id].child[0] == id) - (nodes[next_id].child[2] == id);
      id = next_id;
      i--;
   } 

   // Control.
   if (i != -(*wingsz)) {
      fprintf(stderr, "error in 'trie_getrow': final node != root.\n");
      return -1;
   }

   return 0;
}


uint
trie_insert
(
 trie_t ** triep,
 char   *  path,
 int       pathlen
)
// TODO: UPDATE.
// SYNOPSIS:                                                              
//   Inserts the specified path in the trie and stores the end value and
//   the dfa state at the leaf (last node of the path). If the path already
//   exists, its leaf values will be overwritten.
//                                                                        
// PARAMETERS:                                                            
//   trie     : pointer to a memory space containing the address of the trie.
//   path     : The path as an array of chars containing values {0,1,2}
//   pathlen  : length of the path.
//                                                                        
// RETURN:                                                                
//   On success, dfa_insert returns the id of the leaf where the values were
//   stored, -1 is returned if an error occurred.
//
// SIDE EFFECTS:
//   If the trie has reached its limit of allocated nodes, it will be reallocated
//   doubling its size. The address of the trie may have changed after calling dfa_insert.
{
   trie_t * trie  = *triep;
   node_t * nodes = &(trie->nodes[0]);
   uint id = 0;
   uint initial_pos = trie->pos;

   int i;
   for (i = 0; i < pathlen; i++) {
      if (path[i] < 0 || path[i] >= TRIE_CHILDREN) {
         // Bad path, revert trie and return.
         trie->pos = initial_pos;
         return -1;
      }
      // Walk the tree.
      if (nodes[id].child[(int)path[i]] != 0) {
         id = nodes[id].child[(int)path[i]];
         continue;
      }

      // Create new node.
      if (trie->pos >= trie->size) {
         size_t newsize = trie->size * 2;
         *triep = trie = realloc(trie, sizeof(trie_t) + newsize * sizeof(node_t));
         if (trie == NULL) {
            fprintf(stderr, "error in 'trie_insert' (realloc) trie_t: %s\n", strerror(errno));
            return -1;
         }
         // Update pointers.
         nodes = &(trie->nodes[0]);
         // Initialize new nodes.
         trie->size = newsize;
      }
      
      // Consume one node of the trie.
      uint newid = trie->pos;
      nodes[newid].parent = id;
      nodes[newid].child[0] = nodes[newid].child[1] = nodes[newid].child[2] = 0;
      nodes[id].child[(int)path[i]] = newid;
      trie->pos++;

      // Go one level deeper.
      id = newid;
   }

   return id;
}


void
trie_reset
(
 trie_t * trie
)
// SYNOPSIS:                                                              
//   Resets the trie by pruning the root node. The size of the trie, in terms
//   of preallocated nodes is maintained.
//                                                                        
// PARAMETERS:                                                            
//   trie   : Pointer to the trie.
//                                                                        
// RETURN:                                                                
//   void.
//
// SIDE EFFECTS:
//   None.
{
   trie->pos = 1;
   memset(&(trie->nodes[0]), 0, sizeof(node_t));
}

/**************************/
/*** sorting algorithms ***/
/**************************/

int
mergesort_mt
(
 void  * data,
 int     numelm,
 size_t  elmsize,
 int     param,
 int     thrmax,
 int     (*compar)(const void *, const void *, const int)
)
// SYNOPSIS:                                                              
//   Recursive multithreaded merge sort for generic arrays.
//
// PARAMETERS:                                                            
//   data:   an array of seq_t.                    
//   numels: number of elements in data.                 
//   thrmax: number of threads.                                       
//                                                                        
// RETURN:                                                                
//   
//                                                                        
// SIDE EFFECTS:                                                          
//   Pointers to repeated elements are set to NULL.
{
   // Copy to buffer.
   void *buffer = malloc(numelm * elmsize);
   memcpy(buffer, data, numelm * elmsize);

   // Prepare args struct.
   sortargs_t args;
   args.buf0   = data;
   args.buf1   = buffer;
   args.size   = numelm;
   args.offset = elmsize;
   // There are two alternating buffers for the merge step.
   // 'args.b' alternates on every call to 'nukesort()' to
   // keep track of which is the source and which is the
   // destination. It has to be initialized to 0 so that
   // sorted elements end in 'data' and not in 'buffer'.
   args.b      = 0;
   args.thread = 0;
   args.param  = param;
   args.compar = compar;

   // Allocate a number of threads that is a power of 2.
   while ((thrmax >> (args.thread + 1)) > 0) args.thread++;

   _mergesort(&args);

   free(buffer);

   return 0;
}

void *
_mergesort
(
 void * args
)
// SYNOPSIS:
//   Recursive part of 'seqsort'.
//
// ARGUMENTS:
//   args: a sortargs_t struct.
//
// RETURN:
//   
//
// SIDE EFFECTS:
//   Sorts the array of 'seq_t' specified in 'args'.
{
   sortargs_t * sortargs = (sortargs_t *) args;
   if (sortargs->size < 2) return 0;

   size_t offset = sortargs->offset;

   // Next level params.
   sortargs_t arg1 = *sortargs, arg2 = *sortargs;
   arg1.size /= 2;
   arg2.size = arg1.size + arg2.size % 2;

   // Increment pointer positions.
   arg2.buf0 = (void *) (((char *) arg2.buf0) + arg1.size * offset);
   arg2.buf1 = (void *) (((char *) arg2.buf1) + arg1.size * offset);
   arg1.b = arg2.b = (arg1.b + 1) % 2;

   // Either run threads or DIY.
   if (arg1.thread) {
      // Decrease one level.
      arg1.thread = arg2.thread = arg1.thread - 1;
      // Create threads.
      pthread_t thread1, thread2;
      pthread_create(&thread1, NULL, _mergesort, &arg1);
      pthread_create(&thread2, NULL, _mergesort, &arg2);
      // Wait for threads.
      pthread_join(thread1, NULL); // TODO: Catch retval and return if retval == -1.
      pthread_join(thread2, NULL);
   }
   else {
      _mergesort(&arg1);
      _mergesort(&arg2);
   }

   // Separate data and buffer (b specifies which is buffer).
   char * l = (sortargs->b ? arg1.buf0 : arg1.buf1);
   char * r = (sortargs->b ? arg2.buf0 : arg2.buf1);
   char * buf = (sortargs->b ? arg1.buf1 : arg1.buf0);

   int i = 0;
   int j = 0;
   int idx = 0;
   int cmp = 0;

   // Merge sets
   while (idx < sortargs->size) {
      // Right buffer is exhausted. Copy left buffer...
      if (j == arg2.size) {
         memcpy(buf + idx*offset, l + i*offset, (arg1.size-i)*offset);
         break;
      }
      // ... or vice versa.
      if (i == arg1.size) {
         memcpy(buf + idx*offset, r + j*offset, (arg2.size-j)*offset);
         break;
      }
      // Do the comparison.
      cmp = arg1.compar(l + i*offset, r + j*offset, arg1.param);
      if (cmp < 0) memcpy(buf + (idx++)*offset, l + (i++)*offset, offset);
      else         memcpy(buf + (idx++)*offset, r + (j++)*offset, offset);
   }

   return NULL;
}

void
radix_sort
(
 long * a,      // Values to sort.
 long * b,      // Aux buffer.
 long   n,      // Length of a.
 long   maxval  // Maximum value in a.
)
{
   const int rs_bits = 16;
   const int rs_size = 1 << rs_bits;
   const int rs_mask = rs_size - 1;
 
   int ref = 0, new = 1, it = 0;
   long cnt[rs_size];
   long prf[rs_size];
   long * s[2];

   // Count iterations per value.
   while ((maxval >> (rs_bits*it)) & rs_mask) it++;

   s[0] = a;
   s[1] = b;
   for (long j = 0; j < it; j++) {
      int shift = rs_bits * j;
      // Reset count and prefix.
      memset(cnt, 0, rs_size*sizeof(long));
      prf[0] = 0;
      // Count radix RS_BITS.
      for (long i = 0; i < n; i++) cnt[(s[ref][i] >> shift) & rs_mask]++;
      // Prefix.
      for (int  i = 1; i < rs_size; i++) prf[i] = prf[i-1] + cnt[i-1];
      // Sorted.
      for (long i = 0; i < n; i++) s[new][prf[(s[ref][i] >> shift) & rs_mask]++] = s[ref][i];
      // Swap buffers.
      ref = (ref+1)%2;
      new = (new+1)%2;
   }
   
   // Move data to a.
   if (ref == 1) memcpy(s[0], s[1], n*sizeof(long));
}
