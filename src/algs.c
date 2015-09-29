#include "algs.h"

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
 const char  * q
)
{
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

// Matchlist functions.

matchlist_t *
matchlist_new
(
 int elements
)
{
   if (elements < 1) elements = 1;
   matchlist_t * list = malloc(sizeof(matchlist_t) + elements*sizeof(match_t));
   if (list == NULL) return NULL;
   list->pos  = 0;
   list->size = elements;
   return list;
}

int
matchlist_add
(
 matchlist_t ** listp,
 match_t        match
)
{
   matchlist_t * list = *listp;

   // Check whether stack is full.
   if (list->pos >= list->size) {
      int newsize = list->size * 2;
      *listp = list = realloc(list, sizeof(matchlist_t) + newsize * sizeof(match_t));
      if (list == NULL) return -1;
      list->size = newsize;
   }

   // Add new match to the stack.
   list->match[list->pos++] = match;

   return 0;
}
