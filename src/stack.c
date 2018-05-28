#include "stack.h"

// Define data structures.

struct stack_t {
   size_t   elm_size;
   size_t   num_elm;
   size_t   max_elm;
   void   * ptr;
};


// Interface functions source.

stack_t *
stack_new
(
  size_t elm_size,
  size_t max_elm
)
{
   // Declare variables.
   stack_t * stack  = NULL;
   void    * buffer = NULL;

   // Alloc stack and set defaults.
   stack = malloc(sizeof(stack_t));
   error_test_mem(stack);

   stack->elm_size = elm_size;
   stack->num_elm  = num_elm;
   stack->max_elm  = max_elm;
   stack->ptr      = NULL;

   // Alloc buffer and set in stack.
   buffer = malloc(max_elm * elm_size);
   error_test_mem(buffer);

   stack->ptr = buffer;

   return stack;
   
 failure_return:
   stack_free(stack);
   return NULL;
}

void
stack_free
(
  stack_t * stack
)
{
   if (stack != NULL) {
      if (stack->ptr != NULL) {
         free(stack->ptr);
      }
      free(stack);
   }

   return;
}

void *
stack_pop
(
  stack_t * stack
)
{
   // Empty stacks return NULL, this is not an ERROR!
   if (stack == NULL || stack->num_elm == 0) {
      return NULL;
   }

   return (void *)(((uint8_t *)stack->ptr) + stack->elm_size * --stack->num_elm);
}

void *
stack_get
(
  size_t     idx,
  stack_t  * stack
)
{
   // Check arguments.
   error_test_msg(stack == NULL, "argument 'stack' is NULL.");
   error_test_msg(idx >= stack->num_elm, "index out of bounds.");

   // Get item pointer and return.
   return (void *)(((uint8_t *)stack->ptr) + stack->elm_size * idx);

 failure_return:
   return NULL;
}

int
stack_push
(
  void     * elm_ptr,
  stack_t  * stack
)
{
   // Check arguments.
   error_test_msg(stack == NULL, "argument 'stack' is NULL.");

   // Check stack current size.
   if (stack->num_elm >= stack->max_elm) {
      size_t new_size   = stack->max_elm * 2;
      void * new_buffer = realloc(stack->ptr, new_size*stack->elm_size);
      error_test_mem(new_buffer);

      stack->ptr     = new_buffer;
      stack->max_elm = new_size;
   }

   // Add element to stack.
   void * dst_ptr = (void *)(((uint8_t *)stack->ptr) + stack->elm_size * stack->num_elm);
   memcpy(dst_ptr, elm_ptr, stack->elm_size);
   stack->num_elm++;

   return 0;

 failure_return:
   return -1;
}

int
stack_push_array
(
  void     * base_ptr,
  size_t     elm_cnt,
  stack_t  * stack
)
{
   // Check arguments.
   if (elm_cnt == 0)
      return 0;
   error_test_msg(base_ptr == NULL, "argument 'base_ptr' is NULL.");
   error_test_msg(stack == NULL, "argument 'stack' is NULL.");

   // Check stack current size.
   if (stack->num_elm + elm_cnt >= stack->max_elm) {
      size_t new_size   = stack->max_elm * 2;
      while (stack->mum_elm + elm_cnt >= new_size) {
         new_size *= 2;
      }
      void * new_buffer = realloc(stack->ptr, new_size*stack->elm_size);
      error_test_mem(new_buffer);

      stack->ptr     = new_buffer;
      stack->max_elm = new_size;
   }

   // Add element to stack.
   void * dst_ptr = (void *)(((uint8_t *)stack->ptr) + stack->elm_size * stack->num_elm);
   memcpy(dst_ptr, base_ptr, elm_cnt * stack->elm_size);
   stack->num_elm += elm_cnt;

   return 0;

 failure_return:
   return -1;
}


// Helper functions source.

int64_t
stack_num_elm
(
  stack_t  * stack
)
{
   error_test_msg(stack == NULL, "argument 'stack' is NULL.");
   
   return (int64_t)stack->num_elm;

 failure_return:
   return -1;
}

int64_t
stack_max_elm
(
  stack_t  * stack
)
{
   error_test_msg(stack == NULL, "argument 'stack' is NULL.");
   
   return (int64_t)stack->max_elm;

 failure_return:
   return -1;
}

int64_t
stack_elm_size
(
  stack_t  * stack
)
{
   error_test_msg(stack == NULL, "argument 'stack' is NULL.");
   
   return (int64_t)stack->elm_size;

 failure_return:
   return -1;
}
