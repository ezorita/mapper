#include "gstack.h"

// Define data structures.

struct gstack_t {
   size_t   elm_size;
   size_t   num_elm;
   size_t   max_elm;
   void   * ptr;
};


// Interface functions source.

gstack_t *
gstack_new
(
  size_t elm_size,
  size_t max_elm
)
{
   // Declare variables.
   gstack_t * gstack  = NULL;
   void     * buffer = NULL;

   // Check arguments.
   error_test_msg(elm_size == 0, "argument 'elm_size' must be greater than 0.");
   error_test_msg(max_elm == 0, "argument 'max_elm' must be greater than 0.");

   // Alloc gstack and set defaults.
   gstack = malloc(sizeof(gstack_t));
   error_test_mem(gstack);

   gstack->elm_size = elm_size;
   gstack->num_elm  = 0;
   gstack->max_elm  = max_elm;
   gstack->ptr      = NULL;

   // Alloc buffer and set in gstack.
   buffer = malloc(max_elm * elm_size);
   error_test_mem(buffer);

   gstack->ptr = buffer;

   return gstack;
   
 failure_return:
   gstack_free(gstack);
   return NULL;
}

void
gstack_free
(
  gstack_t * gstack
)
{
   if (gstack != NULL) {
      if (gstack->ptr != NULL) {
         free(gstack->ptr);
      }
      free(gstack);
   }

   return;
}

void *
gstack_pop
(
  gstack_t * gstack
)
/*
** Returns a copy of the element.
*/
{
   // Empty gstacks return NULL (this is not an error).
   if (gstack == NULL || gstack->num_elm == 0) {
      return NULL;
   }
   
   // Declare variables.
   void * elm = NULL;

   // Make a copy of the last element.   
   elm = malloc(gstack->elm_size);
   error_test_msg_errno(elm == NULL, "memory error.", GSTACK_ERRNO_MEMORY);

   void * elm_src = (void *)(((uint8_t *)gstack->ptr) + gstack->elm_size * --gstack->num_elm);
   memcpy(elm, elm_src, gstack->elm_size);

   return elm;

 failure_return:
   free(elm);
   return NULL;
}

void *
gstack_get
(
  size_t      idx,
  gstack_t  * gstack
)
/*
** Returns a pointer to the element in the gstack.
** This is only valid until the element is removed from the gstack!
*/
{
   // Check arguments.
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   error_test_msg(idx >= gstack->num_elm, "index out of bounds.");

   // Get item pointer and return.
   return (void *)(((uint8_t *)gstack->ptr) + gstack->elm_size * idx);

 failure_return:
   return NULL;
}

int
gstack_push
(
  void      * elm_ptr,
  gstack_t  * gstack
)
{
   // Check arguments.
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   error_test_msg(elm_ptr == NULL, "argument 'elm_ptr' is NULL.");

   // Check gstack current size.
   if (gstack->num_elm >= gstack->max_elm) {
      size_t new_size   = gstack->max_elm * 2;
      void * new_buffer = realloc(gstack->ptr, new_size*gstack->elm_size);
      error_test_mem(new_buffer);

      gstack->ptr     = new_buffer;
      gstack->max_elm = new_size;
   }

   // Add element to gstack.
   void * dst_ptr = (void *)(((uint8_t *)gstack->ptr) + gstack->elm_size * gstack->num_elm);
   memcpy(dst_ptr, elm_ptr, gstack->elm_size);
   gstack->num_elm++;

   return 0;

 failure_return:
   return -1;
}

int
gstack_push_array
(
  void      * base_ptr,
  size_t      elm_cnt,
  gstack_t  * gstack
)
{
   // Check arguments.
   if (elm_cnt == 0)
      return 0;
   error_test_msg(base_ptr == NULL, "argument 'base_ptr' is NULL.");
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");

   // Check gstack current size.
   if (gstack->num_elm + elm_cnt >= gstack->max_elm) {
      size_t new_size   = gstack->max_elm * 2;
      while (gstack->num_elm + elm_cnt >= new_size) {
         new_size *= 2;
      }
      void * new_buffer = realloc(gstack->ptr, new_size*gstack->elm_size);
      error_test_mem(new_buffer);

      gstack->ptr     = new_buffer;
      gstack->max_elm = new_size;
   }

   // Add element to gstack.
   void * dst_ptr = (void *)(((uint8_t *)gstack->ptr) + gstack->elm_size * gstack->num_elm);
   memcpy(dst_ptr, base_ptr, elm_cnt * gstack->elm_size);
   gstack->num_elm += elm_cnt;

   return 0;

 failure_return:
   return -1;
}

int
gstack_popr
(
  gstack_t  * gstack
)
/*
** Removes top element without returning it.
*/
{
   // Check arguments.
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   
   if (gstack->num_elm > 0)
      gstack->num_elm--;

   return 0;
   
 failure_return:
   return -1;
}

int
gstack_clear
(
  gstack_t  * gstack
)
/*
** Resets num_elm index to 0.
*/
{
   // Check arguments.
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   
   gstack->num_elm = 0;

   return 0;
   
 failure_return:
   return -1;
}

// Helper functions source.

int64_t
gstack_num_elm
(
  gstack_t  * gstack
)
{
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   
   return (int64_t)gstack->num_elm;

 failure_return:
   return -1;
}

int64_t
gstack_max_elm
(
  gstack_t  * gstack
)
{
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   
   return (int64_t)gstack->max_elm;

 failure_return:
   return -1;
}

int64_t
gstack_elm_size
(
  gstack_t  * gstack
)
{
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");
   
   return (int64_t)gstack->elm_size;

 failure_return:
   return -1;
}
