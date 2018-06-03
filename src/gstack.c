#include "gstack.h"

// Define data structures.

struct gstack_t {
   size_t    num_elm;
   size_t    max_elm;
   void      (*elm_free)(void *);
   void   ** elm;

};

// Interface functions source.

gstack_t *
gstack_new
(
  size_t   max_elm,
  void     (*elm_free)(void *)
)
{
   // Declare variables.
   gstack_t * gstack  = NULL;

   // Check arguments.
   error_test_msg(max_elm == 0, "argument 'max_elm' must be greater than 0.");

   // Alloc gstack and set defaults.
   gstack = malloc(sizeof(gstack_t));
   error_test_mem(gstack);

   gstack->num_elm  = 0;
   gstack->max_elm  = max_elm;
   gstack->elm_free = elm_free;
   gstack->elm      = NULL;

   // Alloc buffer and set in gstack.
   gstack->elm = malloc(max_elm * sizeof(void *));
   error_test_mem(gstack->elm);

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
      if (gstack->elm != NULL) {
	 if (gstack->elm_free != NULL) {
	    gstack_clear(gstack);
	 }
         free(gstack->elm);
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
   
   return gstack->elm[--gstack->num_elm];
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
   return gstack->elm[idx];

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
      void * new_buffer = realloc(gstack->elm, new_size*sizeof(void *));
      error_test_mem(new_buffer);

      gstack->elm     = new_buffer;
      gstack->max_elm = new_size;
   }

   // Add element to gstack.
   gstack->elm[gstack->num_elm++] = elm_ptr;
   return 0;

 failure_return:
   return -1;
}

int
gstack_push_array
(
  void     ** base_ptr,
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
      void * new_buffer = realloc(gstack->elm, new_size*sizeof(void *));
      error_test_mem(new_buffer);

      gstack->elm     = new_buffer;
      gstack->max_elm = new_size;
   }

   // Add element to gstack.
   memcpy(gstack->elm + gstack->num_elm, base_ptr, elm_cnt * sizeof(void *));
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
** Removes and frees top element without returning it.
*/
{
   // Check arguments.
   error_test_msg(gstack == NULL, "argument 'gstack' is NULL.");

   // Free and remove top element.
   if (gstack->num_elm > 0) {
      gstack->elm_free(gstack->elm[--gstack->num_elm]);
   }

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

   // Free remaining elements.
   for (int i = 0; i < gstack->num_elm; i++) {
      gstack->elm_free(gstack->elm[i]);
   }
   
   // Reset num elements.
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
