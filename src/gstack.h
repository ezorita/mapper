#include <stdlib.h>
#include <stdint.h>
#include "errhandler.h"

#ifndef _GSTACK_H
#define _GSTACK_H

#define GSTACK_DEFAULT_SIZE 1024

// Typedef structures.
typedef struct gstack_t gstack_t;

// Type interfaces.
struct gstack_t;

// Interface functions.
gstack_t    * gstack_new           (size_t max_elm, void (*elm_free)(void *));
void          gstack_free          (gstack_t * gstack);
void        * gstack_pop           (gstack_t * gstack);
void        * gstack_get           (size_t idx, gstack_t * gstack);
int           gstack_push          (void * elm_ptr, gstack_t * gstack);
int           gstack_push_array    (void ** base_ptr, size_t elm_cnt, gstack_t * gstack);
int           gstack_popr          (gstack_t * gstack);
int           gstack_transfer_all  (gstack_t * dst, gstack_t * src);
int           gstack_clear         (gstack_t * gstack);

// Helper functions.
int64_t       gstack_num_elm       (gstack_t * gstack);
int64_t       gstack_max_elm       (gstack_t * gstack);

#endif
