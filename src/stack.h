#include <stdlib.h>
#include <stdint.h>
#include "errhandler.h"

#ifndef _STACK_H
#define _STACK_H

#define STACK_DEFAULT_SIZE 1024

// Typedef structures.
typedef struct stack_t stack_t;

// Type interfaces.
struct stack_t;

// Interface functions.
stack_t     * stack_new         (size_t elm_size, size_t max_elm);
void          stack_free        (stack_t * stack);
void        * stack_pop         (stack_t * stack);
void        * stack_get         (size_t idx, stack_t * stack);
int           stack_push        (void * elm_ptr, stack_t * stack);
int           stack_push_array  (void * base_ptr, size_t elm_cnt, stack_t * stack);

// Helper functions.
int64_t       stack_num_elm     (stack_t * stack);
int64_t       stack_max_elm     (stack_t * stack);
int64_t       stack_elm_size    (stack_t * stack);

#endif
