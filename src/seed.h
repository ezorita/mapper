#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "errhandler.h"
#include "gstack.h"
#include "seqread.h"
#include "index.h"


#ifndef _SEED_H
#define _SEED_H

#define SEEDSTACK_DEFAULT_SIZE 64

// Typedef structures.
typedef struct seed_t seed_t;

// Type interfaces.
struct seed_t;

// Interface functions.
seed_t     * seed_new      (int64_t beg, int64_t end, bwtquery_t * bwtq, seqread_t * read);
void         seed_free     (void * seed);
gstack_t   * seed_stack    (size_t max_elm);
seed_t     * seed_pop      (gstack_t * stack);
seed_t     * seed_get      (size_t index, gstack_t * stack);
seed_t     * seed_next_mem (seed_t * last_mem, seqread_t * read, index_t * index);
gstack_t   * seed_mems     (seqread_t * read, index_t * index);

// Helper functions.
int64_t      seed_beg      (seed_t * seed);
int64_t      seed_end      (seed_t * seed);
bwtquery_t * seed_bwtq     (seed_t * seed);
seqread_t  * seed_read     (seed_t * seed);

#endif
