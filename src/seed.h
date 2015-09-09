#include <stdint.h>
#include <string.h>
#include "indexquery.h"

#ifndef _SEED_H
#define _SEED_H

#define SEEDSTACK_SIZE 64

typedef struct seedopt_t   seedopt_t;
typedef struct seed_t      seed_t;
typedef struct seedstack_t seedstack_t;


struct seedopt_t {
   int32_t min_len;
   int32_t max_len;
   int32_t min_loci;
   int32_t max_loci;
};

struct seed_t {
   int32_t qry_pos;
   bwpos_t ref_pos;
};

struct seedstack_t {
   size_t pos;
   size_t size;
   seed_t seed[];
};

seedstack_t  * seed           (char *, seedopt_t, index_t *);
seedstack_t  * seedstack_new  (int);
int            seedstack_push (seed_t, seedstack_t **);

#endif
