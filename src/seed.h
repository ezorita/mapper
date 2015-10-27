#include <stdint.h>
#include <string.h>
#include "algs.h"
#include "indexquery.h"

#ifndef _SEED_H
#define _SEED_H

#define SEEDSTACK_SIZE 64

typedef struct seedopt_t   seedopt_t;
typedef struct seed_t      seed_t;
typedef struct seedstack_t seedstack_t;
typedef struct hit_t       hit_t;

struct seedopt_t {
   int32_t min_len;
   int32_t max_len;
   int32_t min_loci;
   int32_t max_loci;
   int32_t aux_loci;
   int32_t reseed_len;
   uint64_t thr_seed;
};

struct seed_t {
   int32_t bulk;
   int32_t qry_pos;
   bwpos_t ref_pos;
};

struct seedstack_t {
   size_t pos;
   size_t size;
   seed_t seed[];
};

struct hit_t {
   uint64_t locus;
   uint32_t qrypos;
   uint16_t depth;
   uint16_t bulk;
};


seedstack_t  * naive_smem     (char *, seedopt_t, index_t *);
// Seed functions.
int            seed_mem       (char *, int32_t, int32_t, seedstack_t **, seedopt_t, index_t *);
int            seed_wings     (char *, int32_t, int32_t, seedstack_t **, seedopt_t, index_t *);
int            reseed_mem     (char *, seedstack_t **, seedopt_t, index_t *);
int            reseed_smem_rec     (char *, size_t, size_t, seedstack_t **, seedopt_t, index_t *);
int            get_smem       (seedstack_t **, seedstack_t *);
int            seed_thr       (char *, size_t, seedstack_t **, seedopt_t, index_t *);
hit_t        * compute_hits   (seedstack_t * seeds, index_t * index, uint64_t * hit_cnt);
int            chain_seeds    (seedstack_t *, int, matchlist_t *, index_t *, int, double);
int            find_repeats   (seedstack_t * seeds, matchlist_t ** repeats, int32_t repeat_thr);
// Seedstack functions.
seedstack_t  * seedstack_new    (int);
int            seedstack_push   (seed_t, seedstack_t **);
int            seedstack_append (seedstack_t **, seedstack_t *);
// Mergesort compar functions.
int            compar_hit_locus  (const void * a, const void * b, const int param);
int            compar_seed_start (const void * a, const void * b, const int param);
// Aux functions.
char *         rev_comp         (char * seq);

#endif
