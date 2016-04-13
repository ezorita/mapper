#include <stdint.h>
#include <string.h>
#include "algs.h"
#include "definitions.h"
#include "bwtquery.h"
#include "blocksearch.h"
#include "index.h"

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
   int32_t thr_loci;
   int16_t reseed_len;
   int16_t sensitive_mem;
};

struct seed_t {
   int32_t errors;
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
   uint16_t errors;
};

// Annotated seed.
int            find_uniq_seeds (int, uint8_t *, uint8_t *, index_t *, seedstack_t **);
hit_t          find_uniq_seed  (int, int, uint8_t *, uint8_t *, index_t *);
int            find_thr_seed   (int, int, int, uint8_t *, index_t *);

int            mem_unique      (seed_t, index_t *, hit_t *);
// Longest MEM strategy.
seed_t         longest_mem     (char *, index_t *);
// Naive algorithms. (To compare with)
seedstack_t  * naive_smem     (char *, seedopt_t, index_t *);
// Seed functions.
int            seed_mem       (uint8_t *, int, int, index_t *, seedstack_t **);
int            next_mem       (uint8_t *, int, int, seed_t  *, index_t *);
int            extend_lr      (uint8_t *, int, int, fmdpos_t *, index_t *);
int            reseed_mem     (uint8_t *, seed_t, int, int, index_t *, seedstack_t **);
int            seed_mem_bp    (uint8_t *, int, int, int, index_t *, seedstack_t **);
int            extend_lr_bp   (uint8_t *, int, int, int, fmdpos_t *, index_t *);
int            seed_thr       (uint8_t *, int, int, int, index_t *, seedstack_t **);
int            extend_lr_thr  (uint8_t *, int *, int, int, fmdpos_t *, index_t *);
int            seed_interv    (uint8_t *, int, int, int, int, index_t *, seedstack_t **);
int            seed_mismatch  (char *, size_t, int, seedstack_t **, index_t *);
int            find_mismatch  (char *, size_t, int, seedstack_t **, index_t *);
int            seed_the_right_way (char *, seedstack_t **, index_t *);
int            seed_block_all (char *, size_t, seedstack_t **, index_t *, int, int);
int            seed_wings     (char *, int32_t, int32_t, seedstack_t **, seedopt_t, index_t *);
int            get_smem       (seedstack_t **, seedstack_t *);
hit_t        * compute_hits   (seedstack_t * seeds, index_t * index, uint64_t * hit_cnt);
hit_t        * chain_seeds    (seedstack_t *, int, index_t *, int, double, size_t *);
int            find_repeats   (seedstack_t * seeds, matchlist_t ** repeats, int32_t repeat_thr);
// Seedstack functions.
seedstack_t  * seedstack_new    (int);
int            seedstack_push   (seed_t, seedstack_t **);
int            seedstack_append (seedstack_t **, seedstack_t *);
// Mergesort compar functions.
int            compar_hit_locus  (const void * a, const void * b, const int param);
int            compar_seed_start (const void * a, const void * b, const int param);
int            compar_hit_locus_qsort (const void * a, const void * b);
int            compar_hit_errors_qsort (const void * a, const void * b);
// Aux functions.
char *         rev_comp         (char * seq);

int
count_mismatch
(
 char         * seq,
 size_t         slen,
 index_t      * index
 );


#endif
