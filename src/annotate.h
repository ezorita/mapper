#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "index.h"
#include "indexquery.h"
#include "algs.h"
#include "xxhash.h"

#define COMPUTE_INDELS 0

typedef struct {
   size_t   count;
   size_t   index_bits;
   uint64_t index_bitmask;
   size_t   probe_bits;
   uint64_t probe_bitmask;
   uint8_t  probe_datamask;
   size_t   elm_size;
   uint8_t  c[];
} htable_t;

typedef struct {
   bwpos_t  pos;
   uint32_t depth;
#if COMPUTE_INDELS == 1
   uint32_t path;
   int8_t   row[2*MAXTAU+3];
#else
   int8_t   score;
#endif
} pebble_t;

typedef struct {
   size_t pos;
   size_t size;
   pebble_t pebble[];
} pstack_t;

typedef struct {
   uint8_t   * query;
   int         kmer;
   int         trail;
   int         tau;
   index_t   * index;
   pstack_t ** pebbles;
   pstack_t ** hits;
} arg_t;

typedef struct {
   int               kmer;
   int               tau;
   bwpos_t           beg;
   uint64_t          end;
   uint8_t         * seq;
   uint64_t        * computed;
   uint8_t         * counts;
   int             * done;
   bwpos_t         * cache;
   index_t         * index;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   uint64_t        * kmers;
} annjob_t;

int         ppush         (pstack_t **, pebble_t);
pstack_t *  pstack_new    (size_t);
int         seq_lcp       (char *, char *, int);
int         find_next     (bwpos_t *, uint8_t *, char *, int, bwpos_t *, index_t *);
int         reverse_duplicate (char *, int);
int         contains_n        (char *, int);
int         annotate      (int, int, index_t *,int);
void *      annotate_mt   (void*);
int         poucet        (pebble_t, arg_t);
uint64_t    poucet_mismatch (pebble_t, arg_t);
int         dash          (pebble_t, const arg_t);
uint64_t    dash_mismatch (pebble_t, const arg_t);
htable_t *  htable_new    (int, int, size_t);
int         htable_insert (uint8_t*, uint32_t, uint8_t, uint32_t, htable_t*, pthread_mutex_t*);
