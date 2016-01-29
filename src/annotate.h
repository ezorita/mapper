#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "index.h"
#include "indexquery.h"
#include "algs.h"
#include "xxhash.h"
#include "blocksearch.h"

#define COMPUTE_INDELS 0
/*
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
*/
typedef struct {
   uint64_t mask;
   uint8_t  bits;
   uint8_t  table[];
} htable_t;

typedef struct {
   int               kmer;
   int               tau;
   int               repeat_thr;
   uint64_t          beg;
   uint64_t          end;
   uint64_t        * computed;
   uint64_t        * collision;
   uint8_t         * repeat_bf;
   int             * done;
   index_t         * index;
   pthread_mutex_t * mutex;
   pthread_mutex_t * ht_mutex;
   pthread_cond_t  * monitor;
   uint64_t        * kmers;
   htable_t        * htable;
} annjob_t;

int         reverse_duplicate (uint8_t *, int);
int         contains_n        (uint8_t *, int);
int         annotate          (int, int, int, index_t *,int);
void *      annotate_mt       (void*);
htable_t *  htable_new        (uint8_t);
/*
htable_t *  htable_new        (int, int, size_t);
int         htable_insert     (uint8_t*, uint32_t, uint8_t, uint32_t, htable_t*, pthread_mutex_t*);
*/
