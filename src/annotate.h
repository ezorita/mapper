#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "index.h"
#include "indexquery.h"
#include "algs.h"
#include "xxhash.h"
#include "blocksearch.h"

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
   int               kmer;
   int               tau;
   uint64_t          beg;
   uint64_t          end;
   uint64_t        * computed;
   uint8_t         * counts;
   int             * done;
   index_t         * index;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   uint64_t        * kmers;
} annjob_t;

int         reverse_duplicate (char *, int);
int         contains_n        (char *, int);
int         annotate          (int, int, index_t *,int);
void *      annotate_mt       (void*);
htable_t *  htable_new        (int, int, size_t);
int         htable_insert     (uint8_t*, uint32_t, uint8_t, uint32_t, htable_t*, pthread_mutex_t*);
