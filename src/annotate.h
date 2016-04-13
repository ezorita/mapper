#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "bwtquery.h"
#include "algs.h"
#include "blocksearch.h"

#ifndef _ANNOTATE_H
#define _ANNOTATE_H

#define COMPUTE_INDELS 0
#define STORE_ANNOTATION 1
#define STORE_SEEDTABLE  2

typedef struct {
   htable_t * htable;
   uint8_t  * bitfield;
   size_t     ann_size;
   uint64_t   ann_set;
   uint64_t   sht_coll;
   uint64_t   sht_set;
} annotation_t;

typedef struct {
   int               kmer;
   int               tau;
   int               seed_tau;
   int               repeat_thr;
   int               mode;
   uint64_t          beg;
   uint64_t          end;
   uint64_t        * unique;
   uint64_t        * computed;
   uint64_t        * collision;
   uint64_t        * htable_pos;
   uint8_t         * repeat_bf;
   int             * done;
   index_t         * index;
   pthread_mutex_t * mutex;
   pthread_mutex_t * ht_mutex;
   pthread_cond_t  * monitor;
   uint64_t        * kmers;
   htable_t        * htable;
} annjob_t;

int             reverse_duplicate (uint8_t *, int);
int             contains_n        (uint8_t *, int);
annotation_t    annotate          (int, int, int, int, index_t *,int, int);
void          * annotate_mt       (void*);


#endif
