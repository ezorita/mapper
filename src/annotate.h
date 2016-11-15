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
   uint8_t  * bitfield;
   size_t     ann_size;
   uint64_t   ann_set;
} annotation_t;

typedef struct {
   int               kmer;
   int               tau;
   uint64_t          beg;
   uint64_t          end;
   uint64_t        * unique;
   uint64_t        * computed;
   uint8_t         * repeat_bf;
   int             * done;
   index_t         * index;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   uint64_t        * kmers;
} annjob_t;

int             reverse_duplicate (uint8_t *, int);
int             contains_n        (uint8_t *, int);
annotation_t    annotate          (int, int, index_t *,int);
void          * annotate_mt       (void*);


#endif
