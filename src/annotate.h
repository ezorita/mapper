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
   uint32_t          kmer;
   uint16_t          tau;
   uint16_t          wsize;
   uint64_t          beg;
   uint64_t          end;
   uint64_t        * computed;
   int             * done;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   index_t         * index;
   uint8_t         * info;
} annjob_t;

void            job_ranges_rec    (fmdpos_t, int, int, int, int, int*, annjob_t*, index_t*);
int             next_seq          (int32_t, int32_t, int64_t, int64_t*, uint8_t*, fmdpos_t*, index_t*);
void            store_hits        (annjob_t*, pstree_t*);
annotation_t    annotate          (int, int, index_t*, int);
void          * annotate_mt       (void*);


#endif
