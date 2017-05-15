#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "bwtquery.h"
#include "algs.h"
#include "blocksearch.h"

#ifndef _ANNOTATE_H
#define _ANNOTATE_H

#define NO_INFO 0xFFFF

typedef struct {
   uint8_t    kmer;
   uint8_t    tau;
   uint8_t  * info;
   size_t     size;
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

// Annotation core functions.
annotation_t    annotate                    (int, int, index_t*, int);
void          * annotate_mt                 (void*);

// Multithreading scheduler functions.
void            job_ranges_rec              (fmdpos_t, int, int, int, int, int*, annjob_t*, index_t*);

// Thread helper functions.
int             next_seq                    (int32_t, int32_t, int64_t, int64_t*, uint8_t*, fmdpos_t*, index_t*);
void            store_hits                  (annjob_t*, pathstack_t*, fmdpos_t);
void            merge_alignments            (uint8_t *, uint8_t *, int);
void            compute_aln_positions       (uint64_t *, uint8_t  *, int, int, int);
void            send_progress_to_scheduler  (int64_t, int64_t*, annjob_t*);

#endif
