#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>

#include "index_sar.h"
#include "index_bwt.h"
#include "blocksearch.h"

#ifndef _INDEX_ANN_H
#define _INDEX_ANN_H

/* Magic number:
** bits 32-63 (mapper): 0x0fcb0fcb
** bits 16-31 (index_ann): 0x0005
** bits 00-15 (file type/version): 0x001
*/ 
#define ANN_FILE_MAGICNO      0x0fcb0fcb00050001

// Annotation definitions.
#define ANN_NO_INFO 0xFFFF

// Typedef structures.
typedef struct ann_t      ann_t;
typedef struct locinfo_t  locinfo_t;

// Type interface.
struct ann_t;
struct locinfo_t {
   int32_t dist;
   int32_t neigh_cnt;
   int32_t align_cnt;
   int32_t align_pos[];
};

// Annotation build/free functions.
ann_t       *  ann_build        (int kmer, int tau, bwt_t * bwt, sar_t * sar, int threads);
void           ann_free         (ann_t * ann);

// Query functions.
locinfo_t   *  ann_query        (int64_t pos, ann_t * ann);

// I/O functions.
int            ann_file_write   (char * filename, ann_t * ann);
ann_t       *  ann_file_read    (char * filename);

#endif
