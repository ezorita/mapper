#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "definitions.h"
#include "index.h"

#ifndef _INDEXQUERY_H
#define _INDEXQUERY_H

#define LUT_KMER_SIZE     9
#define OCC_MARK_INTERVAL 14
#define OCC_WORD_SIZE     64
#define OCC_MARK_BITS     (OCC_MARK_INTERVAL * OCC_WORD_SIZE)
#define LCP_MARK_INTERVAL 16
#define LCP_WORD_SIZE     64
#define LCP_MARK_BITS     (LCP_MARK_INTERVAL * LCP_WORD_SIZE)
#define LCP_MIN_DEPTH     0

#define SSV_DIR_FWD       1
#define SSV_DIR_BWD       -1

typedef struct bwpos_t    bwpos_t;
typedef struct fmdpos_t   fmdpos_t;

struct bwpos_t {
   int32_t depth;
   int64_t sp;
   int64_t ep;
};

struct fmdpos_t {
   int64_t fp;
   int64_t rp;
   int64_t  sz;
};


uint64_t get_sa            (uint64_t pos, uint64_t * sa, int bits);
int      get_occ           (int64_t pos, uint64_t * occ, int64_t * val);
uint64_t get_occ_nt        (int64_t pos, uint64_t * occ, int nt);
fmdpos_t extend_bw         (int nt, fmdpos_t pos, index_t * index);
fmdpos_t extend_fw         (int nt, fmdpos_t pos, index_t * index);
int      suffix_extend     (int nt, bwpos_t pos, bwpos_t * newpos, index_t * index);
int      suffix_shrink     (bwpos_t pos, bwpos_t * newpos, index_t * index);
int      suffix_ssv_search (uint64_t pos, bwpos_t * newpos, index_t * index);
int      suffix_ssv        (bwpos_t pos, bwpos_t * newpos, index_t * index);
int      suffix_string     (char * suf, int slen, uint64_t minloci, bwpos_t * newpos, index_t * index);

// Deprecated functions.
//int      suffix_lut        (char * suf, int slen, bwpos_t * newpos, index_t * index);

#endif
