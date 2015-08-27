#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "definitions.h"

#define LUT_KMER_SIZE     9
#define OCC_MARK_INTERVAL 14
#define OCC_WORD_SIZE     64
#define OCC_MARK_BITS     (OCC_MARK_INTERVAL * OCC_WORD_SIZE)
#define LCP_MARK_INTERVAL 16
#define LCP_WORD_SIZE     64
#define LCP_MARK_BITS     (LCP_MARK_INTERVAL * LCP_WORD_SIZE)
#define LCP_MIN_DEPTH     15

#define SSV_DIR_FWD       1
#define SSV_DIR_BWD       -1

typedef struct chr_t      chr_t;
typedef struct index_t    index_t;
typedef struct bwpos_t    bwpos_t;
typedef struct fmdpos_t   fmdpos_t;
typedef struct lcpval_t   lcpval_t;
typedef struct lcpdata_t  lcpdata_t;
typedef struct list32_t   list32_t;

struct lcpval_t {
   uint8_t lcp;
   int8_t  offset;
};

struct lcpdata_t {
   uint64_t   size;
   lcpval_t * lcp;
};

struct list32_t {
   uint64_t   size;
   int32_t  * val;
};

struct chr_t {
   int     nchr;
   long  * start;
   char ** name;
};

struct index_t {
   // Gen file.
   uint64_t    size;
   char      * genome;
   // OCC.
   uint64_t  * c;
   uint64_t    occ_mark_int;
   uint64_t  * occ;
   // Suffix Array.
   uint64_t    sa_bits;
   uint64_t  * sa;
   // Lookup table.
   uint64_t    lut_kmer;
   uint64_t  * lut;
   // LCP index.
   uint64_t    lcp_mark_int;
   uint64_t    lcp_min_depth;
   uint64_t  * lcp_sample_idx;
   uint64_t  * lcp_extend_idx;
   lcpdata_t * lcp_sample;
   list32_t  * lcp_extend;
   // Chromosome index.
   chr_t     * chr;
};

struct bwpos_t {
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
bwpos_t  suffix_extend     (int nt, bwpos_t pos, index_t * index);
int      suffix_shrink     (bwpos_t pos, bwpos_t * newpos, index_t * index);
int      suffix_ssv_search (uint64_t pos, bwpos_t * newpos, index_t * index);
int      suffix_ssv        (bwpos_t pos, bwpos_t * newpos, index_t * index);
