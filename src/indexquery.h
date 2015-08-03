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

typedef struct chr_t    chr_t;
typedef struct index_t  index_t;
typedef struct bwpos_t  bwpos_t;
typedef struct fmdpos_t fmdpos_t;



struct chr_t {
   int     nchr;
   long  * start;
   char ** name;
};

struct index_t {
   void     * gen_file;
   void     * occ_file;
   void     * sa_file;
   uint64_t   size;
   uint64_t * c;
   char     * genome;
   uint64_t   sa_bits;
   uint64_t * sa;
   uint64_t * occ;
   chr_t    * chr;
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


uint64_t get_sa     (uint64_t pos, uint64_t * sa, int bits);
int      get_occ    (int64_t pos, uint64_t * occ, int64_t * val);
uint64_t get_occ_nt (int64_t pos, uint64_t * occ, int nt);
fmdpos_t extend_bw  (int nt, fmdpos_t pos, index_t * index);
fmdpos_t extend_fw  (int nt, fmdpos_t pos, index_t * index);
bwpos_t  bw_search  (int nt, bwpos_t pos, index_t * index);
bwpos_t  bw_shrink  (bwpos_t pos, index_t * index);
