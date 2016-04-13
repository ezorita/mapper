#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "definitions.h"

#ifndef _BWTQUERY_H
#define _BWTQUERY_H

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

uint64_t get_sa            (uint64_t pos, sar_t * sa);
int      get_sa_range      (uint64_t start, uint64_t size, uint64_t * out, sar_t * sa);
int      get_occ           (int64_t pos, uint64_t * occ, int64_t * val);
uint64_t get_occ_nt        (int64_t pos, uint64_t * occ, int nt);
fmdpos_t extend_bw         (int nt, fmdpos_t pos, bwt_t * bwt);
fmdpos_t extend_fw         (int nt, fmdpos_t pos, bwt_t * bwt);
int      extend_bw_all     (fmdpos_t pos, fmdpos_t * newpos, bwt_t * bwt);
int      extend_fw_all     (fmdpos_t pos, fmdpos_t * newpos, bwt_t * bwt);
int      suffix_extend     (int nt, bwpos_t pos, bwpos_t * newpos, bwt_t * bwt);
int      suffix_extend_all (bwpos_t pos, bwpos_t * newpos, bwt_t * bwt);
int      suffix_string     (char * suf, int slen, uint64_t minloci, bwpos_t * newpos, bwt_t * index);

#endif
