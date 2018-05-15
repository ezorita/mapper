#define _POSIX_C_SOURCE  200809L
#include <stdio.h>
#include <string.h>
#include <glob.h>
#include "index_sym.h"
#include "index_txt.h"
#include "index_sar.h"
#include "index_bwt.h"
#include "index_ann.h"

#ifndef _INDEX_H
#define _INDEX_H

#define INDEX_FASTA_BUFFER 100

typedef struct index_t index_t;

struct index_t {
   char     * fname_base;
   sym_t    * sym;
   txt_t    * txt;
   sar_t    * sar;
   bwt_t    * bwt;
   int32_t    ann_cnt;
   ann_t   ** ann;
};


index_t   * index_read     (char * filename_base);
index_t   * index_build    (char * ref_txt, char * filename_base);
int         index_ann_new  (int kmer, int tau, int threads, index_t * index);
void        index_free     (index_t * index);

#endif
