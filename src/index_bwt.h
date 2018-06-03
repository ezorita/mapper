#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "errhandler.h"
#include "index_sar.h"

#ifndef _INDEX_BWT_H
#define _INDEX_BWT_H

/* Magic number:
** bits 32-63 (mapper): 0x0fcb0fcb
** bits 16-31 (index_bwt): 0x0003
** bits 00-15 (file type/version): 0x001
*/ 
#define BWT_FILE_MAGICNO      0x0fcb0fcb00030001

// Parameter macros.
#define BWT_QUERY_PREFIX 0
#define BWT_QUERY_SUFFIX 1

// Default BWT structure parameters.
#define BWT_OCC_MARK_INTV_DEF 14
#define BWT_OCC_WORD_SIZE_DEF 64
#define BWT_OCC_MARK_BITS_DEF (BWT_OCC_MARK_INTV_DEF * BWT_OCC_WORD_SIZE_DEF)

// Typedef structures.
typedef struct bwt_t      bwt_t;
typedef struct bwtquery_t bwtquery_t;

// Type interface.
struct bwt_t;
struct bwtquery_t;

// bwtquery_t function interface.
bwtquery_t  *  bwt_new_query   (bwt_t * bwt);
bwtquery_t  *  bwt_dup_query   (bwtquery_t * q);
bwtquery_t **  bwt_new_vec     (bwt_t * bwt);
bwtquery_t **  bwt_dup_vec     (bwtquery_t ** qv);
int            bwt_free_vec    (bwtquery_t ** qv);

// Index query function interface.
int            bwt_query       (int sym, int end, bwtquery_t * q, bwtquery_t * qo);
int            bwt_query_all   (int end, bwtquery_t * q, bwtquery_t ** qv);
int            bwt_prefix      (int sym, bwtquery_t * q, bwtquery_t * qo);
int            bwt_prefix_all  (bwtquery_t * q, bwtquery_t ** qv);

// bwt build functions.
bwt_t       *  bwt_build       (txt_t * txt, sar_t * sar);
bwt_t       *  bwt_build_opt   (txt_t * txt, sar_t * sar, uint64_t mark_intv);
void           bwt_free        (bwt_t * bwt);

// Helper functions.
txt_t       *  bwt_get_text    (bwt_t * bwt);
int64_t        bwt_start       (bwtquery_t * q);
int64_t        bwt_size        (bwtquery_t * q);
int64_t        bwt_depth       (bwtquery_t * q);
int64_t        bwt_rcstart     (bwtquery_t * q);
bwt_t       *  bwt_get_bwt     (bwtquery_t * q);

// I/O functions.
int            bwt_file_write  (char * filename, bwt_t * bwt);
bwt_t       *  bwt_file_read   (char * filename, txt_t * txt);

#endif
