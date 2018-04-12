#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#ifndef _INDEX_BWT_H
#define _INDEX_BWT_H

#define BWT_QUERY_PREFIX 0
#define BWT_QUERY_SUFFIX 1

// Typedef structures.
typedef struct bwt_t      bwt_t;
typedef struct bwtquery_t bwtquery_t;

// Type interface.
struct bwt_t;
struct bwtquery_t;

// bwtquery_t function interface.
bwtquery_t  *  bwt_new_query     (bwt_t * bwt);
bwtquery_t  *  bwt_dup_query     (bwtquery_t * q);
bwtquery_t **  bwt_new_vec       (bwt_t * bwt);
bwtquery_t **  bwt_dup_vec       (bwtquery_t ** qv);
int            bwt_free_vec      (bwtquery_t ** qv);

// Index query function interface.
int           bwt_query         (int sym, int end, bwtquery_t * q);
int           bwt_query_all     (int end, bwtquery_t * q, bwtquery_t ** qv);
int           bwt_prefix        (int sym, bwtquery_t * q);
int           bwt_prefix_all    (bwtquery_t * q, bwtquery_t ** qv);

#endif
