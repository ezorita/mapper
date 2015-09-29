#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "algs.h"
#include "align.h"
#include "index.h"

#ifndef _FILTER_H
#define _FILTER_H

typedef struct filteropt_t filteropt_t;

struct filteropt_t {
   long dist_accept;
   long max_align_per_read;
   double read_ref_ratio;
   double align_accept_eexp;
   double overlap_max_tolerance;
   double align_seed_filter_thr;
   double align_filter_ident;
   double align_filter_eexp;
   double mapq_evalue_ratio;
};


// Aligning
int             align_seeds      (char *, matchlist_t *, matchlist_t **, index_t *, filteropt_t, alignopt_t);
double          e_value          (int L, int m, long gsize);

// Post-processing functions.
matchlist_t **  merge_intervals  (matchlist_t *, double, int32_t *);
int             compute_mapq     (matchlist_t **, int, double, seq_t, index_t *);

// Compar functions.
int             compar_seedhits  (const void * a, const void * b, const int param);
int             compar_matcheexp (const void * a, const void * b, const int param);
int             compar_intvstart (const void * a, const void * b, const int param);

#endif
