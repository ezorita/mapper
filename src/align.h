#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _ALIGN_H
#define _ALIGN_H

// SW-alignment
#define SCORE_MATCH    1
#define SCORE_DELETE   -1
#define SCORE_INSERT   -1

// Align direction
#define ALIGN_FORWARD   1
#define ALIGN_BACKWARD -1

// Breakpoint detection parameters
#define BREAKPOINT_THR   10.0
#define READ_MATCH_PROB  0.75
#define RAND_MATCH_PROB  0.50
#define READ_ERROR_PROB  (1.0 - READ_MATCH_PROB)
#define RAND_ERROR_PROB  (1.0 - RAND_MATCH_PROB)

// Dynamic allocation
#define ALLOC_BLOCK_SIZE 100

#define align_max(a,b) ((a) > (b) ? (a) : (b))
#define align_min(a,b) ((a) < (b) ? (a) : (b))

typedef struct align_t align_t;
typedef struct alignopt_t alignopt_t;
struct align_t {
   long start;
   long end;
   long score;
   long pathlen;
};

struct alignopt_t {
   double border_slope;
   double border_y0;
   double bp_thr;
   int    bp_period;
   int    bp_repeats;
   double read_error;
   double rand_error;
   double read_subst;
   double rand_subst;

};

align_t nw_align (const char * query, const char * ref, const int len_q, const int dir_q, int align_min, const int dir_r, const alignopt_t opt);

#endif
