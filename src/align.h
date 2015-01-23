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

#define ALIGN_WIDTH         100
#define MIN_ALIGNMENT_LEN   1

// Breakpoint detection parameters
#define BREAKPOINT_THR   10.0
#define READ_MATCH_PROB  0.85
#define RAND_MATCH_PROB  0.50
#define READ_ERROR_PROB  (1.0 - READ_MATCH_PROB)
#define RAND_ERROR_PROB  (1.0 - RAND_MATCH_PROB)

// Dynamic allocation
#define ALLOC_BLOCK_SIZE 100

#define align_max(a,b) ((a) > (b) ? (a) : (b))
#define align_min(a,b) ((a) < (b) ? (a) : (b))

typedef struct align_t align_t;

struct align_t {
   long start;
   long end;
   long score;
   long pathlen;
};

align_t      nw_align         (char * query, char * ref, int len_q, int dir_q, int dir_r);
align_t      sw_align         (char * read, int rdlen, char * ref, int rflen);

#endif
