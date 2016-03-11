#include <string.h>
#include <stdint.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef _ALIGN_H
#define _ALIGN_H

// Algorithm select.
#define MAX_NAIVE_WIDTH 8

// BF-alignment
#define WORD_SIZE 64
#define MASK_MSB 0x8000000000000000

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

typedef struct align_t    align_t;
typedef struct alignopt_t alignopt_t;
typedef struct path_t     path_t;

struct align_t {
   long start;
   long end;
   long score;
   long pathlen;
};

struct alignopt_t {
   double bp_thr;
   double bp_max_thr;
   int    bp_diagonal;
   int    bp_period;
   int    bp_resolution;
   int    bp_repeats;
   double read_error;
   double rand_error;
   double width_ratio;
   int    width_min;
   double logAe;
   double logBe;
};

struct path_t {
   uint32_t score;
   uint32_t row;
   uint32_t col;
   uint32_t max_indel;
};

int    align        (int qlen,char* query,int rlen,char* ref,int dir_r,int dir_q,int max,alignopt_t o,path_t * path);
int    dbf_align    (int qlen,char* query,int rlen,char* ref,int dir_r,int dir_q,int max,alignopt_t o,path_t * path);
int    naive_align  (int qlen,char* query,int rlen,char* ref,int dir_r,int dir_q,int max,alignopt_t o,path_t * path);
path_t align_bp       (int qlen,char* query,int rlen,char* ref,int min_len,int dir_r,int dir_q,alignopt_t opt);
path_t dbf_align_bp   (int qlen,char* query,int rlen,char* ref,int min_len,int dir_r,int dir_q,alignopt_t opt);
path_t naive_align_bp (int qlen,char* query,int rlen,char* ref,int min_len,int dir_r,int dir_q,alignopt_t opt);
int    mapq_align   (const char*, const char*, const char*, const int, const int, int*);

#endif
