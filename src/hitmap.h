#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "algs.h"
#include "poucet.h"
#include "align.h"

#ifndef _HITMAP_H
#define _HITMAP_H

// Default size params
#define HITMAP_SIZE   1024
#define HITSTACK_SIZE 16
#define SORTBUF_SIZE  4096
#define PBSTACK_SIZE  1024

// Search algorithm parameters
#define KMER_SIZE      17
#define LAST_THRESHOLD 7
#define SUBSEQID_BITS  24
#define HIT_MAX_LOCI   1000

#define WINDOW_SIZE    200

#define SCORE_BITS     16
#define SCORE_MASK     0x000000000000FFFF

// Private structure definition

typedef struct sub_t      sub_t;
typedef struct match_t    match_t;
typedef struct matchlist_t matchlist_t;
typedef struct sublist_t  sublist_t;

struct sub_t {
   char             * seq;
   struct vstack_t ** hitmap;
};

struct match_t {
   long ref_s;
   long ref_e;
   int read_s;
   int read_e;
   int hits;
   int score;
};

struct matchlist_t {
          int    size;
          int    pos;
   struct match_t  match[];
};

struct sublist_t {
          int   size;
   struct sub_t sub[];
};


int           hitmap           (int tau, index_t index, chr_t * chr, seqstack_t * seqs);
int           hitmap_analysis  (vstack_t * hitmap, matchlist_t * loci, int kmer_size, int tau, int maxdist);
int           map_hits         (pstack_t ** hits, vstack_t ** hitmap, index_t * index, int tau, int id);
sublist_t   * process_subseq   (seq_t * seqs, int numseqs, int k, vstack_t ** hitmaps);
void          mergesort_match  (match_t * data, match_t * aux, int size, int b);
int           subseqsort       (sub_t * data, int numels, int slen, int thrmax);
int           compar_seqsort   (const void * a, const void * b, const int val);


#endif
