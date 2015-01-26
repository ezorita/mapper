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
#define KMER_SIZE      22
#define LAST_THRESHOLD 7
#define HIT_MAX_LOCI   20

// Hitmap analysis parameters
#define READ_TO_GENOME_RATIO 2
#define MATCHLIST_SIZE       1000

// Sequence quality parameters
#define SEQ_MINLEN     100
#define SEQ_MINID      0.6

// Sequence ID definitions.
#define KMERID_BITS     24
#define KMERID_MASK     0x0000000000FFFFFF

// Macros
#define seqid_read(a)  ((a >> SEQID_BITS) & SEQID_MASK)
#define seqid_kmer(a)  (a & SEQID_MASK)

// Private structure definition

typedef struct sub_t       sub_t;
typedef struct match_t     match_t;
typedef struct matchlist_t matchlist_t;
typedef struct sublist_t   sublist_t;

struct sub_t {
   long               seqid;
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
   double ident;
   int dir;
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
int           hitmap_push      (vstack_t ** hitmap, index_t * index, long * fm_ptr, int id);
sublist_t   * process_subseq   (seq_t * seqs, int numseqs, int k, vstack_t ** hitmaps);
void          mergesort_match  (match_t * data, match_t * aux, int size, int b);
int           subseqsort       (sub_t * data, int numels, int slen, int thrmax);
int           compar_seqsort   (const void * a, const void * b, const int val);
int           compar_matchlen  (const void * a, const void * b, const int param);
int           compar_matchid   (const void * a, const void * b, const int param);
int           compar_matchstart(const void * a, const void * b, const int param);

#endif
