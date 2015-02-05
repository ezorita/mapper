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
#define HITBUF_SIZE   32
#define REPEATS_SIZE  4
#define MNODE_SIZE    4

// Output format parameters
#define PRINT_REPEATS_NUM 5

// Search algorithm parameters
#define KMER_SIZE      22
#define LAST_THRESHOLD 7
#define HIT_MAX_LOCI   20

// Hitmap analysis parameters
#define READ_TO_GENOME_RATIO 2
#define MIN_DISTANCE_ACCEPT  10
#define MATCHLIST_SIZE       10000

// Sequence quality parameters
#define SEQ_MINLEN     75
#define SEQ_MINID      0.7
#define INTERVAL_MINID 0.8

// Read assembly parameters
#define REPEAT_OVERLAP    0.9
#define OVERLAP_THR       20
#define CONTIGS_OVERLAP   0.5

// Sequence ID definitions.
#define KMERID_BITS     24
#define KMERID_MASK     0x0000000000FFFFFF

// Flags definition.
#define WARNING_OVERLAP 0x00000001
#define FLAG_FUSED      0x00000002

// Macros
#define seqid_read(a)  ((a >> SEQID_BITS) & SEQID_MASK)
#define seqid_kmer(a)  (a & SEQID_MASK)
#define hm_max(a,b)    ((a) > (b) ? (a) : (b))
#define hm_min(a,b)    ((a) < (b) ? (a) : (b))

// Private structure definition

typedef struct sub_t       sub_t;
typedef struct match_t     match_t;
typedef struct matchlist_t matchlist_t;
typedef struct sublist_t   sublist_t;
typedef struct mnode_t     mnode_t;
typedef struct hmargs_t    hmargs_t;

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
   int dir;
   int flags; // This fills the alingment. Maybe for future use.
   double ident;
   matchlist_t * repeats;
};

struct matchlist_t {
          int       size;
          int       pos;
   struct match_t * match[];
};

struct sublist_t {
          int   size;
   struct sub_t sub[];
};

struct hmargs_t {
   int verbose;
};



int           hitmap           (int tau, index_t * index, chr_t * chr, seqstack_t * seqs, hmargs_t args);
int           poucet_search    (sublist_t * subseqs, pstack_t ** pebbles, pstack_t ** hits, trie_t ** trie, index_t * index, int tau, int kmer_size, int max_trail, int verbose);
int           hitmap_analysis  (vstack_t * hitmap, matchlist_t * loci, int kmer_size, int maxdist);
int           map_hits         (pstack_t ** hits, vstack_t ** hitmap, index_t * index, int tau, int id);
int           align_seeds      (seq_t seq, matchlist_t * seeds, matchlist_t ** seqmatches, index_t * index);
int           hitmap_push      (vstack_t ** hitmap, index_t * index, long * fm_ptr, int id);
sublist_t   * process_subseq   (seq_t * seqs, int numseqs, int k, vstack_t ** hitmaps);
void          mergesort_match  (match_t * data, match_t * aux, int size, int b);
int           subseqsort       (sub_t * data, int numels, int slen, int thrmax);
void          fuse_matches     (matchlist_t ** listp, int slen);
int           find_repeats     (matchlist_t * list);
matchlist_t * combine_matches  (matchlist_t * list);
int           feedback_gaps    (int kmer_size, int seqnum, seq_t seq, matchlist_t * intervals, sublist_t * subseqs, vstack_t ** hitmap);
void          print_intervals  (matchlist_t * intervals, chr_t * chr, int max_repeats);
int           fill_gaps        (matchlist_t ** intervp, matchlist_t * matches, double contigs_overlap);
matchlist_t * matchlist_new    (int elements);
int           matchlist_add    (matchlist_t ** listp, match_t * match);

// mergesort_mt compar functions.
int           compar_seqsort   (const void * a, const void * b, const int val);
int           compar_matchlen  (const void * a, const void * b, const int param);
int           compar_matchid   (const void * a, const void * b, const int param);
int           compar_matchstart(const void * a, const void * b, const int param);
int           compar_matchend  (const void * a, const void * b, const int param);
int           compar_matchsize (const void * a, const void * b, const int param);
int           compar_refstart  (const void * a, const void * b, const int param);

#endif
