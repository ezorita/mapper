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
   int flags;
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
   int maxtau;
   int verbose;
   int repeat_print_num;
   int kmer_size;
   int seed_max_loci;
   int read_ref_ratio;
   int dist_accept;
   int max_align_per_read;
   int match_min_len;
   double match_min_id;
   double feedback_id_thr;
   double repeat_min_overlap;
   double overlap_tolerance;
   double overlap_max_tolerance;
   double align_likelihood_thr;
   double read_match_prob;
   double rand_match_prob;
};



int           hitmap           (index_t * index, chr_t * chr, seqstack_t * seqs, hmargs_t args);
int           poucet_search    (sublist_t * subseqs, pstack_t ** pebbles, pstack_t ** hits, trie_t ** trie, index_t * index, int tau, int kmer_size, int max_trail, int max_loci_per_hit, int verbose);
int           hitmap_analysis  (vstack_t * hitmap, matchlist_t * loci, int kmer_size, int maxdist, hmargs_t hmargs);
int           map_hits         (pstack_t ** hits, vstack_t ** hitmap, index_t * index, int tau, int id, int max_loci);
int           align_seeds      (seq_t seq, matchlist_t * seeds, matchlist_t ** seqmatches, index_t * index, hmargs_t hmargs);
sublist_t   * process_subseq   (seq_t * seqs, int numseqs, int k, vstack_t ** hitmaps);
void          mergesort_match  (match_t * data, match_t * aux, int size, int b);
int           subseqsort       (sub_t * data, int numels, int slen, int thrmax);
void          fuse_matches     (matchlist_t ** listp, int slen, hmargs_t hmargs);
int           find_repeats     (matchlist_t * list, double overlap);
matchlist_t * combine_matches  (matchlist_t * list, double overlap_tolerance);
int           feedback_gaps    (int kmer_size, int seqnum, seq_t seq, matchlist_t * intervals, sublist_t * subseqs, vstack_t ** hitmap, hmargs_t hmargs);
void          print_intervals  (matchlist_t * intervals, chr_t * chr, int max_repeats);
int           fill_gaps        (matchlist_t ** intervp, matchlist_t * matches, int seq_len, int seq_minlen, double gap_coverage, double max_overlap);
matchlist_t * matchlist_new    (int elements);
int           matchlist_add    (matchlist_t ** listp, match_t * match);
void          free_match       (match_t * match);

// mergesort_mt compar functions.
int           compar_seqsort   (const void * a, const void * b, const int val);
int           compar_matchid   (const void * a, const void * b, const int param);
int           compar_readstart (const void * a, const void * b, const int param);
int           compar_readend   (const void * a, const void * b, const int param);
int           compar_matchspan (const void * a, const void * b, const int param);
int           compar_refstart  (const void * a, const void * b, const int param);

#endif
