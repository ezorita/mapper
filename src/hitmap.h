#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "algs.h"
#include "poucet.h"
#include "align.h"

#ifndef _HITMAP_H
#define _HITMAP_H

// Seeding strategies
#define SEED_PARALLEL 0
#define SEED_SINGLE   1

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
   double e_exp;
   matchlist_t * repeats;
};

struct matchlist_t {
          int       size;
          int       pos;
   struct match_t * match[];
};

struct sublist_t {
   int size;
   int pos;
   struct sub_t sub[];
};

struct hmargs_t {
   int verbose;
   int search_rounds;
   int sequence_blocks;
   // Seeding strategy.
   int * tau;
   int * kmer_size;
   int * kmer_offset;
   char * qthr;
   int repeat_print_num;
   // Seeding options.
   int seed_max_loci;
   int seed_abs_max_loci;
   // Hitmap analysis options.
   int dist_accept;
   int max_align_per_read;
   int feedback_gap_minlen;
   double read_ref_ratio;
   // Alignment filter options.
   double align_full_seed_thr;
   double align_filter_ident; // Minimum identity to compute E-value.
   double align_filter_eexp; // Maximum log10(E) to consider an alignment acceptable.
   double align_accept_eexp; // Maximum log10(E) to directly accept an alignment and do not align again that region.
   double align_seed_filter_thr; // If an alignment with N hits has passed 'align_accept_ratio', alignments with at least N*align_seed_filter_thr hits will be performed anyway.
   // Post-processing options.
   double feedback_eexp_thr;
   double repeat_min_overlap;
   double overlap_tolerance;
   double overlap_max_tolerance;
   double fuse_min_spanratio;
   // Alignment algorithm options.
   alignopt_t align;
};



int           hitmap           (index_t * index, seqstack_t * seqs, hmargs_t args);
int           poucet_search    (sublist_t * subseqs, pstack_t ** pebbles, pstack_t ** hits, index_t * index, int tau, int kmer_size, int max_trail, int max_loci_per_hit, int abs_max_loci_per_hit, int verbose);
int           hitmap_analysis  (vstack_t * hitmap, matchlist_t * loci, int kmer_size, int maxdist, hmargs_t hmargs);
int           map_hits         (pstack_t ** hits, vstack_t ** hitmap, index_t * index, int tau, long id, int max_loci, int abs_max_loci);
int           align_seeds      (seq_t seq, matchlist_t * seeds, matchlist_t ** seqmatches, index_t * index, hmargs_t hmargs);
sublist_t   * process_subseq   (seq_t * seqs, int numseqs, int k, int offset, char qthr, vstack_t ** hitmaps);
void          fuse_matches     (matchlist_t ** listp, int slen, hmargs_t hmargs);
int           find_repeats     (matchlist_t * list, double overlap);
matchlist_t * combine_matches  (matchlist_t * list, double overlap_tolerance);
int           feedback_gaps    (int seqnum, seq_t seq, matchlist_t * intervals, sublist_t ** subseqs, vstack_t ** hitmap, hmargs_t hmargs, int next_kmer_size, int next_kmer_offset, char next_qthr);
void          print_intervals  (char * tagname, matchlist_t * intervals, index_t * index, int max_repeats);
int           fill_gaps        (matchlist_t ** intervp, matchlist_t * matches, int seq_len, double gap_coverage, double max_overlap);
double        e_value          (int L, int errors, long gsize);
matchlist_t * matchlist_new    (int elements);
int           matchlist_add    (matchlist_t ** listp, match_t * match);
void          free_match       (match_t * match);

// mergesort_mt compar functions.
int           compar_seqsort   (const void * a, const void * b, const int val);
int           compar_long      (const void * a, const void * b, const int param);
int           compar_matchid   (const void * a, const void * b, const int param);
int           compar_readstart (const void * a, const void * b, const int param);
int           compar_readend   (const void * a, const void * b, const int param);
int           compar_seedhits  (const void * a, const void * b, const int param);
int           compar_matchspan (const void * a, const void * b, const int param);
int           compar_matcheexp (const void * a, const void * b, const int param);
int           compar_refstart  (const void * a, const void * b, const int param);

#endif
