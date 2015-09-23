#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <execinfo.h>
#include <signal.h>
#include "algs.h"
#include "seed.h"
#include "indexquery.h"
#include "align.h"

#ifndef _MAPPER_H
#define _MAPPER_H

#define BUFFER_SIZE   100
#define CHRSTR_SIZE   50
#define SEQSTACK_SIZE 1024

#define REPEATS_SIZE 16

#define map_min(x,y) ((x) < (y) ? (x) : (y))
#define map_max(x,y) ((x) > (y) ? (x) : (y))

// Flags definition.
#define WARNING_OVERLAP 0x00000001
#define FLAG_FUSED      0x00000002


typedef struct idxfiles_t  idxfiles_t;
typedef struct hit_t       hit_t;
typedef struct match_t     match_t;
typedef struct matchlist_t matchlist_t;
typedef struct mapopt_t    mapopt_t;

// Input types.

typedef enum {
   FASTA,
   FASTQ,
   RAW
} format_t;


struct idxfiles_t {
   // Original file pointers.
   void      * gen_file;
   void      * occ_file;
   void      * sa_file;
   void      * lcp_file;
//   void      * lut_file;
   chr_t     * chr;
};

struct hit_t {
   uint64_t locus;
   uint32_t qrypos;
   uint16_t depth;
   uint16_t bulk;
};

struct match_t {
   long ref_s;
   long ref_e;
   int read_s;
   int read_e;
   int hits;
   int score;
   int flags;
   int interval;
   double ident;
   double e_exp;
   matchlist_t * repeats;
};

struct matchlist_t {
   int       size;
   int       pos;
   match_t   match[];
};

struct mapopt_t {
   long dist_accept;
   long max_align_per_read;
   double read_ref_ratio;
   double align_accept_eexp;
   double overlap_max_tolerance;
   double align_seed_filter_thr;
   double align_filter_ident;
   double align_filter_eexp;
   alignopt_t align;
   seedopt_t seed;
};

// File/Index management.
idxfiles_t    * index_open      (char * file);
index_t       * index_format    (idxfiles_t * files);
chr_t         * read_CHRindex   (char * filename);
seqstack_t    * read_file       (FILE * inputf);

// Seeding
hit_t         * compute_hits    (seedstack_t * seeds, index_t * index, uint64_t * hit_cnt);
int             match_seeds     (seedstack_t *, int, matchlist_t *, index_t *, int, double);
matchlist_t   * matchlist_new    (int elements);
int             matchlist_add    (matchlist_t ** listp, match_t match);

// Aligning
int             align_seeds      (char *, matchlist_t *, matchlist_t **, index_t *, mapopt_t);
int             align_simple     (char *, matchlist_t *, matchlist_t **, index_t *, mapopt_t);
double          e_value          (int L, int m, long gsize);

// Post-processing functions.
matchlist_t **  merge_intervals  (matchlist_t *, double, int32_t *);

// Compar functions.
int             compar_hit_locus (const void * a, const void * b, const int param);
int             compar_seedhits  (const void * a, const void * b, const int param);
int             compar_matcheexp (const void * a, const void * b, const int param);
int             compar_intvstart (const void * a, const void * b, const int param);

#endif
