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
#include <pthread.h>
#include "dc3.h"

// Definitions
#define NUM_BASES 6
#define NUM_COMBS NUM_BASES*NUM_BASES*NUM_BASES
#define MAXSEQLEN 1024
#define MAXTAU    512
#define MAXSEQOUT 50

// Initial buffer/stack sizes
#define BUFFER_SIZE   100
#define GENOME_SIZE   100000
#define CHRSTR_SIZE   100
#define SEQSTACK_SIZE 1024
#define PBSTACK_SIZE  1024
#define HITSTACK_SIZE 16
#define HITMAP_SIZE   1024
#define DSTACK_SIZE   16
#define TRIE_SIZE     1024
#define SORTBUF_SIZE  4096

// Query read buffers
#define QUERYBUF_SIZE  100
#define MAXHEADER_SIZE 100
#define MINIMUM_STREAK 10

// Search algorithm parameters
#define MAX_TRAIL      30
#define EOS            -1
#define TRIE_CHILDREN  3

#define KMER_SIZE      14
#define LAST_THRESHOLD 7
#define SUBSEQID_BITS  24
#define HIT_MAX_LOCI   1000

#define WINDOW_SIZE    200

#define SCORE_BITS     16
#define SCORE_MASK     0x000000000000FFFF

// SW-alignment
#define SCORE_MATCH    1
#define SCORE_DELETE   -1
#define SCORE_INSERT   -1

// Inline definition
#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

char revert[256]  = {[0 ... 255] = 0, [0] = '@', [1] = 'A', [2] = 'C', [3] = 'G', [4] = 'N', [5] = 'T'};

char rcode[256] = {[0 ... 255] = 0,
                   ['a'] = 'T', ['c'] = 'G', ['g'] = 'C', ['t'] = 'A', ['u'] = 'A', ['n'] = 'N',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A', ['N'] = 'N' };

// Input types.

typedef enum {
   FASTA,
   FASTQ,
   RAW
} format_t;

// Data types

typedef unsigned int uint;
typedef struct index_t    index_t;
typedef struct seq_t      seq_t;
typedef struct sub_t      sub_t;
typedef struct list_t     list_t;
typedef struct chr_t      chr_t;
typedef struct pebble_t   pebble_t;
typedef struct match_t    match_t;
typedef struct matchlist_t matchlist_t;
typedef struct node_t     node_t;
typedef struct trie_t     trie_t;
typedef struct arg_t      arg_t;
typedef struct sortargs_t sortargs_t;
typedef struct seqstack_t seqstack_t;
typedef struct sublist_t  sublist_t;
typedef struct vstack_t   vstack_t;
typedef struct pstack_t   pstack_t;

// Data structures

struct seq_t {
   char * tag;
   char * seq;
   char * rseq;
};

struct sub_t {
   char             * seq;
   struct vstack_t ** hitmap;
};

struct list_t {
   long   max;
   long * val;
};

struct chr_t {
   int     nchr;
   long  * start;
   char ** name;
};

struct pebble_t {
   long sp;
   long ep;
   long rowid; // [bits 63..16] node id. [bits 15..0] score.
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

struct index_t {
   long            gsize;
   long          * c;
   char          * genome;
   long          * pos;
   struct list_t * occ;
};

struct node_t {
   uint child[TRIE_CHILDREN];
   uint parent;
}; 

struct trie_t {
   uint          pos;
   uint          size;
   struct node_t nodes[];
};

// Stacks

struct sublist_t {
          int   size;
   struct sub_t sub[];
};

struct seqstack_t {
   long         pos;
   long         size;
   struct seq_t seq[];
};

struct vstack_t {
   long pos;
   long size;
   long val[];
};

struct pstack_t {
   long            pos;
   long            size;
   struct pebble_t pebble[];
};

// Parameter structs.

struct arg_t {
   char             * query;
   int                tau;
   int                trail;
   int                qlen;
   struct index_t   * index;
   struct trie_t   ** triep;
   struct pstack_t ** pebbles;
   struct pstack_t ** hits;
};

struct sortargs_t {
   struct sub_t * buf0;
   struct sub_t * buf1;
   int            size;
   int            b;
   int            thread;
   int            slen;
};


// Query functions.
pebble_t      sw_align         (char * read, int rdlen, char * ref, int rflen);
int           hitmap_analysis  (vstack_t * hitmap, matchlist_t * loci, int kmer_size, int tau, int maxdist);
int           map_hits         (pstack_t ** hits, vstack_t ** hitmap, index_t * index, int tau, int id);
sublist_t   * process_subseq   (seq_t * seqs, int numseqs, int k, vstack_t ** hitmaps);
int           poucet           (const long sp, const long ep, const int wingsz, const uint* prow, const int depth, char* path, arg_t * arg);
void          dash             (long sp, long ep, const int depth, const int align, const char* path, const arg_t* arg);
int           query_index      (char* query, long gsize, long* c, long* ptr, list_t* occs);
int         * translate_query  (char* query);
long          bisect_search    (long start, long end, long* set, long value);

// Index functions.
int           format_FMindex   (long* index, index_t* fmindex);
chr_t       * read_CHRindex    (char* filename);
ssize_t       write_index      (char* filename);
char        * compact_genome   (char* filename, long* genomesize);
void          bwt_index        (char* text, long tlen, long** pos, vstack_t** occ);
long        * compute_c        (char* genome, long gsize);

// Stack functions.
vstack_t    * new_stack        (long size);
int           push             (vstack_t ** stackp, long value);
int           pushvec          (vstack_t ** stackp, long * vector, int vecsize);
pstack_t    * new_pstack       (long size);
void          ppush            (pstack_t ** stackp, pebble_t pebble);

// Trie functions.
trie_t      * trie_new         (int initial_size);
int           trie_getrow      (trie_t * trie, uint nodeid, int refval, int* wingsz, uint *nwrow);
uint          trie_insert      (trie_t ** triep, char * path, int pathlen);
void          trie_reset       (trie_t * trie);

// Misc functions.
seqstack_t  * read_file        (FILE * inputf, const int reverse, const int verbose);
seqstack_t  * new_seqstack     (int size);
int           seqsort          (sub_t * data, int numels, int slen, int thrmax);
void        * nukesort         (void * args);
void          mergesort_long   (long * data, long * aux, int size, int b);
void          mergesort_match  (match_t * data, match_t * aux, int size, int b);
void          radix_sort       (long * a, long * b, long n, long maxval);
