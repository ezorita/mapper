#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include "definitions.h"

#ifndef _BWM_ALGS_H
#define _BWM_ALGS_H

typedef struct index_t    index_t;
typedef struct seq_t      seq_t;
typedef struct list_t     list_t;
typedef struct chr_t      chr_t;
typedef struct pebble_t   pebble_t;
typedef struct node_t     node_t;
typedef struct trie_t     trie_t;
typedef struct arg_t      arg_t;
typedef struct sortargs_t sortargs_t;
typedef struct seqstack_t seqstack_t;
typedef struct vstack_t   vstack_t;
typedef struct pstack_t   pstack_t;

typedef struct bwarg_t    bwarg_t;
typedef struct bwpos_t    bwpos_t;

// Definitions
#define TRIE_SIZE     1024
#define TRIE_CHILDREN  3

// Data types
typedef unsigned int uint;

// Data structures

struct seq_t {
   char * tag;
   char * seq;
   char * rseq;
   char * q;
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

struct bwpos_t {
   long sp;
   long ep;
};

struct pebble_t {
   long sp;
   long ep;
   char nwrow[2*MAXTAU+3];
};

struct index_t {
   void     * gen_file;
   void     * occ_file;
   void     * sa_file;
   uint64_t   gsize;
   uint64_t * c;
   char     * genome;
   uint64_t * pos;
   uint64_t * occ;
   chr_t    * chr;
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

struct bwarg_t {
   struct index_t * index;
   struct bwpos_t * pebbles;
   struct bwpos_t * hits;
};

struct arg_t {
   char             * query;
   int                tau;
   int                trail;
   int                qlen;
   struct index_t   * index;
   struct pstack_t ** pebbles;
   struct pstack_t ** hits;
};

struct sortargs_t {
   void * buf0;
   void * buf1;
   int    size;
   size_t offset;
   int    b;
   int    thread;
   int    param;
   int    (*compar)(const void *, const void *, const int);
};


// Function headers.

// Sorting algorithms.
int           mergesort_mt     (void * data, int numels, size_t elmsize, int param, int thrmax, int (*compar)(const void*, const void*, const int));
void        * _mergesort       (void * args);
void          radix_sort       (long * a, long * b, long n, long maxval);

// General algorithms.
uint64_t      compute_occ      (uint64_t ptr, uint64_t * occ);
long          bisect_search    (long start, long end, long* set, long value);
uint64_t      compact_array    (uint64_t * array, uint64_t len, int bits);

// Stack functions.
vstack_t    * new_stack        (long size);
int           push             (vstack_t ** stackp, long value);
int           pushvec          (vstack_t ** stackp, long * vector, int vecsize);
pstack_t    * new_pstack       (long size);
int           ppush            (pstack_t ** stackp, pebble_t pebble);
seqstack_t  * new_seqstack     (int size);
int           seq_push         (seqstack_t ** stackp, const char* tag, const char* seq, const char* q, const int reverse);

// Trie functions.
trie_t      * trie_new         (int initial_size);
int           trie_getrow      (trie_t * trie, uint nodeid, int refval, int* wingsz, uint *nwrow);
uint          trie_insert      (trie_t ** triep, char * path, int pathlen);
void          trie_reset       (trie_t * trie);

#endif
