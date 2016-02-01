#define _GNU_SOURCE
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include "definitions.h"

#ifndef _MAPPER_ALGS_H
#define _MAPPER_ALGS_H

// Flags definition.
#define WARNING_OVERLAP 0x00000001
#define FLAG_REPEAT     0x00000002


// Structure typedef.
typedef struct index_t     index_t;
typedef struct seq_t       seq_t;
typedef struct list_t      list_t;
typedef struct chr_t       chr_t;
typedef struct node_t      node_t;
typedef struct sortargs_t  sortargs_t;
typedef struct seqstack_t  seqstack_t;
typedef struct vstack_t    vstack_t;
typedef struct match_t     match_t;
typedef struct matchlist_t matchlist_t;
typedef struct htable_t    htable_t;

// Macros
#define min(x,y) ((x) < (y) ? (x) : (y))
#define max(x,y) ((x) > (y) ? (x) : (y))

// Data types
typedef unsigned int uint;

// Data structures

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
   double mapq;
};

struct matchlist_t {
   int       size;
   int       pos;
   match_t   match[];
};

struct seq_t {
   char * tag;
   char * seq;
   char * q;
};

struct list_t {
   long   max;
   long * val;
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

// Hash tables.

struct htable_t {
   uint64_t mask;
   uint8_t  bits;
   uint8_t  table[];
};

// Parameter structs.

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
long          bisect_search    (long start, long end, long* set, long value);
double        binom            (int l, int k);

// Stack functions.
vstack_t    * new_stack        (long size);
int           push             (vstack_t ** stackp, long value);
int           pushvec          (vstack_t ** stackp, long * vector, int vecsize);
seqstack_t  * new_seqstack     (int size);
int           seq_push         (seqstack_t ** stackp, const char* tag, const char* seq, const char* q);

// 2-bit Hash table.
htable_t    * htable_new       (uint8_t);
int           htable_get       (htable_t *, uint64_t);
int           htable_set       (htable_t *, uint64_t, uint8_t);

// Matchlist functions.
matchlist_t   * matchlist_new    (int elements);
int             matchlist_add    (matchlist_t ** listp, match_t match);


#endif
