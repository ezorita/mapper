#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "bwtquery.h"
#include "definitions.h"

#ifndef _BLOCKSEARCH_H
#define _BLOCKSEARCH_H
#define PATHSTACK_DEF_SIZE 64
#define MAX_K 254
#define ALIGN_WORDS 4
#define ALIGN_WORD_SIZE 64

typedef struct {
   fmdpos_t pos;
   uint64_t align[ALIGN_WORDS];
   int      score;
} spath_t;

typedef struct {
   size_t   pos;
   size_t   size;
   spath_t  path[];
} pathstack_t;

typedef struct pstree_t pstree_t;

struct pstree_t {
   pathstack_t * stack;
   pstree_t    * next_l;
   pstree_t    * next_r;
};

// Block-search functions.
void          blocksc_trail          (uint8_t *, fmdpos_t *, int, int, int, index_t *, pstree_t *);
void          blocksearch_trail_rec  (uint8_t *, int, int, int, int, index_t *, pstree_t *);
// Search functions.
int           scsearch_fw            (spath_t, uint8_t *, int, int, int, int, int, int, index_t *, pathstack_t **);
int           seqsearch_fw           (spath_t, uint8_t *, int, int, int, int, int, index_t *, pathstack_t **);
int           seqsearch_bw           (spath_t, uint8_t *, int, int, int, int, int, index_t *, pathstack_t **);
int           seqdash_fw             (spath_t, uint8_t *, int, int, index_t *, pathstack_t **);
int           seqdash_bw             (spath_t, uint8_t *, int, int, index_t *, pathstack_t **);
// Stack tree functions
pathstack_t * pathstack_new          (int);
int           path_push              (spath_t, pathstack_t **);
pstree_t    * alloc_stack_tree       (int);
pstree_t    * alloc_stack_tree_rec   (int);
void          free_stack_tree        (pstree_t *);

// Misc.
void          mpos_set               (spath_t *, int);


#endif
