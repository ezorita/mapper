#include <stdint.h>
#include <stdlib.h>
#include "bwtquery.h"
#include "definitions.h"

#ifndef _BLOCKSEARCH_H
#define _BLOCKSEARCH_H
#define PATHSTACK_DEF_SIZE 64

typedef struct {
   fmdpos_t pos;
   int      score;
} spath_t;

typedef struct {
   size_t pos;
   size_t size;
   spath_t path[];
} pathstack_t;

typedef struct pstree_t pstree_t;

struct pstree_t {
   pathstack_t * stack;
   pstree_t    * next_l;
   pstree_t    * next_r;
};

void          blocksearch            (uint8_t *, int, int, index_t *, pathstack_t **);
void          blocksearch_rec        (uint8_t *, int, int, index_t *, pathstack_t **);
void          blocksearch_trail      (uint8_t *, int, int, int, index_t *, pstree_t *);
void          blocksearch_trail_rec  (uint8_t *, int, int, int, int, index_t *, pstree_t *);
int           seqsearch_fw           (spath_t, uint8_t *, int, int, int, int, int, index_t *, pathstack_t **);
int           seqsearch_bw           (spath_t, uint8_t *, int, int, int, int, int, index_t *, pathstack_t **);
int           seqdash_fw             (spath_t *, uint8_t *, int, int, index_t *);
int           seqdash_bw             (spath_t *, uint8_t *, int, int, index_t *);
int           path_push              (spath_t, pathstack_t **);
pathstack_t * pathstack_new          (int);
pstree_t    * alloc_stack_tree       (int);
pstree_t    * alloc_stack_tree_rec   (int);
void          free_stack_tree        (pstree_t *);
int           compar_path_score      (const void *, const void *);
#endif
