#include <stdint.h>
#include <stdlib.h>
#include "index.h"
#include "indexquery.h"

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

void            blocksearch     (uint8_t *, int, int, index_t *, pathstack_t **);
void            blocksearch_rec (uint8_t *, int, int, index_t *, pathstack_t **);
int             seqsearch_fw    (spath_t, uint8_t *, int, int, int, int, int, index_t *, pathstack_t **);
int             seqsearch_bw    (spath_t, uint8_t *, int, int, int, int, int, index_t *, pathstack_t **);
int             seqdash_fw      (spath_t *, uint8_t *, int, int, index_t *);
int             seqdash_bw      (spath_t *, uint8_t *, int, int, index_t *);
int             path_push       (spath_t, pathstack_t **);
pathstack_t   * pathstack_new   (int);
