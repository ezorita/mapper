#include <stdlib.h>
#include <stdint.h>
#include "index.h"
#include "indexquery.h"
#include "algs.h"

typedef struct {
   fmdpos_t pos;
   uint32_t depth;
   uint32_t path;
   int8_t   row[2*MAXTAU+3];
} pebble_t;

typedef struct {
   size_t pos;
   size_t size;
   pebble_t pebble[];
} pstack_t;

typedef struct {
   uint8_t   * query;
   int         kmer;
   int         trail;
   int         tau;
   index_t   * index;
   pstack_t ** pebbles;
} arg_t;

typedef struct {
   int               kmer;
   int               tau;
   fmdpos_t          beg;
   uint64_t          end;
   uint8_t         * seq;
   uint64_t        * computed;
   uint8_t         * counts;
   int             * done;
   index_t         * index;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;

} annjob_t;

int         ppush         (pstack_t **, pebble_t);
pstack_t *  pstack_new    (size_t);
int         seq_lcp       (char *, char *, int);
int         find_next     (fmdpos_t, fmdpos_t *, uint8_t *, int, index_t *);
int         annotate      (int, int, uint8_t *, index_t *,int);
void *      annotate_mt   (void*);
int         poucet        (pebble_t, arg_t);
int         poucet_mismatch (pebble_t, arg_t);
int         dash          (pebble_t, const arg_t);
int         dash_mismatch (pebble_t, const arg_t);
