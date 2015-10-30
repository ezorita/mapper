#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "index.h"
#include "algs.h"

#ifndef _FORMAT_H
#define _FORMAT_H

typedef struct formatopt_t formatopt_t;

struct formatopt_t {
   int mapq_thr;
   int print_first;
   double eval_thr;
};

int print_and_free (seq_t seq, matchlist_t *, index_t * index, formatopt_t opt);

#endif
