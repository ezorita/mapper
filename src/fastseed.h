#include "definitions.h"
#include "algs.h"
#include "align.h"

#ifndef _BWM_FASTSEED_H
#define _BWM_FASTSEED_H

#define ALIGNMENT_EXTRA_L  100
#define ALIGNMENT_EXTRA_R  100

int          fastseed         (seqstack_t* seqs, index_t index, chr_t* chr, int k);
int          fastseed_filter  (seqstack_t* seqs, index_t index, chr_t* chr, int k);
int          query_index      (char* query, long gsize, long* c, long* ptr, list_t* occs);
void         translate_query  (char* query, int* qval, int qlen);
#endif
