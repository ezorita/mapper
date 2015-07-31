#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <string.h>
#include "definitions.h"
#include "divsufsort.h"
#include "indexquery.h"

#define STACK_LCP_INITIAL_SIZE 1024
#define BUFFER_SIZE   100
#define GENOME_SIZE   100000


static const char translate[256] = {[0 ... 255] = 4,
                           ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['n'] = 4,
                           ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['N'] = 4, ['$'] = 5 };

static const char revcomp[256] = {[0 ... 255] = 'N',
                   ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A' };

static const char uppercase[256] = {[0 ... 255] = 'e',
                           ['a'] = 'a', ['c'] = 'b', ['g'] = 'c', ['t'] = 'd', ['n'] = 'e',
                           ['A'] = 'a', ['C'] = 'b', ['G'] = 'c', ['T'] = 'd', ['N'] = 'e', ['$'] = '$'};

typedef struct lcpstack_t  lcpstack_t;
typedef struct lcpsample_t lcpsample_t;
typedef struct lcpcorner_t lcpcorner_t;
typedef struct cstack_t    cstack_t;

struct lcpcorner_t {
   int     lcp;
   int64_t pos;
};

struct cstack_t {
   uint64_t    pos;
   uint64_t    size;
   lcpcorner_t c[];
};

struct lcpstack_t {
   uint64_t pos;
   uint64_t size;
   int8_t   lcp[];
};

struct lcpsample_t {
   uint8_t lcp;
   uint64_t ep;
};

int          write_index      (char * filename);
int64_t    * compute_sa       (char * genome, uint64_t gsize);
uint64_t   * compute_occ      (char * genome, uint64_t * sa, uint64_t gsize, uint64_t * occ_size, uint64_t * wbwt);
uint64_t   * compute_c        (uint64_t * occ);
uint64_t   * compute_lut      (uint64_t * c, uint64_t * occ, int depth);
void         recursive_lut    (uint64_t * c, uint64_t * occ, uint64_t * lut, uint64_t p, int d, int maxd, int * path);
int          lcp_stack_push8  (lcpstack_t ** lcp_stack, uint8_t lcp, int8_t offset);
int          lcp_stack_push32 (lcpstack_t ** lcp_stack, uint8_t lcp, int32_t offset);
void         lcp_index_sample (uint64_t * lcp_index, uint64_t pos, int size32);
lcpstack_t * compute_lcp      (uint64_t * index_size, uint64_t ** lcp_index, int min_depth, int sar_bits, uint64_t * sar, uint64_t idxsize, char * genome);
int64_t      recursive_lcp    (fmdpos_t pos, int depth, index_t * index, uint64_t * lcp_index, lcpstack_t ** lcp, int mindepth);
char *       compact_genome   (char * filename, uint64_t * genomesize);
uint64_t     compact_array    (uint64_t * array, uint64_t len, int bits);
int64_t      naive_lcp        (int min_depth, uint64_t idxsize, uint64_t * lcpindex, lcpstack_t ** lcp, uint64_t * sar, int sar_bits, char * genome);
lcpcorner_t  corner_pop       (cstack_t * stack);
int          corner_push      (cstack_t ** stack, lcpcorner_t c);
cstack_t   * cstack_new       (int size);
int          seq_lcp          (char * seq_a, char * seq_b);
