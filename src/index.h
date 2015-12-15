#ifndef _INDEX_H
#define _INDEX_H

#define VERBOSE_DEBUG 0

static const char translate[256] = {[0 ... 255] = 4,
                           ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['n'] = 4, ['$'] = 5,
                           ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['N'] = 4};

static const char revcomp[256] = {[0 ... 255] = 'N',
                   ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A' };

static const char uppercase[256] = {[0 ... 255] = 'e',
                           ['a'] = 'a', ['c'] = 'b', ['g'] = 'c', ['t'] = 'd', ['n'] = 'e', ['$'] = '$',
                           ['A'] = 'a', ['C'] = 'b', ['G'] = 'c', ['T'] = 'd', ['N'] = 'e'};

typedef struct index_t    index_t;
typedef struct chr_t      chr_t;
typedef struct lcpdata_t  lcpdata_t;
typedef struct list64_t   list64_t;
typedef struct lcpval_t   lcpval_t;

struct lcpval_t {
   uint8_t lcp;
   int8_t  offset;
};

struct lcpdata_t {
   uint64_t   size;
   lcpval_t * lcp;
};

struct list64_t {
   uint64_t   size;
   int64_t  * val;
};

struct chr_t {
   int     nchr;
   long  * start;
   char ** name;
};

struct index_t {
   // Gen file.
   uint64_t    size;
   char      * genome;
   // OCC.
   uint64_t  * c;
   uint64_t    occ_mark_int;
   uint64_t  * occ;
   // Suffix Array.
   uint64_t    sa_bits;
   uint64_t  * sa;
   // Lookup table.
   //   uint64_t    lut_kmer;
   //   uint64_t  * lut;
   // LCP index.
   uint64_t    lcp_mark_int;
   uint64_t    lcp_min_depth;
   uint64_t  * lcp_sample_idx;
   uint64_t  * lcp_extend_idx;
   lcpdata_t * lcp_sample;
   list64_t  * lcp_extend;
   // Repeat annotation.
   uint8_t   * repeats;
   // Chromosome index.
   chr_t     * chr;
};

#endif
