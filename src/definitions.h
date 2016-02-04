#include <sys/mman.h>

#ifndef _BWM_DEFINITIONS_H
#define _BWM_DEFINITIONS_H

#define VERBOSE_DEBUG 0

static const char bases[5] = "ACGTN";

static const char translate[256] = {[0 ... 255] = 4,
                           ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3, ['n'] = 4, ['$'] = 5,
                           ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3, ['N'] = 4};

static const char translate_rc[256] = {[0 ... 255] = 4,
                           ['a'] = 3, ['c'] = 2, ['g'] = 1, ['t'] = 0, ['n'] = 4, ['$'] = 5,
                           ['A'] = 3, ['C'] = 2, ['G'] = 1, ['T'] = 0, ['N'] = 4};

static const char revcomp[256] = {[0 ... 255] = 'N',
                   ['a'] = 't', ['c'] = 'g', ['g'] = 'c', ['t'] = 'a', ['u'] = 'a',
                   ['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['U'] = 'A' };

static const char uppercase[256] = {[0 ... 255] = 'e',
                           ['a'] = 'a', ['c'] = 'b', ['g'] = 'c', ['t'] = 'd', ['n'] = 'e', ['$'] = '$',
                           ['A'] = 'a', ['C'] = 'b', ['G'] = 'c', ['T'] = 'd', ['N'] = 'e'};


// Definitions
#define NUM_BASES 5
#define NUM_COMBS NUM_BASES*NUM_BASES*NUM_BASES
#define MAXSEQLEN 10000
#define MAXTAU    2
#define MAXSEQOUT 50
// Index params.
#define MMAP_FLAGS        (MAP_PRIVATE | MAP_POPULATE)

// General sequence definitions.
#define EOS            -1

typedef struct index_t    index_t;
typedef struct bwt_t      bwt_t;
typedef struct sar_t      sar_t;
typedef struct chr_t      chr_t;

struct chr_t {
   int     nchr;
   long  * start;
   char ** name;
};

struct bwt_t {
   uint64_t   occ_mark_int;
   uint64_t * c;
   uint64_t * occ;
};

struct sar_t {
   uint64_t   bits;
   uint64_t * sa;
};

struct index_t {
   uint64_t    size;
   // Gen file.
   char      * genome;
   // OCC.
   bwt_t     * bwt;
   // Suffix Array.
   sar_t     * sar;
   // Chromosome index.
   chr_t     * chr;
   // Annotation.
   void      * seeds;
   uint8_t   * repeats;
};


#endif
