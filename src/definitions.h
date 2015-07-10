#ifndef _BWM_DEFINITIONS_H
#define _BWM_DEFINITIONS_H

// Definitions
#define NUM_BASES 5
#define NUM_COMBS NUM_BASES*NUM_BASES*NUM_BASES
#define MAXSEQLEN 10000
#define MAXTAU    2
#define MAXSEQOUT 50
// Index params.
#define OCC_MARK_INTERVAL 14
#define OCC_WORD_SIZE     64
#define OCC_MARK_BITS     (OCC_MARK_INTERVAL * OCC_WORD_SIZE)
#define MMAP_FLAGS        (MAP_PRIVATE | MAP_POPULATE)


#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))


// General sequence definitions.
#define EOS            -1

#endif
