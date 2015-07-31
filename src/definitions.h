#ifndef _BWM_DEFINITIONS_H
#define _BWM_DEFINITIONS_H

// Definitions
#define NUM_BASES 5
#define NUM_COMBS NUM_BASES*NUM_BASES*NUM_BASES
#define MAXSEQLEN 10000
#define MAXTAU    2
#define MAXSEQOUT 50
// Index params.
#define MMAP_FLAGS        (MAP_PRIVATE | MAP_POPULATE)


#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))


// General sequence definitions.
#define EOS            -1

#endif
