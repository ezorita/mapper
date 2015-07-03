#ifndef _BWM_DEFINITIONS_H
#define _BWM_DEFINITIONS_H

// Definitions
#define NUM_SYMBOLS 6
#define NUM_BASES 5
#define NUM_COMBS NUM_BASES*NUM_BASES*NUM_BASES
#define MAXSEQLEN 10000
#define MAXTAU    2
#define MAXSEQOUT 50

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))


// General sequence definitions.
#define EOS            -1

#endif
