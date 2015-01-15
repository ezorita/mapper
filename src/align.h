#ifndef _ALIGN_H
#define _ALIGN_H

// SW-alignment
#define SCORE_MATCH    1
#define SCORE_DELETE   -1
#define SCORE_INSERT   -1

// This will be substituted in the future by a properly computed table.
// Now let's say that the max query len is 500 nt and our ID threshold is 0.6, so max error is 200nt.
#define ALIGN_WIDTH         100
#define MIN_ALIGNMENT_LEN   50

#define align_max(a,b) ((a) > (b) ? (a) : (b))
#define align_min(a,b) ((a) < (b) ? (a) : (b))

typedef struct align_t align_t;

struct align_t {
   long start;
   long max;
   long score;
   long ident;
};

align_t      sw_align         (char * read, int rdlen, char * ref, int rflen);

#endif
