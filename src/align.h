#ifndef _ALIGN_H
#define _ALIGN_H

// SW-alignment
#define SCORE_MATCH    1
#define SCORE_DELETE   -1
#define SCORE_INSERT   -1

#define align_max(a,b) ((a) > (b) ? (a) : (b))
#define align_min(a,b) ((a) < (b) ? (a) : (b))

typedef struct align_t align_t;

struct align_t {
   long start;
   long max;
   long score;
};

align_t      sw_align         (char * read, int rdlen, char * ref, int rflen);

#endif
