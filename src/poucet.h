#include "definitions.h"
#include "algs.h"
#ifndef _BWM_POUCET_H
#define _BWM_POUCET_H

// Poucet algorithm definitions.
#define MAX_TRAIL           30
#define PATH_SCORE_BITS     16
#define PATH_SCORE_MASK     0x000000000000FFFF


// Poucet function headers.
int           poucet           (const long sp, const long ep, const char* prow, const int depth, char* path, arg_t * arg);
void          dash             (long sp, long ep, const int depth, const int align, const char* path, const arg_t* arg);

#endif
