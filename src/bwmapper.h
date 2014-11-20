#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <execinfo.h>
#include <signal.h>
#include "dc3.h"

// Definitions
#define NUM_BASES 6
#define NUM_COMBS NUM_BASES*NUM_BASES*NUM_BASES
#define BUFFER_SIZE 100
#define GENOME_SIZE 100000
#define CHRSTR_SIZE 100
#define max(a,b) (a > b ? a : b)

char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

typedef struct {
   long pos;
   long size;
   long offset;
   long val[];
} vstack_t;

typedef struct {
   long   max;
   long * val;
} list_t;

typedef struct {
   int     nchr;
   long *  start;
   char ** name;
} chr_t;

// Stack functions.
vstack_t  * new_stack       (long size);
void        push            (vstack_t ** stackp, long value);

// Query functions.
long      * read_FMindex    (char* filename, long* gsize, long** Cp, long** posp, list_t** occp);
chr_t       read_CHRindex   (char* filename);
vstack_t  * query_index     (char* query, long gsize, long* c, long* pos, list_t* occs);
int       * translate_query (char* query);
long        bisect_search   (long start, long end, long* set, long value);

// Index functions.
ssize_t     write_index     (char* filename);
char      * compact_genome  (char* filename, long* genomesize);
void        bwt_index       (char* text, long tlen, long** pos, vstack_t** occ);
long      * compute_c       (char* genome, long gsize);




