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
#define max(a,b) (a > b ? a : b)

char translate[256] = {[0 ... 255] = 3, ['@'] = 0,
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
   long   pos;
   long   size;
   void * job[];
} jstack_t;

void        bwt_index       (char* text, long tlen, long** pos, vstack_t** occ);
void        prefix_sort     (char * text, long tlen, long ** cp, vstack_t ** sorted);
void        push            (vstack_t ** stackp, long value);
vstack_t  * new_stack       (long size);
long      * read_index      (char* filename, long* gsize, long** Cp, long** posp, list_t** occp);
void        write_index     (char* filename);
char      * compact_genome  (char* filename, long* genomesize);
vstack_t  * query_index     (char* query, long gsize, long* c, long* pos, list_t* occs);
int       * translate_query (char* query);
long        bisect_search   (long start, long end, long* set, long value);
jstack_t  * new_jstack      (long size);
void        push_job        (jstack_t** stackp, void* job);
void      * pop_job         (jstack_t* stack);
long      * compute_c       (char* genome, long gsize);
