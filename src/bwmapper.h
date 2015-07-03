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
#include "definitions.h"
#include "dc3.h"
#include "algs.h"
#include "hitmap.h"

#ifndef _BWMAPPER_H
#define _BWMAPPER_H

// Index params.
#define OCC_MARK_INTERVAL 14
#define OCC_WORD_SIZE     64
#define OCC_MARK_BITS     (OCC_MARK_INTERVAL * OCC_WORD_SIZE)
#define OCC_WORD_MASK     0xFFFFFFFFFFFFFFFF
#define MMAP_FLAGS        (MAP_PRIVATE | MAP_POPULATE)
// Initial buffer/stack sizes
#define BUFFER_SIZE   100
#define GENOME_SIZE   100000
#define CHRSTR_SIZE   100
#define SEQSTACK_SIZE 1024
#define DSTACK_SIZE   16

// Input types.

typedef enum {
   FASTA,
   FASTQ,
   RAW
} format_t;

// Index functions.
int           format_FMindex   (long* index, index_t* fmindex);
chr_t       * read_CHRindex    (char* filename);
ssize_t       write_index      (char* filename);
char        * compact_genome   (char* filename, long* genomesize);
void          bwt_index        (char* text, long tlen, long** pos, uint64_t** occ, uint64_t* occ_size);
long        * compute_c        (char* genome, long gsize);

// Misc functions.
seqstack_t  * read_file        (FILE * inputf, const int reverse);

#endif
