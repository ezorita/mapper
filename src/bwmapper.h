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
index_t     * load_index       (char* file);
chr_t       * read_CHRindex    (char* filename);
ssize_t       write_index      (char* filename);
char        * compact_genome   (char* filename, uint64_t* genomesize);
void          bwt_index        (char* text, uint64_t tlen, uint64_t** pos, uint64_t** occ, uint64_t* occ_size, uint64_t ** C);
uint64_t    * compute_c        (char* genome, long gsize);

// Misc functions.
seqstack_t  * read_file        (FILE * inputf, const int reverse);

#endif
