#define _GNU_SOURCE
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <errno.h>
#include <string.h>
#include "definitions.h"
#include "divsufsort.h"
#include "index.h"
#include "bwtquery.h"

#ifndef _INDEXBUILD_H
#define _INDEXBUILD_H

#define STACK_LCP_INITIAL_SIZE 1024
#define NAIVE_LCP_MAX_SAMPLES 1000000
#define BUFFER_SIZE   100
#define GENOME_SIZE   100000

int          write_index      (char * filename, int def_kmer, int def_tau, int threads);
int64_t    * compute_sa       (char * genome, uint64_t gsize);
uint64_t   * compute_occ      (char * genome, uint64_t * sa, uint64_t gsize, uint64_t * occ_size, uint64_t * wildcard_cnt);
uint64_t   * compute_c        (uint64_t * occ, uint64_t wildcard_cnt);
uint64_t     compact_array    (uint64_t * array, uint64_t len, int bits);
char       * compact_genome   (char * filename, uint64_t * genomesize);

#endif
