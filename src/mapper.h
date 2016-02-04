#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>
#include <pthread.h>
#include "algs.h"
#include "seed.h"
#include "bwtquery.h"
#include "indexbuild.h"
#include "filter.h"
#include "format.h"
#include "annotate.h"
#include "blocksearch.h"
#include "interface.h"

#ifndef _MAPPER_H
#define _MAPPER_H

#define SEQSTACK_SIZE 1024

typedef struct idxfiles_t  idxfiles_t;
typedef struct mtjob_t     mtjob_t;
typedef struct mtcontrol_t mtcontrol_t;
typedef struct param_map_t param_map_t;

// Input types.

typedef enum {
   FASTA,
   FASTQ,
   RAW
} format_t;

struct mtjob_t {
   size_t            count;
   seq_t           * seq;
   index_t         * index;
   param_map_t     * opt;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   mtcontrol_t     * control;
};

struct param_map_t {
   seedopt_t   seed;
   filteropt_t filter;
   alignopt_t  align;
   formatopt_t format;
   int threads;
};

struct mtcontrol_t {
   uint64_t mapped;
   int      active;
};


// Thread management.
int             mt_scheduler (seqstack_t *, param_map_t *, index_t *, int, size_t);
void          * mt_worker    (void * args);

// File/Index management.
idxfiles_t    * index_open      (char * file);
index_t       * index_format    (idxfiles_t * files);
chr_t         * read_CHRindex   (char * filename);
seqstack_t    * read_file       (FILE * inputf);

#endif
