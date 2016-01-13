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
#include <pthread.h>
#include "algs.h"
#include "seed.h"
#include "indexquery.h"
#include "filter.h"
#include "format.h"
#include "annotate.h"

#ifndef _MAPPER_H
#define _MAPPER_H

#define BUFFER_SIZE   100
#define CHRSTR_SIZE   50
#define SEQSTACK_SIZE 1024

typedef struct idxfiles_t  idxfiles_t;
typedef struct mapopt_t    mapopt_t;
typedef struct mtjob_t     mtjob_t;
typedef struct mtcontrol_t mtcontrol_t;

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
   mapopt_t        * opt;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   mtcontrol_t     * control;
};

struct mtcontrol_t {
   uint64_t mapped;
   int      active;
};

struct idxfiles_t {
   // Original file pointers.
   void      * gen_file;
   size_t      gen_len;
   void      * occ_file;
   size_t      occ_len;
   void      * sa_file;
   size_t      sa_len;
   void      * lcp_file;
   size_t      lcp_len;
   void      * ann_file;
   size_t      ann_len;
//   void      * lut_file;
   chr_t     * chr;
};


struct mapopt_t {
   seedopt_t   seed;
   filteropt_t filter;
   alignopt_t  align;
   formatopt_t format;
   int threads;
};

// Thread management.
int             mt_scheduler (seqstack_t *, mapopt_t *, index_t *, int, size_t);
void          * mt_worker    (void * args);

// File/Index management.
idxfiles_t    * index_open      (char * file);
index_t       * index_format    (idxfiles_t * files);
chr_t         * read_CHRindex   (char * filename);
seqstack_t    * read_file       (FILE * inputf);

#endif
