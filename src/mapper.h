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
#include "indexquery.h"
#include "algs.h"

#ifndef _MAPPER_H
#define _MAPPER_

#define BUFFER_SIZE   100
#define CHRSTR_SIZE   50
#define SEQSTACK_SIZE 1024

typedef struct idxfiles_t idxfiles_t;

// Input types.

typedef enum {
   FASTA,
   FASTQ,
   RAW
} format_t;


struct idxfiles_t {
   // Original file pointers.
   void      * gen_file;
   void      * occ_file;
   void      * sa_file;
   void      * lcp_file;
   void      * lut_file;
   chr_t     * chr;
};

idxfiles_t    * index_open    (char * file);
index_t       * index_format  (idxfiles_t * files);
chr_t         * read_CHRindex (char * filename);
seqstack_t    * read_file     (FILE * inputf);

#endif
