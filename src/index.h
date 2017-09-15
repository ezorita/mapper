#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <glob.h>
#include "definitions.h"
#include "annotate.h"


#ifndef _INDEX_H
#define _INDEX_H

#define CHRSTR_SIZE   50
#define BUFFER_SIZE   100

char       * add_suffix            (const char *, const char *);
annlist_t  * ann_index_read        (char *);
int          ann_index_write       (annlist_t *, char *);
int          ann_new_filename      (int, int, char*);
void         ann_print_index       (annlist_t *);
int          index_add_annotation  (int, int, int, index_t *, char *);
index_t    * index_load_base       (char *);
char       * index_load_gen        (char *);
bwt_t      * index_load_bwt        (char *);
sar_t      * index_load_sar        (char *);
chr_t      * index_load_chr        (char *);
int          ann_read              (ann_t, uint64_t, int *);
ann_t      * ann_find              (int, annlist_t *);
int          compar_ann_k_asc      (const void *, const void *);
#endif
