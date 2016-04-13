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
#include "definitions.h"
#include "annotate.h"

#ifndef _INDEX_H
#define _INDEX_H

#define CHRSTR_SIZE   50
#define BUFFER_SIZE   100

char       * add_suffix            (const char *, const char *);
annlist_t  * ann_index_read        (char *);
shtlist_t  * sht_index_read        (char *);
int          ann_index_write       (annlist_t *, char *);
int          sht_index_write       (shtlist_t *, char *);
int          ann_index_load        (annlist_t *, char *);
int          sht_index_load        (shtlist_t *, char *);
int          ann_find_slot         (annlist_t *);
int          sht_find_slot         (shtlist_t *);
void         ann_print_index       (annlist_t *);
void         sht_print_index       (shtlist_t *);
int          index_add_annotation  (int, int, int, int, int, int, index_t *, char *);
index_t    * index_load_base       (char *);
char       * index_load_gen        (char *);
bwt_t      * index_load_bwt        (char *);
sar_t      * index_load_sar        (char *);
chr_t      * index_load_chr        (char *);
int          ann_read              (ann_t, uint64_t, int *);
int          ann_find              (int, annlist_t *, ann_t *);
#endif
