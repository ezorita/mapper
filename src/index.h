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

typedef struct ann_t       ann_t;
typedef struct sht_t       sht_t;
typedef struct annlist_t   annlist_t;
typedef struct shtlist_t   shtlist_t;


#define CHRSTR_SIZE   50
#define BUFFER_SIZE   100

struct ann_t {
   uint8_t    id;
   int        k;
   int        d;
   uint64_t   size;
   uint64_t   unique;
   uint8_t  * data;
};


struct annlist_t {
   uint8_t  count;
   ann_t    ann[];
};

struct sht_t {
   uint8_t    id;
   uint8_t    bits;
   int        k;
   int        d;
   int        repeat_thr;
   uint64_t   set_count;
   uint64_t   collision;
   htable_t * htable;
};

struct shtlist_t {
   uint8_t  count;
   sht_t    sht[];
};

char       * add_suffix            (const char *, const char *);
annlist_t  * ann_index_read        (char *);
shtlist_t  * sht_index_read        (char *);
int          ann_index_write       (annlist_t *, char *);
int          sht_index_write       (shtlist_t *, char *);
int          ann_index_load        (annlist_t *, char *);
int          sht_index_load        (shtlist_t *, char *);
int          ann_find_slot         (annlist_t *);
int          sht_find_slot         (shtlist_t *);
int          index_add_annotation  (int, int, int, int, int, index_t *, char *);
index_t    * index_load_base       (char *);
char       * index_load_gen        (char *);
bwt_t      * index_load_bwt        (char *);
sar_t      * index_load_sar        (char *);
chr_t      * index_load_chr        (char *);

#endif
