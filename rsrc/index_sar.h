#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <glob.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>

#include "index_txt.h"
#include "divsufsort.h"

#ifndef _INDEX_SAR_H
#define _INDEX_SAR_H

/* Magic number:
** bits 32-63 (mapper): 0x0fcb0fcb
** bits 16-31 (index_txt): 0x0004
** bits 00-15 (file type/version): 0x001
*/ 
#define SAR_FILE_MAGICNO      0x0fcb0fcb00040001

// Typedef structures.
typedef struct sar_t sar_t;

// Type interfaces.
struct sar_t;

// Interface functions.
sar_t     *  sar_build        (txt_t * txt);
void         sar_free         (sar_t * sar);
int64_t      sar_get          (int64_t pos, sar_t * sar);
int          sar_get_range    (int64_t beg, int64_t len, int64_t * vec, sar_t * sar);

// I/o functions.
int          sar_file_write   (char * filename, sar_t * sar);
sar_t     *  sar_file_read    (char * filename);
#endif
