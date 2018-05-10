#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "index_sym.h"

#ifndef _INDEX_TXT_H
#define _INDEX_TXT_H

/* Magic number:
** bits 32-63 (mapper): 0x0fcb0fcb
** bits 16-31 (index_txt): 0x0002
** bits 00-15 (file type/version): 0x001
*/ 
#define TXT_FILE_MAGICNO      0x0fcb0fcb00020001

#define TXT_BUFFER_SIZE       1024
#define TXT_SEQ_SIZE          32

// Typedef structures.
typedef struct txt_t txt_t;

// Type interfaces.
struct txt_t;


// Interface functions.
txt_t    * txt_new             (sym_t * sym);
void       txt_free            (txt_t * txt);
int32_t    txt_sym             (int64_t pos, txt_t * txt);
uint8_t  * txt_sym_range       (int64_t beg, int64_t len, txt_t * txt);
int        txt_append          (char * text, txt_t * txt);
int        txt_append_wildcard (txt_t * txt);
int        txt_commit_seq      (char * seqname, txt_t * txt);
int        txt_commit_rc       (txt_t * txt);

// Helper functions.
int64_t    txt_length          (txt_t * txt);
sym_t    * txt_get_symbols     (txt_t * txt);
int64_t    txt_wildcard_count  (txt_t * txt);

// Sequence helper functions.
int64_t    txt_seq_count       (txt_t * txt);
int64_t    txt_seq_start       (int32_t seq, txt_t * txt);
int64_t    txt_seq_length      (int32_t seq, txt_t * txt);
char     * txt_seq_name        (int32_t seq, txt_t * txt);
char     * txt_pos_to_str      (int64_t pos, txt_t * txt);
int64_t    txt_str_to_pos      (char *  str, txt_t * txt);

// I/O functions.
int        txt_file_write     (char * filename, txt_t * txt);
txt_t    * txt_file_read      (char * filename, sym_t * sym);


#endif
