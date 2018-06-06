#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "errhandler.h"

#ifndef _INDEX_SYM_H
#define _INDEX_SYM_H

/* Magic number:
** bits 32-63 (mapper): 0x0fcb0fcb
** bits 16-31 (index_sym): 0x0001
** bits 00-15 (file type/version): 0x001
*/ 
#define SYM_FILE_MAGICNO      0x0fcb0fcb00010001

#define SYM_TABLE_SIZE        256
#define SYM_MAX_ALPHABET_SIZE 255

// Typedef structures
typedef struct sym_t sym_t;

// Type interface
struct sym_t;

// Interface functions.
sym_t    * sym_new            (char ** alphabet, char ** complement, uint8_t sym_default);
sym_t    * sym_new_dna        (void);
void       sym_free           (sym_t * sym);
int        sym_set_complement (char ** complement, sym_t  * sym);
char       sym_character      (uint8_t s, sym_t * sym);
int32_t    sym_complement     (uint8_t s, sym_t * sym);
int32_t    sym_index          (char c, sym_t * sym);
uint8_t  * sym_str_index      (char * str, sym_t * sym);
int        sym_is_canonical   (char c, sym_t * sym);

// Helper functions.
int        sym_count          (sym_t * sym);

// I/O functions.
int        sym_file_write     (char * filename, sym_t * sym);
sym_t    * sym_file_read      (char * filename);



#endif
