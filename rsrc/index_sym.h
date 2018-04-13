#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#ifndef _INDEX_SYM_H
#define _INDEX_SYM_H

#define SYM_TABLE_SIZE        256
#define SYM_MAX_ALPHABET_SIZE 255

// Typedef structures
typedef struct sym_t sym_t;

// Type interface
struct sym_t;

// Interface functions.
sym_t    * sym_new            (char ** alphabet, char ** complement, uint8_t sym_default);
void       sym_free           (sym_t * sym);
int        sym_set_complement (char ** complement, sym_t  * sym);
char       sym_character      (uint8_t s, sym_t * sym);
uint8_t    sym_index          (char c, sym_t * sym);
int        sym_is_canonical   (char c, sym_t * sym);

#endif
