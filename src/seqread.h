#include <stdlib.h>
#include <stdint.h>
#include "errhandler.h"
#include "index_sym.h"
#include "gstack.h"

#ifndef _SEQREAD_H
#define _SEQREAD_H

// Typedef structures.
typedef struct seqread_t seqread_t;

// Type interfaces.
struct seqread_t;

// Interface functions.
seqread_t     * seqread_new     (char * tag, char * seq, char * qscore, sym_t * sym);
void            seqread_free    (void * ptr);
gstack_t      * seqread_stack   (size_t max_elm);
seqread_t     * seqread_pop     (gstack_t * stack);
seqread_t     * seqread_get     (size_t index, gstack_t * stack);

// Helper functions.
int64_t         seqread_len     (seqread_t * read);
char          * seqread_seq     (seqread_t * read);
char          * seqread_tag     (seqread_t * read);
char          * seqread_qscore  (seqread_t * read);
uint8_t       * seqread_sym     (seqread_t * read);

#endif
