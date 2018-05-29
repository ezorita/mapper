#include <stdlib.h>
#include <stdint.h>
#include "errhandler.h"
#include "index_sym.h"

#ifndef _SEQREAD_H
#define _SEQREAD_H

// Typedef structures.
typedef struct seqread_t seqread_t;

// Type interfaces.
struct seqread_t;

// Interface functions.
seqread_t     * seqread_new     (char * tag, char * seq, char * qscore, sym_t * sym);
void            seqread_free    (seqread_t * read);

// Helper functions.
int64_t         seqread_len     (seqread_t * read);
char          * seqread_seq     (seqread_t * read);
char          * seqread_tag     (seqread_t * read);
char          * seqread_qscore  (seqread_t * read);
uint8_t       * seqread_sym     (seqread_t * read);

#endif
