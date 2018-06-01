#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include "errhandler.h"
#include "gstack.h"
#include "seqread.h"

#ifndef _IO_H
#define _IO_H

#define IO_DEFAULT_BUFSIZE 256*1024*1024 // 256MB

// Typedef structures.
typedef struct iostream_t iostream_t;

// Type interfaces.
struct iostream_t;

// Interface functions.
iostream_t     * io_stream_open        (char * path);
iostream_t     * io_stream_open_buf    (char * path, size_t bufsize);
int              io_stream_close       (iostream_t * stream);
void             io_stream_free        (iostream_t * stream);
gstack_t       * io_stream_read_seq    (iostream_t * stream);

#endif
