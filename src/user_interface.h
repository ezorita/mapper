#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "errhandler.h"
#include "version.h"
#include "mapper.h"
#include "index.h"

#ifndef _INTERFACE_H
#define _INTERFACE_H

typedef struct opt_map_t   opt_map_t;
typedef struct opt_add_t   opt_add_t;

char * USAGE_MAP;
char * USAGE_INDEX;
char * USAGE_BUILD;
char * USAGE_ADD;
char * USAGE_VIEW;
char * ERROR_MAX_ARG;
char * ERROR_INSUF_ARG;
char * ERROR_COMMAND;
char * ERROR_INCORRECT_OPT;

struct opt_map_t {
   int print_first;
   int eval_thr;
   int mapq_thr;
   int threads;
};

struct opt_add_t {
   int k;
   int d;
   int threads;
};

int    ui_parse  (int argc, char * argv[]);

#endif
