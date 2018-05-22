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
