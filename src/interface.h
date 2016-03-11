#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

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
   int sd;
   int threads;
   int repeat_thr;
   int mode;
};

void say_map_usage(void);
void say_add_usage(void);
void say_index_usage(void);
void say_build_usage(void);
int  parse_opt_build(int, char **, opt_add_t *, char **);
int  parse_opt_add(int, char **, opt_add_t *, char **);
int  parse_opt_map(int, char **, opt_map_t *, char **, char **);

#endif
