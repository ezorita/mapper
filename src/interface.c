#include "interface.h"

char *USAGE_MAP =
   "\n"
   "usage:\n"
   "  mapper [options] index-file input-file\n"
   "\n"
   "  options:\n"
   "    -t --threads: number of concurrent threads (default 1)\n"
   "    -a --all:     print all the hits found in each interval.\n"
   "    -q --quality: do not print reads with mapping quality below [quality].\n"
   "    -e --eval:    do not print reads with significance below [e-value].\n"
   "\n"
   " to build/edit/extend an index-file type:\n"
   "  \"mapper index -h\"\n";

char *USAGE_ADD =
   "\n"
   "usage:\n"
   " mapper index add [options] index-file\n"
   "\n"
   "  options:\n"
   "    -k --kmer:       sequence length. [required]\n"
   "    -d --distance:   sequence mismatches for annotation. [required]\n"
   "    -s --seed_dist:  set different distance for seed table.\n"
   "    -r --repeat-thr: repeat threshold. (default 100)\n"
   "    -t --threads:    number of threads. (default 1)\n"
   "    --seed-save:     store the seed hash table.\n"
   "    --seed-only:     store only the seed hash table.\n";

char *USAGE_INDEX =
   "\n"
   "usage:\n"
   " mapper index [command]\n"
   "\n"
   "  commands:\n"
   "    build:  build a new index from scratch.\n"
   "    add:    add annotation data to an existing index.\n"
   "    delete: delete annotation data from an existing index.\n"
   "    view:   view index state and annotations.\n"
   "\n"
   "  for more help type: index command\n";

char *USAGE_BUILD =
   "\n"
   "usage:\n"
   " mapper index build [options] genome.fasta\n"
   "\n"
   "  creates a new index for the specified genome file.\n"
   "  the index will be created along with a default\n"
   "  annotation and seed table.\n"
   "\n"
   "  options:\n"
   "    -k --kmer:       sequence length for default annotation. (default 25)\n"
   "    -d --distance:   sequence mismatches for default annotation. (default 1)\n"
   "    -r --repeat-thr: repeat threshold for default annotation. (default 20)\n"
   "    -t --threads:    number of threads. (default 1)\n"
   "\n"
   "  files that will be generated:\n"
   "    genome.fasta.gen:  genome sequence\n"
   "    genome.fasta.sar:  suffix array\n"
   "    genome.fasta.bwt:  bwt index\n"
   "    genome.fasta.ann:  annotation index\n"
   "    genome.fasta.sht:  seed table index\n"
   "    genome.fasta.chr:  chromosome index\n"
   "  annotation and seed table data:\n"
   "    genome.fasta.ann.#\n"
   "    genome.fasta.sht.#\n";

void say_map_usage(void) { fprintf(stderr, "%s\n", USAGE_MAP); }
void say_add_usage(void) { fprintf(stderr, "%s\n", USAGE_ADD); }
void say_index_usage(void) { fprintf(stderr, "%s\n", USAGE_INDEX); }
void say_build_usage(void) { fprintf(stderr, "%s\n", USAGE_BUILD); }

int
parse_opt_build
(
 int          argc,
 char       * argv[],
 opt_add_t  * opt,
 char      ** gfile
 )
{
   int arg_k = -1, arg_d = -1, arg_t = -1, arg_r = -1;
   int c;
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"threads",    required_argument,      0, 't'},
         {"kmer",       required_argument,      0, 'k'},
         {"distance",   required_argument,      0, 'd'},
         {"repeat-thr", required_argument,      0, 'r'},
         {"help",       no_argument,            0, 'h'},
         {0, 0, 0, 0}
      };

      c = getopt_long(argc, argv, "t:k:d:r:h", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
      case 't':
         if (arg_t < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "[error] threads option (-t) must be a positive number.\n");
               return -1;
            }
            arg_t = v;
         } else {
            fprintf(stderr,"[error] threads option (-t) set more than once.\n");
         }
         break;

      case 'k':
         if (arg_k < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "[error] kmer option (-k) must be a positive number.\n");
               return -1;
            }
            arg_k = v;
         } else {
            fprintf(stderr,"[error] kmer option (-k) set more than once.\n");
         }
         break;

      case 'd':
         if (arg_d < 0) {
            int v = atoi(optarg);
            if (v < 0) {
               fprintf(stderr, "[error] distance option (-d) must be a non-negative number.\n");
               return -1;
            }
            arg_d = v;
         } else {
            fprintf(stderr,"[error] distance option (-d) set more than once.\n");
         }
         break;

      case 'r':
         if (arg_r < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "[error] repeat_thr option (-r) must be a positive number.\n");
               return -1;
            }
            arg_r = v;
         } else {
            fprintf(stderr,"[error] repeat_thr option (-r) set more than once.\n");
         }
         break;

      case 'h':
         say_build_usage();
         return 1;

      }
   }

   if (optind + 2 != argc-1) {
      fprintf(stderr, "[error] not enough parameters.\n");
      say_add_usage();
      return -1;
   }

   *gfile = argv[argc-1];

   opt->k          = (arg_k < 0 ? 25 : arg_k);
   opt->d          = (arg_d < 0 ? 1 : arg_d);
   opt->threads    = (arg_t < 0 ? 1 : arg_t);
   opt->repeat_thr = (arg_r < 0 ? 100 : arg_r);

   return 0;
}


int
parse_opt_add
(
 int          argc, 
 char       * argv[],
 opt_add_t  * opt,
 char      ** ifile
)
{
   int arg_k = -1, arg_d = -1, arg_t = -1, arg_r = -1, arg_s = -1;
   static int arg_m = 1;
   int c;
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"seed-save",  no_argument,       &arg_m,  3 },
         {"seed-only",  no_argument,       &arg_m,  2 },
         {"threads",    required_argument,      0, 't'},
         {"kmer",       required_argument,      0, 'k'},
         {"distance",   required_argument,      0, 'd'},
         {"seed_dist",  required_argument,      0, 's'},
         {"repeat-thr", required_argument,      0, 'r'},
         {"help",       no_argument,            0, 'h'},
         {0, 0, 0, 0}
      };

      c = getopt_long(argc, argv, "t:k:d:r:s:h", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
      case 't':
         if (arg_t < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "[error] threads option (-t) must be a positive number.\n");
               return -1;
            }
            arg_t = v;
         } else {
            fprintf(stderr,"[error] threads option (-t) set more than once.\n");
         }
         break;

      case 'k':
         if (arg_k < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "[error] kmer option (-k) must be a positive number.\n");
               return -1;
            }
            arg_k = v;
         } else {
            fprintf(stderr,"[error] kmer option (-k) set more than once.\n");
         }
         break;

      case 'd':
         if (arg_d < 0) {
            int v = atoi(optarg);
            if (v < 0) {
               fprintf(stderr, "[error] distance option (-d) must be a non-negative number.\n");
               return -1;
            }
            arg_d = v;
         } else {
            fprintf(stderr,"[error] distance option (-d) set more than once.\n");
         }
         break;

      case 's':
         if (arg_s < 0) {
            int v = atoi(optarg);
            if (v < 0) {
               fprintf(stderr, "[error] seed_dist option (-s) must be a non-negative number.\n");
               return -1;
            }
            arg_s = v;
         } else {
            fprintf(stderr,"[error] seed_dist option (-s) set more than once.\n");
         }
         break;

      case 'r':
         if (arg_r < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "[error] repeat_thr option (-r) must be a positive number.\n");
               return -1;
            }
            arg_r = v;
         } else {
            fprintf(stderr,"[error] repeat_thr option (-r) set more than once.\n");
         }
         break;

      case 'h':
         say_add_usage();
         return 1;

      }
   }

   if (optind + 2 != argc-1) {
      fprintf(stderr, "[error] not enough parameters.\n");
      say_add_usage();
      return -1;
   }

   *ifile = argv[argc-1];

   if (arg_k < 0 || arg_d < 0) {
      fprintf(stderr, "[error] kmer (-k) and distance (-d) are required parameters.\n");
      return -1;
   }

   opt->k          = arg_k;
   opt->d          = arg_d;
   opt->sd         = (arg_s < 0 ? arg_d :arg_s);
   opt->threads    = (arg_t < 0 ? 1 : arg_t);
   opt->repeat_thr = (arg_r < 0 ? 20 : arg_r);
   opt->mode       = arg_m;

   return 0;
}

int
parse_opt_map
(
 int         argc,
 char      * argv[],
 opt_map_t * opt,
 char     ** qfile,
 char     ** ifile
)
{
   int arg_t = 0, arg_a = 0, arg_e = 0, arg_q = 0;
   int c;
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"all",      no_argument,       0, 'a'},
         {"threads",  required_argument, 0, 't'},
         {"eval",     required_argument, 0, 'e'},
         {"quality",  required_argument, 0, 'q'},
         {0, 0, 0, 0}
      };

      c = getopt_long(argc, argv, "at:e:q:",
            long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'a':
         if (arg_a) {
            fprintf(stderr, "error: option -a set more than once.\n");
            say_map_usage();
            exit(EXIT_FAILURE);
         }
         arg_a = 1;
         opt->print_first = 0;
         break;
      case 't':
         if (arg_t) {
            fprintf(stderr, "error: option -t set more than once.\n");
            say_map_usage();
            exit(EXIT_FAILURE);
         }
         arg_t = 1;
         opt->threads = atoi(optarg);
         if (opt->threads < 1) {
            fprintf(stderr, "[error] threads option (-t) must be a positive number.\n");
            exit(EXIT_FAILURE);
         }
         break;
      case 'e':
         if (arg_e) {
            fprintf(stderr, "error: option -e set more than once.\n");
            say_map_usage();
            exit(EXIT_FAILURE);
         }
         arg_e = 1;
         int eval = atoi(optarg);
         if (eval > 0) eval = -eval;
         opt->eval_thr = eval;
         break;
      case 'q':
         if (arg_q) {
            fprintf(stderr, "error: option -q set more than once.\n");
            say_map_usage();
            exit(EXIT_FAILURE);
         }
         arg_q = 1;
         int q = atoi(optarg);
         if (q < 0 || q > 60) {
            fprintf(stderr, "error: quality value must be in interval [0,60].\n");
            exit(EXIT_FAILURE);
         }
         opt->mapq_thr = q;
         break;
      }
   }

   if (optind != argc-2) {
      fprintf(stderr, "error: not enough parameters.\n");
      say_map_usage();
      exit(EXIT_FAILURE);
   }

   // Store index and query files.
   *ifile = argv[optind];
   *qfile = argv[optind+1];

   return 0;
}
