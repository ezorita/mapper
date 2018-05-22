#include "user_interface.h"

char *USAGE_MAP =
   "\n"
   "usage:\n"
   "  mapper [options] index-file input-file\n"
   "\n"
   "  options:\n"
   "    -v --version: prints version.\n"
   /*
   "    -t --threads: number of concurrent threads (default 1)\n"
   "    -a --all:     print all hits found in each interval.\n"
   "    -q --quality: do not print reads with mapping quality below [quality].\n"
   "    -e --eval:    do not print reads with significance below [e-value].\n"
   */
   "\n"
   " to build/edit/extend an index file, type:\n"
   "  \"mapper index -h\"\n";

char *USAGE_ADD =
   "\n"
   "usage:\n"
   " mapper index add [options] index-file\n"
   "\n"
   "  options:\n"
   "    -k --kmer:       sequence length. [required]\n"
   "    -d --distance:   sequence mismatches. [required]\n"
   "    -t --threads:    number of threads. (default 1)\n";

char *USAGE_INDEX =
   "\n"
   "usage:\n"
   " mapper index [command]\n"
   "\n"
   "  commands:\n"
   "    build:  build a new index from scratch.\n"
   "    add:    add annotation data to an existing index.\n"
   //   "    view:   view index annotations.\n"
   "\n"
   "  for more help type: index command\n";

char *USAGE_BUILD =
   "\n"
   "usage:\n"
   " mapper index build [options] genome.fasta\n"
   "\n"
   "  creates a new index for the specified genome file.\n"
   "  the index will be created without neighbor annotations.\n"
   "  To add neighbor annotations to an existing index,\n"
   "  see 'mapper index add -h'.\n"
   "\n"
   "  options:\n"
   "    -o --output:     output filename prefix (default is same as input),\n"
   /*
   "    -k --kmer:       sequence length for default annotation. (default 25)\n"
   "    -d --distance:   sequence mismatches for default annotation. (default 1)\n"
   "    -t --threads:    number of threads. (default 1)\n"
   */
   "\n"
   "  files that will be generated:\n"
   "    genome.fasta.sym:  symbol dictionary\n"
   "    genome.fasta.txt:  genome sequence\n"
   "    genome.fasta.sar:  suffix array\n"
   "    genome.fasta.bwt:  bwt index\n"
   "  annotation data:\n"
   "    genome.fasta.ann.#\n";

void say_map_usage(void) { fprintf(stderr, "%s\n", USAGE_MAP); }
void say_add_usage(void) { fprintf(stderr, "%s\n", USAGE_ADD); }
void say_index_usage(void) { fprintf(stderr, "%s\n", USAGE_INDEX); }
void say_build_usage(void) { fprintf(stderr, "%s\n", USAGE_BUILD); }
void say_version(void) { fprintf(stderr, "mapper version: %s\n", MAPPER_VERSION); }

int  ui_index_build(int, char **, char **, char **);
int  ui_index_add(int, char **, opt_add_t *, char **);
int  ui_map(int, char **, opt_map_t *, char **, char **);


int
ui_parse
(
 int     argc,
 char  * argv[]
)
{
   // Declare variables.
   index_t * index = NULL;

   // No arguments.
   if (argc == 1) {
      say_map_usage();
      return EXIT_SUCCESS;
   }

   /*
   ** ./mapper index
   */
   if (strcmp(argv[1],"index") == 0) {
      // Without params, print help.    
      if (argc == 2 || (argc == 3 && (!strcmp(argv[2],"-h") || !strcmp(argv[2],"--help")))) {
         say_index_usage();
         return EXIT_SUCCESS;
      } 

      /*
      ** ./mapper index add
      */
      else if (strcmp(argv[2], "add") == 0) {
         if (argc == 3) {
            say_add_usage();
            return EXIT_SUCCESS;
         }
         // Parse params.
         char * index_file;
         opt_add_t opt;
         int rval;

         // Return if return value != 0.
         if ((rval = ui_index_add(argc, argv, &opt, &index_file))) {
            return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
         }

         // Load index.
         index = index_read(index_file);
         error_test(index == NULL);
         // Compute annotation.
         error_test(index_ann_new(opt.k, opt.d, opt.threads, index) == -1);
         // Free index mem.
         index_free(index);
         
         return EXIT_SUCCESS;
      } 

      /*
      ** ./mapper index build
      */
      else if (strcmp(argv[2], "build") == 0) {
         if (argc == 3) {
            say_build_usage();
            return EXIT_SUCCESS;
         }
         // Parse params.
         char * genome_file;
         char * output_file;
         int rval;
         if ((rval = ui_index_build(argc, argv, &genome_file, &output_file))) {
            return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
         }

         // Compute index.
         output_file = output_file == NULL ? genome_file : output_file;
         fprintf(stderr, "[index/build] input file: %s, output file: %s\n", genome_file, output_file);
         index = index_build(genome_file, output_file);
         error_test(index == NULL);
         fprintf(stderr, "[index/build] done\n");
         // Free index.
         index_free(index);
         
         return EXIT_SUCCESS;
      }
   } else {
      int rval;
      opt_map_t opt;
      char * index_file;
      char * query_file;

      if ((rval = ui_map(argc, argv, &opt, &index_file, &query_file))) {
         return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
      }

      // Load index.
      index = index_read(index_file);
      error_test(index == NULL);
      // Map query file.
      error_test(mapper(query_file, index) == -1);
      // Free index.
      index_free(index);
         
      return EXIT_SUCCESS;
   }

 failure_return:
   index_free(index);
   return EXIT_FAILURE;
}


int
ui_index_build
(
 int            argc,
 char         * argv[],
 char        ** gfile,
 char        ** ofile
 )
{
   int c;
   *ofile = NULL;
  
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"output",     required_argument,      0, 'o'},
         {"help",       no_argument,            0, 'h'},
         {0, 0, 0, 0}
      };

      c = getopt_long(argc, argv, "ho:", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
      case 'o':
         *ofile = optarg;
         break;
      case 'h':
         say_build_usage();
         return 0;
      }
   }

   if (optind + 2 == argc-1) {
      *gfile = argv[argc-1];
   } else {
      fprintf(stderr, "error: incorrect parameters.\n");
      say_build_usage();
      return -1;
   }

   return 0;
}


int
ui_index_add
(
 int          argc, 
 char       * argv[],
 opt_add_t  * opt,
 char      ** ifile
)
{
   int arg_k = -1, arg_d = -1, arg_t = -1;
   int c;
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"threads",    required_argument,      0, 't'},
         {"kmer",       required_argument,      0, 'k'},
         {"distance",   required_argument,      0, 'd'},
         {"help",       no_argument,            0, 'h'},
         {0, 0, 0, 0}
      };

      c = getopt_long(argc, argv, "t:k:d:h", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
      case 't':
         if (arg_t < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "error: threads option (-t) must be a positive number.\n");
               return -1;
            }
            arg_t = v;
         } else {
            fprintf(stderr,"error: threads option (-t) set more than once.\n");
         }
         break;

      case 'k':
         if (arg_k < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "error: kmer option (-k) must be a positive number.\n");
               return -1;
            }
            arg_k = v;
         } else {
            fprintf(stderr,"error: kmer option (-k) set more than once.\n");
         }
         break;

      case 'd':
         if (arg_d < 0) {
            int v = atoi(optarg);
            if (v < 0) {
               fprintf(stderr, "error: distance option (-d) must be a non-negative number.\n");
               return -1;
            }
            arg_d = v;
         } else {
            fprintf(stderr,"error: distance option (-d) set more than once.\n");
         }
         break;

      case 'h':
         say_add_usage();
         return 1;

      }
   }

   if (optind + 2 != argc-1) {
      fprintf(stderr, "error: not enough parameters.\n");
      say_add_usage();
      return -1;
   }

   *ifile = argv[argc-1];

   if (arg_k < 0 || arg_d < 0) {
      fprintf(stderr, "error: kmer (-k) and distance (-d) are required parameters.\n");
      return -1;
   }

   opt->k          = arg_k;
   opt->d          = arg_d;
   opt->threads    = (arg_t < 0 ? 1 : arg_t);

   return 0;
}

int
ui_map
(
 int         argc,
 char      * argv[],
 opt_map_t * opt,
 char     ** ifile,
 char     ** qfile
)
{
   int arg_t = 0, arg_a = 0, arg_e = 0, arg_q = 0;
   int c;
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         /*
         {"all",      no_argument,       0, 'a'},
         {"threads",  required_argument, 0, 't'},
         {"eval",     required_argument, 0, 'e'},
         {"quality",  required_argument, 0, 'q'},
         */
         {"help",     no_argument,       0, 'h'},
         {"version",  no_argument,       0, 'v'},
         {0, 0, 0, 0}
      };

      c = getopt_long(argc, argv, "ahvt:e:q:",
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
            fprintf(stderr, "error: threads option (-t) must be a positive number.\n");
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

      case 'v':
         say_version();
         return 1;

      case 'h':
         say_map_usage();
         return 1;
      }
   }

   if (optind != argc-2) {
      fprintf(stderr, "error: not enough parameters.\n");
      say_map_usage();
      return -1;
   }

   // Store index and query files.
   *ifile = argv[optind];
   *qfile = argv[optind+1];

   return 0;
}
