#include "user_interface.h"

char *ERROR_MAX_ARG = "error: too many arguments.\n";
char *ERROR_INSUF_ARG = "error: not enough arguments.\n";
char *ERROR_COMMAND = "error: incorrect command.\n";
char *ERROR_INCORRECT_OPT = "error: incorrect option - ";

char *OPT_THREAD_POSITIVE = "error: threads option (-t) must be a positive number.\n";
char *OPT_THREAD_REPEAT   = "error: threads option (-t) set more than once.\n";
char *OPT_KMER_POSITIVE   = "error: kmer option (-k) must be a positive number.\n";
char *OPT_KMER_REPEAT     = "error: kmer option (-k) set more than once.\n";
char *OPT_DIST_POSITIVE   = "error: distance option (-d) must be a non-negative number.\n";
char *OPT_DIST_REPEAT     = "error: distance option (-d) set more than once.\n;";
char *OPT_OUTPUT_REPEAT   = "error: output option (-o) set more than once.\n";

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
   "    view:   view index information.\n"
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

char *USAGE_VIEW =
   "\n"
   "usage:\n"
   " mapper index view <index-name>\n"
   "\n"
   "  displays index information.\n";

void say_map_usage(void) { fprintf(stderr, "%s\n", USAGE_MAP); }
void say_add_usage(void) { fprintf(stderr, "%s\n", USAGE_ADD); }
void say_index_usage(void) { fprintf(stderr, "%s\n", USAGE_INDEX); }
void say_build_usage(void) { fprintf(stderr, "%s\n", USAGE_BUILD); }
void say_view_usage(void) { fprintf(stderr, "%s\n", USAGE_VIEW); }
void say_version(void) { fprintf(stderr, "mapper version: %s\n", MAPPER_VERSION); }

int
ui_parse
(
 int     argc,
 char  * argv[]
)
{
   // Declare variables.
   index_t * index = NULL;
   char * index_path = NULL;

   // Check arguments.
   error_test_msg(argc <= 0, "argument 'argc' must be greater than 0.");
   error_test_msg(argv == NULL, "argument 'argv' is NULL.");

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

         index_path = get_filename_base(index_file);
         // Load index.
         index = index_read(index_path);
         error_test(index == NULL);
         // Compute annotation.
         error_test(index_ann_new(opt.k, opt.d, opt.threads, index) == -1);
         // Free index mem.
         index_free(index);
         free(index_path);
         
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
         char * genome_file = NULL;
         char * output_file = NULL;
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
      
      /*
      ** ./mapper index view
      */
      else if (strcmp(argv[2], "view") == 0) {
         if (argc == 3) {
            say_view_usage();
            return EXIT_SUCCESS;
         } else if (argc != 4) {
            fprintf(stderr, "%s", ERROR_MAX_ARG);
            say_view_usage();
            return EXIT_FAILURE;
         }
         
         if (strcmp(argv[3],"-h") == 0 || strcmp(argv[3],"--help") == 0) {
            say_view_usage();
            return EXIT_SUCCESS;
         }

         // Parse params.
         char * index_file = argv[3];
         index_path = get_filename_base(index_file);
         // Read index file.
         index = index_read(index_path);
         error_test(index == NULL);
         // Show index info.
         ui_index_info(index);
         // Free index.
         index_free(index);
         free(index_path);

         return EXIT_SUCCESS;
      } else {
         fprintf(stderr, "%s", ERROR_COMMAND);
            say_index_usage();
            return EXIT_FAILURE;
      }
   } else {
      int rval;
      opt_map_t opt;
      char * index_file;
      char * query_file;

      if ((rval = ui_map(argc, argv, &opt, &index_file, &query_file))) {
         return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
      }

      index_path = get_filename_base(index_file);

      // Load index.
      index = index_read(index_path);
      error_test(index == NULL);
      // Map query file.
      error_test(mapper(query_file, index) == -1);
      // Free index.
      index_free(index);
      free(index_path);
         
      return EXIT_SUCCESS;
   }

 failure_return:
   index_free(index);
   free(index_path);
   return EXIT_FAILURE;
}

int
ui_index_info
(
 index_t  * index
)
{
   // Check arguments.
   error_test_msg(index == NULL, "argument 'index' is NULL.");

   // Index basic info.
   fprintf(stderr, "[basic info]\n");
   fprintf(stderr, " index basename:   %s\n", index->fname_base);
   fprintf(stderr, " index structures:\n");
   fprintf(stderr, "  symbols info:    %s\n", (index->sym == NULL ? "NO" : "YES"));
   fprintf(stderr, "  reference text:  %s\n", (index->txt == NULL ? "NO" : "YES"));
   fprintf(stderr, "  suffix array:    %s\n", (index->sar == NULL ? "NO" : "YES"));
   fprintf(stderr, "  FM index:        %s\n", (index->bwt == NULL ? "NO" : "YES"));
   fprintf(stderr, "  annotations:     %s\n", (index->ann_cnt == 0 ? "NO" : "YES"));

   // Symbols info.
   if (index->sym != NULL) {
      int s_cnt = sym_count(index->sym);
      fprintf(stderr, "\n[index symbols]\n");
      fprintf(stderr, " path:             %s.sym\n", index->fname_base);
      fprintf(stderr, " symbol count:     %d\n", s_cnt);
      fprintf(stderr, " alphabet:         { ");
      for (int i = 0; i < s_cnt; i++) {
         fprintf(stderr, "%c ", sym_character(i, index->sym));
      }
      fprintf(stderr, "}\n");
      fprintf(stderr, " complement rel.:  { ");
      for (int i = 0; i < s_cnt; i++) {
         fprintf(stderr, "%c->%c ", sym_character(i, index->sym), \
                 sym_character(sym_complement(i,index->sym), index->sym));
      }
      fprintf(stderr, "}\n");
   }

   // Text info.
   if (index->txt != NULL) {
      fprintf(stderr, "\n[reference text]\n");
      fprintf(stderr, " path:             %s.txt\n", index->fname_base);
      fprintf(stderr, " bidirectional:    %s\n", (txt_has_rc(index->txt) ? "YES" : "NO"));
      fprintf(stderr, " text length:      %ld\n", txt_length(index->txt));
      fprintf(stderr, " sequence count:   %ld\n", txt_seq_count(index->txt));
      fprintf(stderr, " sequences (id, name, length):\n");
      for (int i = 0; i < txt_seq_count(index->txt); i++) {
         fprintf(stderr, "  %d. %s\t%ld\n", i, txt_seq_name(i,index->txt), txt_seq_length(i,index->txt));
      }
   }

   // SAR info.
   if (index->sar != NULL) {
      fprintf(stderr, "\n[suffix array]\n");
      fprintf(stderr, " path:             %s.sar\n", index->fname_base);
   }

   // BWT info.
   if (index->sar != NULL) {
      fprintf(stderr, "\n[FM index]\n");
      fprintf(stderr, " path:             %s.bwt\n", index->fname_base);
   }


   // Annotation info.
   if (index->ann_cnt > 0) {
      fprintf(stderr, "\n[annotations]\n");
      fprintf(stderr, " path:             %s.ann.#.#\n", index->fname_base);
      fprintf(stderr, " annotation count: %d", index->ann_cnt);
      fprintf(stderr, " annotations: (id, kmer, distance):\n");
      for (int i = 0; i < index->ann_cnt; i++) {
         fprintf(stderr, "   %d. (%d,%d)\n", i, ann_get_kmer(index->ann[i]), ann_get_dist(index->ann[i]));
      }
   }

   return 0;

 failure_return:
   return -1;
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
   // Check arguments.
   error_test_msg(argc <= 0, "argument 'argc' must be greater than 0.");
   error_test_msg(argv == NULL, "argument 'argv' is NULL.");
   error_test_msg(gfile == NULL, "argument '*gfile' is NULL.");
   error_test_msg(ofile == NULL, "argument '*ofile' is NULL.");

   int c;
   *ofile = NULL;
  
   while (1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"output",     required_argument,      0, 'o'},
         {"help",       no_argument,            0, 'h'},
         {0, 0, 0, 0}
      };

      opterr = 0;
      c = getopt_long(argc, argv, "ho:", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
      case '?':
         fprintf(stderr, "%s%c\n", ERROR_INCORRECT_OPT, optopt);
         say_build_usage();
         return -1;

      case 'o':
         if (*ofile != NULL) {
            fprintf(stderr, "%s", OPT_OUTPUT_REPEAT);
            return -1;
         }
         *ofile = optarg;
         break;
      case 'h':
         say_build_usage();
         return 1;
      }
   }

   if (optind + 2 == argc-1) {
      *gfile = argv[argc-1];
   } else if (optind + 2 < argc - 1) {
      fprintf(stderr, "%s", ERROR_MAX_ARG);
      say_build_usage();
      return -1;
   } else if (optind + 2 > argc - 1) {
      fprintf(stderr, "%s", ERROR_INSUF_ARG);
      say_build_usage();
      return -1;
   }

   return 0;

 failure_return:
   return -1;
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
   // Check arguments.
   error_test_msg(argc <= 0, "argument 'argc' must be greater than 0.");
   error_test_msg(argv == NULL, "argument 'argv' is NULL.");
   error_test_msg(opt == NULL, "argument 'opt' is NULL.");
   error_test_msg(ifile == NULL, "argument 'ifile' is NULL.");

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

      opterr = 0;
      c = getopt_long(argc, argv, "t:k:d:h", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
      case '?':
         fprintf(stderr, "%s%c\n", ERROR_INCORRECT_OPT, optopt);
         say_add_usage();
         return -1;
      case 't':
         if (arg_t < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "%s", OPT_THREAD_POSITIVE);
               return -1;
            }
            arg_t = v;
         } else {
            fprintf(stderr,"%s", OPT_THREAD_REPEAT);
         }
         break;

      case 'k':
         if (arg_k < 0) {
            int v = atoi(optarg);
            if (v <= 0) {
               fprintf(stderr, "%s", OPT_KMER_POSITIVE);
               return -1;
            }
            arg_k = v;
         } else {
            fprintf(stderr, "%s", OPT_KMER_REPEAT);
         }
         break;

      case 'd':
         if (arg_d < 0) {
            int v = atoi(optarg);
            if (v < 0) {
               fprintf(stderr, "%s", OPT_DIST_POSITIVE);
               return -1;
            }
            arg_d = v;
         } else {
            fprintf(stderr, "%s", OPT_DIST_REPEAT);
         }
         break;

      case 'h':
         say_add_usage();
         return 1;

      }
   }

   if (argc-1 < optind + 2 || arg_k < 0 || arg_d < 0) {
      fprintf(stderr, "%s", ERROR_INSUF_ARG);
      say_add_usage();
      return -1;
   } else if (argc-1 > optind + 2) {
      fprintf(stderr, "%s", ERROR_MAX_ARG);
      say_add_usage();
      return -1;
   }

   *ifile = argv[argc-1];

   opt->k          = arg_k;
   opt->d          = arg_d;
   opt->threads    = (arg_t < 0 ? 1 : arg_t);

   return 0;

 failure_return:
   return -1;
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
   // Check arguments.
   error_test_msg(argc <= 0, "argument 'argc' must be greater than 0.");
   error_test_msg(argv == NULL, "argument 'argv' is NULL.");
   error_test_msg(opt == NULL, "argument 'opt' is NULL.");
   error_test_msg(ifile == NULL, "argument 'ifile' is NULL.");
   error_test_msg(qfile == NULL, "argument 'qfile' is NULL.");

   //   int arg_t = 0, arg_a = 0, arg_e = 0, arg_q = 0;

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

      opterr = 0;
      c = getopt_long(argc, argv, "hv", long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case '?':
         fprintf(stderr, "%s%c\n", ERROR_INCORRECT_OPT, optopt);
         say_map_usage();
         return -1;
         /*
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
            fprintf(stderr, "%s", OPT_THREAD_REPEAT);
            say_map_usage();
            exit(EXIT_FAILURE);
         }
         arg_t = 1;
         opt->threads = atoi(optarg);
         if (opt->threads < 1) {
            fprintf(stderr, "%s", OPT_THREAD_POSITIVE);
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
         */
      case 'v':
         say_version();
         return 1;

      case 'h':
         say_map_usage();
         return 1;
      }
   }

   if (argc-2 < optind) {
      fprintf(stderr, "%s", ERROR_INSUF_ARG);
      say_map_usage();
      return -1;
   } else if (argc-2 > optind) {
      fprintf(stderr, "%s", ERROR_MAX_ARG);
      say_map_usage();
      return -1;
   }

   // Store index and query files.
   *ifile = argv[optind];
   *qfile = argv[optind+1];

   return 0;

 failure_return:
   return -1;
}

char *
get_filename_base
(
 char  * path
)
{
   char * index_path = NULL;
   
   error_test_msg(path == NULL, "argument 'index_path' is NULL.");

   index_path = malloc(strlen(path) + 1);
   error_test_mem(index_path);
   memcpy(index_path, path, strlen(path)+1);


   char * ext = index_path + strlen(index_path) - 4;

   if ((strcmp(ext, ".sym") == 0) || \
       (strcmp(ext, ".txt") == 0) || \
       (strcmp(ext, ".bwt") == 0) || \
       (strcmp(ext, ".sar") == 0) || \
       (strcmp(ext, ".ann") == 0)) {
      index_path[strlen(index_path) - 4] = 0;
   }

   return index_path;
   
 failure_return:
   free(index_path);
   return NULL;
}
