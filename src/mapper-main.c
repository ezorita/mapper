#include "bwmapper.h"
#include <time.h>

void SIGSEGV_handler(int sig) {
   void *array[10];
   size_t size;

   // get void*'s for all entries on the stack
   size = backtrace(array, 10);

   // print out all the frames to stderr
   fprintf(stderr, "Error: signal %d:\n", sig);
   backtrace_symbols_fd(array, size, STDERR_FILENO);
   exit(1);
}


int main(int argc, char * argv[]) {

   // TODO: parametrize
   int opt_reverse  = 1;
   int opt_verbose  = 1;

   // Backtrace handler
   //    signal(SIGSEGV, SIGSEGV_handler); 

   if (argc < 3) {
      fprintf(stderr, "usage: bwmapper {query <file.fastq> | index} <genome file>\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1],"query") == 0) {
      char * filename = malloc(strlen(argv[3])+5);
      strcpy(filename, argv[3]);

      // Open index file.
      strcpy(filename+strlen(argv[3]), ".fmi");
      int fd = open(filename, O_RDONLY);
      if (fd == -1) {
         fprintf(stderr, "error opening '%s': %s\n", filename, strerror(errno));
         exit(EXIT_FAILURE);
      }

      // Parse query filename.
      FILE * queryfile = fopen(argv[2], "r");
      if (queryfile == NULL) {
         fprintf(stderr, "error: could not open file.\n");
         return EXIT_FAILURE;
      }

      // Read file.
      if (opt_verbose) fprintf(stderr, "reading query file...\n");
      seqstack_t * seqs = read_file(queryfile, opt_reverse); // TODO: set reverse and verbose options.
      if (seqs == NULL) {
         return EXIT_FAILURE;
      }

      // Load index.
      int mflags = MAP_PRIVATE | MAP_POPULATE;
      long idxsize = lseek(fd, 0, SEEK_END);
      lseek(fd, 0, SEEK_SET);
      long * indexp = mmap(NULL, idxsize, PROT_READ, mflags, fd, 0);
      if (indexp == NULL) {
         fprintf(stderr, "error opening index file: %s.\n", strerror(errno));
         return EXIT_FAILURE;
      }
      
      // Read FM index format.
      index_t index;
      if(format_FMindex(indexp, &index))
         return EXIT_FAILURE;

      // Read chromosome index.
      strcpy(filename+strlen(argv[3]), ".index");
      chr_t * chr = read_CHRindex(filename);
      if (chr == NULL) return EXIT_FAILURE;

      // Hitmap.
      // Default arguments.

      int rounds = 4;
      int kmers[6] = {22,20,18,22,20,18};
      int tau[6] = {0,0,0,1,0,1};
      char qthr[6] = {0,0,0,0,0,0};

      hmargs_t hmargs = {
         // Format options.
         .verbose = 1,
         .repeat_print_num = 5,
         // Seeding strategy.
         .search_rounds = rounds,
         .tau = tau,
         .kmer_size = kmers,
         .qthr = qthr,
         // Seeding options.
         .seed_max_loci = 20,
         .seed_abs_max_loci = 10000,
         // Hitmap analysis options.
         .read_ref_ratio = 1 + 0.15,
         .dist_accept = 20,
         // Alignment filter options.
         .max_align_per_read = 10000,
         .align_full_seed_thr = 0.01,
         .align_filter_eexp = -6.0,
         .align_accept_eexp = -50.0,
         .align_seed_filter_thr = 0.5,
         // Post-processing options.
         .feedback_eexp_thr = -9.0,
         .feedback_gap_minlen = 75,
         .repeat_min_overlap = 0.9,
         .overlap_tolerance = 0.1,
         .overlap_max_tolerance = 0.5,
         .fuse_min_spanratio = 0.5,
         // Alignment algorithm options.
         .align = (alignopt_t) {
            .border_slope = 0.168,
            .border_y0    = 26,
            .bp_thr       = 30,
            .bp_max_thr   = 25,
            .bp_min_thr   = 5,
            .bp_period    = 20,
            .bp_repeats   = 3,
            .read_error   = 0.10,
            .rand_error   = 0.45,
            .read_subst   = 0.10,
            .rand_subst   = 0.30
         }
      };
      clock_t tstart = clock();
      hitmap(&index, chr, seqs, hmargs);
      double totaltime = ((clock()-tstart)*1.0)/CLOCKS_PER_SEC;
      if (opt_verbose) fprintf(stderr, "query time [%.3fs] / rate [%.3f s/read]\n", totaltime, totaltime/seqs->pos);
      
   }
   else if (strcmp(argv[1],"index") == 0) {
      write_index(argv[2]);
   }
   else 
      fprintf(stderr, "usage: bwmapper {query <seq> | index} <genome file>\n");
   
   return 0;
}
