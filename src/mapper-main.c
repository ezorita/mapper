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
   signal(SIGSEGV, SIGSEGV_handler); 

   if (argc < 3) {
      fprintf(stderr, "usage: bwmapper {query <file.fastq> | index} <genome file>\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1],"query") == 0) {
      // Parse query filename.
      FILE * queryfile = fopen(argv[2], "r");
      if (queryfile == NULL) {
         fprintf(stderr, "error: could not open file: %s (%s).\n", argv[2],strerror(errno));
         return EXIT_FAILURE;
      }

      // Read FM index format.
      if (opt_verbose) fprintf(stderr, "loading index...\n");
      index_t * index = load_index(argv[3]);
      if(index == NULL) return EXIT_FAILURE;

      // Read file.
      if (opt_verbose) fprintf(stderr, "reading query file...\n");
      seqstack_t * seqs = read_file(queryfile, opt_reverse); // TODO: set reverse and verbose options.
      if (seqs == NULL) {
         return EXIT_FAILURE;
      }

      // ILLUMINA STRATEGY

      // Quality strategy (decresing quality).
      /*
      int rounds = 4;
      int kmers[7] = {24,18,22,18};
      //int kmers[4] = {22,20,18,18};
      int tau[7] = {0,0,0,0,0,0,0};
      int off[7] = {80,80,11,4};
      char qthr[7] = {0,0,0,0,0,0,0};
      //      char qthr[7] = {'F','9','0','0',0};
      */
      int rounds = 4;
      int kmers[7] = {22,18,22,18,20};
      //int kmers[4] = {22,20,18,18};
      int tau[7] = {0,0,0,0,0,0};
      int off[7] = {22,18,11,1,1};
      char qthr[7] = {0,0,0,0,0,0,0};
      /*
      int rounds = 1;
      int kmers[7] = {18,18,22,18};
      //int kmers[4] = {22,20,18,18};
      int tau[7] = {0,0,0,0,0,0,0};
      int off[7] = {1,80,11,4};
      char qthr[7] = {0,0,0,0,0,0,0};
      */

      hmargs_t hmargs = {
         // Format options.
         .verbose = 1,
         .repeat_print_num = 5,
         .sequence_blocks = 100000,
         // Seeding strategy.
         .search_rounds = rounds,
         .tau = tau,
         .kmer_size = kmers,
         .kmer_offset = off,
         .qthr = qthr,
         // Seeding options.
         .seed_max_loci = 20,
         .seed_abs_max_loci = 100,
         // Hitmap analysis options.
         .read_ref_ratio = 1 + 0.05,
         .dist_accept = 20,
         // Alignment filter options.
         .max_align_per_read = 10000,
         .align_filter_ident = 0.9,
         .align_filter_eexp = -10.0,
         .align_accept_eexp = -30.0,
         .align_seed_filter_thr = 0.5,
         // Post-processing options.
         .feedback_eexp_thr = -25.0,
         .feedback_gap_minlen = 5,
         .repeat_min_overlap = 0.9,
         .overlap_tolerance = 0.1,
         .overlap_max_tolerance = 0.5,
         .fuse_min_spanratio = 0.5,
         // Alignment algorithm options.
         .align = (alignopt_t) {
            .bp_thr       = 10,
            .bp_max_thr   = 50,
            .bp_resolution = 1,
            .bp_period    = 5, // This is N times the resolution!!!
            .bp_repeats   = 2,
            .read_error   = 0.05,
            .rand_error   = 0.50,
            .width_ratio  = 0
         }
      };

      /*
      // Hitmap.
      // Default arguments.

      // Decreasing strategy.
      // PACBIO STRATEGY
      //int rounds = 3;
      //int kmers[4] = {22,20,18,22};
      //int tau[4] = {0,0,0,1};
      //char qthr[4] = {0,0,0,0};
      // Quality strategy (decresing seed).
      int rounds = 5;
      int kmers[8] = {22,22,20,20,18,22,20,18};
      int tau[8] = {0,0,0,0,0,0,0,1};
      char qthr[8] = {'(','$','(','$',0,0,0,0};
      // Quality strategy (decresing quality).
      //int rounds = 7;
      //int kmers[8] = {22,20,18,22,20,18,18,18};
      //int tau[8] = {0,0,0,0,0,0,0,0};
      //char qthr[8] = {'(','(','(','$','$','$',0,0};
      // Short-seed quality strategy.
      //int rounds = 3;
      //int kmers[3] = {18,18,18};
      //int tau[3] = {0,0,0};
      //char qthr[3] = {'(','$',0};
      // All-in strategy.
      //int rounds = 1;
      //int kmers[1] = {18};
      //int tau[1] = {0};
      //char qthr[1] = {0};


      hmargs_t hmargs = {
         // Format options.
         .verbose = 1,
         .repeat_print_num = 5,
         // Seeding strategy.
         .seed_strategy = SEED_PARALLEL,
         .search_rounds = rounds,
         .tau = tau,
         .kmer_size = kmers,
         .qthr = qthr,
         .sequence_blocks = 100000,
         // Seeding options.
         .seed_max_loci = 20,
         .seed_abs_max_loci = 10000,
         // Hitmap analysis options.
         .read_ref_ratio = 1 + 0.15,
         .dist_accept = 20,
         // Alignment filter options.
         .max_align_per_read = 10000,
         .align_filter_ident = 0.7,
         .align_filter_eexp = -15.0,
         .align_accept_eexp = -100.0,
         .align_seed_filter_thr = 0.5,
         // Post-processing options.
         .feedback_eexp_thr = -25.0,
         .feedback_gap_minlen = 50,
         .repeat_min_overlap = 0.9,
         .overlap_tolerance = 0.1,
         .overlap_max_tolerance = 0.5,
         .fuse_min_spanratio = 0.5,
         // Alignment algorithm options.
         .align = (alignopt_t) {
            .bp_thr       = 10,
            .bp_max_thr   = 50,
            .bp_period    = 10,
            .bp_repeats   = 20,
            .read_error   = 0.20,
            .rand_error   = 0.50,
            .width_ratio  = 0.15
         }
      };
      */     
      clock_t tstart = clock();
      hitmap(index, seqs, hmargs);
      double totaltime = ((clock()-tstart)*1.0)/CLOCKS_PER_SEC;
      if (opt_verbose) fprintf(stderr, "query time [%.3fs] / rate [%.3f ms/read]\n", totaltime, totaltime*1000.0/seqs->pos);
      
   }
   else if (strcmp(argv[1],"index") == 0) {
      write_index(argv[2]);
   }
   else 
      fprintf(stderr, "usage: bwmapper {query <seq> | index} <genome file>\n");
   
   return 0;
}
