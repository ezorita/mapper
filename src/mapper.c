#include "mapper.h"
#include <getopt.h>
#include <time.h>

char *USAGE =
   "\n"
   "Usage:"
   "  mapper [options] input-file genome-index\n"
   "\n"
   "  options:\n"
   "    -t --threads: number of concurrent threads (default 1)\n"
   "    -a --all:     print all the hits found in each interval.\n"
   "    -q --quality: do not print reads with mapping quality below [quality].\n"
   "    -e --eval:    do not print reads with significance below [e-value].\n";


void say_usage(void) { fprintf(stderr, "%s\n", USAGE); }

int parse_opt(int argc, char *argv[], mapopt_t * opt, char ** qfile, char ** ifile) {
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
            say_usage();
            exit(EXIT_FAILURE);
         }
         arg_a = 1;
         opt->format.print_first = 0;
         break;
      case 't':
         if (arg_t) {
            fprintf(stderr, "error: option -t set more than once.\n");
            say_usage();
            exit(EXIT_FAILURE);
         }
         arg_t = 1;
         opt->threads = atoi(optarg);
         break;
      case 'e':
         if (arg_e) {
            fprintf(stderr, "error: option -e set more than once.\n");
            say_usage();
            exit(EXIT_FAILURE);
         }
         arg_e = 1;
         int eval = atoi(optarg);
         if (eval > 0) eval = -eval;
         opt->format.eval_thr = eval;
         break;
      case 'q':
         if (arg_q) {
            fprintf(stderr, "error: option -q set more than once.\n");
            say_usage();
            exit(EXIT_FAILURE);
         }
         arg_q = 1;
         int q = atoi(optarg);
         if (q < 0 || q > 60) {
            fprintf(stderr, "error: quality value must be in interval [0,60].\n");
            exit(EXIT_FAILURE);
         }
         opt->format.mapq_thr = q;
         break;
      }
   }

   if (optind != argc-2) {
      fprintf(stderr, "error: not enough parameters.\n");
      say_usage();
      exit(EXIT_FAILURE);
   }

   // Store index and query files.
   *qfile = argv[optind];
   *ifile = argv[optind+1];

   return 0;
}


int main(int argc, char *argv[])
{
   int opt_verbose = 1;
   size_t seqs_per_thread = 10000;

   // DEFAULT Options.
   seedopt_t seedopt = { // MEM
      .min_len = 18,
      .max_len = 1000,
      .min_loci = 1,
      .max_loci = 200,
      .aux_loci = 1000,
      .thr_seed = 20,
      .reseed_len = 28
   };
   filteropt_t filteropt = {
      .dist_accept = 10,
      .max_align_per_read = 200,
      .read_ref_ratio = 0.05,
      .align_accept_eexp = -100.0,
      .overlap_max_tolerance = 0.5,
      .align_seed_filter_thr = 0.5,
      .align_seed_filter_dif = 19,//38,
      .align_filter_ident = 0.85,
      .align_filter_eexp = 0.0,
      .mapq_evalue_ratio = 0.5
   };
   alignopt_t alignopt = {
      .bp_diagonal  = 1,
      .bp_thr       = 10,
      .bp_max_thr   = 50,
      .bp_resolution = 1,
      .bp_period    = 5,
      .bp_repeats   = 2,
      .read_error   = 0.02,
      .rand_error   = 0.50,
      .width_ratio  = 0.05
   };
   formatopt_t formatopt = {
      .print_first = 2, // 2 for first_only
      .mapq_thr = 0,
      .eval_thr = 0
   };
   mapopt_t mapopt = {
      .seed = seedopt,
      .align = alignopt,
      .filter = filteropt,
      .format = formatopt,
      .threads = 1
   };

   // Parse input parameters.
   if (argc == 1) {
      say_usage();
      exit(EXIT_SUCCESS);
   }
   char * qfile, *ifile;
   parse_opt(argc, argv, &mapopt, &qfile, &ifile);



   // Parse query filename.
   FILE * queryfile = fopen(qfile, "r");
   if (queryfile == NULL) {
      fprintf(stderr, "error: could not open file: %s (%s).\n", argv[1],strerror(errno));
      return EXIT_FAILURE;
   }

   // Read FM index format.
   if (opt_verbose) fprintf(stderr, "loading index...\n");
   idxfiles_t * files = index_open(ifile);
   index_t * index = index_format(files);
   if (index == NULL)
      return EXIT_FAILURE;

   // Read file.
   if (opt_verbose) fprintf(stderr, "reading query file...\n");
   seqstack_t * seqs = read_file(queryfile);
   if (seqs == NULL) {
      return EXIT_FAILURE;
   }
   fclose(queryfile);

   // Run thread scheduler.
   mt_scheduler(seqs, &mapopt, index, mapopt.threads, seqs_per_thread);

   // Free memory.
   free(seqs);
   free(index->lcp_sample);
   free(index->lcp_extend);
   free(index);
   free(files->chr->start);
   for (int i = 0; i < files->chr->nchr; i++) free(files->chr->name[i]);
   free(files->chr->name);
   free(files->chr);
   munmap(files->gen_file, files->gen_len);
   munmap(files->occ_file, files->occ_len);
   munmap(files->lcp_file, files->lcp_len);
   munmap(files->sa_file, files->sa_len);
   free(files);

   return 0;
}

int
mt_scheduler
(
 seqstack_t * seqs,
 mapopt_t   * options,
 index_t    * index,
 int          threads,
 size_t       thread_seq_count
)
{
   // Create control struct.
   mtcontrol_t * control = malloc(sizeof(mtcontrol_t));
   control->mapped = 0;
   control->active = 0;

   // Create mutex and monitor.
   pthread_mutex_t * mutex = malloc(sizeof(pthread_mutex_t));
   pthread_mutex_init(mutex,NULL);
   pthread_cond_t  * monitor = malloc(sizeof(pthread_cond_t));
   pthread_cond_init(monitor,NULL);

   // Process reads.
   clock_t t = clock();

   // Lock mutex.
   for (size_t i = 0; i < seqs->pos; i += thread_seq_count) {
      pthread_mutex_lock(mutex);
      // Verbose.
      if(i%1000 == 0) fprintf(stderr, "progress: %.2f%%\r",i*100.0/seqs->pos);
      // Alloc job (will be freed by thread).
      mtjob_t * job = malloc(sizeof(mtjob_t));
      // Set job.
      job->count   = min(seqs->pos-i,thread_seq_count);
      job->seq     = seqs->seq + i;
      job->index   = index;
      job->opt     = options;
      job->mutex   = mutex;
      job->monitor = monitor;
      job->control = control;
      // Increase thread count.
      control->active += 1;
      // Create thread.

      pthread_t thread;
      if (pthread_create(&thread, NULL, mt_worker, job)) {
         fprintf(stderr, "error 'pthread_create'\n");
         return 1;
      }
      // Deatch thread.
      pthread_detach(thread);
      // Sleep.
      while (control->active >= threads) {
         pthread_cond_wait(monitor, mutex);
      }
      pthread_mutex_unlock(mutex);
   }

   // Wait for remaining threads.
   pthread_mutex_lock(mutex);
   while (control->active > 0) {
      pthread_cond_wait(monitor, mutex);
   }
   pthread_mutex_unlock(mutex);

   // Verbose.
   fprintf(stderr, "progress: 100%% (%.2f%% mapped) in %.3fs\n",control->mapped*100.0/seqs->pos,(clock()-t)*1.0/CLOCKS_PER_SEC);   

   // Free structs.
   free(control);
   free(mutex);
   free(monitor);
   return 0;
}


void *
mt_worker
(
 void * args
)
{
   mtjob_t * job = (mtjob_t *) args;
   // Get params.
   seq_t    * seq   = job->seq;
   mapopt_t * opt   = job->opt;
   index_t  * index = job->index;

   // Alloc data structures.
   matchlist_t * seed_matches = matchlist_new(opt->filter.max_align_per_read);
   matchlist_t * map_matches = matchlist_new(64);
   matchlist_t * repeats = matchlist_new(64);
   seedstack_t * seeds  = seedstack_new(SEEDSTACK_SIZE);
   
   size_t mapped = 0;
   for (size_t i = 0; i < job->count; i++) {
      if (VERBOSE_DEBUG) fprintf(stdout, "seq: %s\n", seq[i].seq);
      size_t slen = strlen(seq[i].seq);
      // Reset seeds.
      seeds->pos = 0;
      // Debug.
      size_t n_mems;
      // Seed MEM.
      seed_mem(seq[i].seq, 0, slen, &seeds, opt->seed, index);
      // DEBUG MEMs.
      if (VERBOSE_DEBUG) {
         fprintf(stdout, "MEMs: %ld\n", seeds->pos);
         for (int i = 0; i < seeds->pos; i++) {
            seed_t seed = seeds->seed[i];
            fprintf(stdout, "seed: p=%d,l=%ld,d=%d,b=%d\n", seed.qry_pos, seed.ref_pos.ep - seed.ref_pos.sp + 1, seed.ref_pos.depth, seed.bulk);
         }
      }
      n_mems = seeds->pos;
      // Reseed.
      reseed_mem(seq[i].seq,&seeds, opt->seed, index);
      // DEBUG reseedMEMs.
      if (VERBOSE_DEBUG) {
         fprintf(stdout, "Reseed MEMs: %ld\n", seeds->pos - n_mems);
         for (int i = n_mems; i < seeds->pos; i++) {
            seed_t seed = seeds->seed[i];
            fprintf(stdout, "seed: p=%d,l=%ld,d=%d,b=%d\n", seed.qry_pos, seed.ref_pos.ep - seed.ref_pos.sp + 1, seed.ref_pos.depth, seed.bulk);
         }
      }
      n_mems = seeds->pos;
      // Reset repeats.
      repeats->pos = 0;
      // Find repeated regions in read.
      find_repeats(seeds, &repeats, opt->seed.max_loci);
      // Adaptive seeds.
      seed_thr(seq[i].seq, slen, &seeds, opt->seed, index);
      // DEBUG SEEDS.
      if (VERBOSE_DEBUG) {
         fprintf(stdout, "Threshold: %ld\n", seeds->pos - n_mems);
         for (int i = n_mems; i < seeds->pos; i++) {
            seed_t seed = seeds->seed[i];
            fprintf(stdout, "seed: p=%d,l=%ld,d=%d,b=%d\n", seed.qry_pos, seed.ref_pos.ep - seed.ref_pos.sp + 1, seed.ref_pos.depth, seed.bulk);
         }
         fprintf(stdout, "Repeats found: %d\n", repeats->pos);
      }
      // Match seeds.
      seed_matches->pos = 0;
      chain_seeds(seeds, slen, seed_matches, index, opt->filter.dist_accept, opt->filter.read_ref_ratio);
      // Align seeds.
      map_matches->pos = 0;
      align_seeds(seq[i].seq, seed_matches, &map_matches, index, opt->filter, opt->align);
      // Merge intervals.
      int32_t n_ints;
      matchlist_t ** intervals;
      intervals = merge_intervals(map_matches, 0.8, &n_ints);
      // Compute map qualities.
      compute_mapq(intervals, n_ints, opt->filter.mapq_evalue_ratio, 0.5, seq[i], index);
      // Repeat penalty.
      if (repeats->pos > 0) filter_repeats(intervals, repeats, n_ints);
      // Print matches.
      mapped += print_and_free(seq[i], intervals, n_ints, index, opt->format);
      // Free structs.
      free(intervals);
   }

   // Free data structures.
   free(seeds);
   free(seed_matches);
   free(map_matches);
   
   // Report to scheduler.
   pthread_mutex_lock(job->mutex);
   job->control->mapped += mapped;
   job->control->active -= 1;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free args.
   free(args);


   return NULL;
}

idxfiles_t *
index_open
(
 char * file
)
{
   // Alloc index struct.
   idxfiles_t * files = malloc(sizeof(idxfiles_t));
   
   // Read chromosome index.
   char * chr_file = malloc(strlen(file)+5);
   strcpy(chr_file, file);
   strcpy(chr_file+strlen(file), ".chr");
   files->chr = read_CHRindex(chr_file);
   if (files->chr == NULL) {
      fprintf(stderr, "error opening '%s': %s\n", chr_file, strerror(errno));      
      return NULL;
   }
   free(chr_file);

   // Open OCC file.
   char * occ_file = malloc(strlen(file)+5);
   strcpy(occ_file, file);
   strcpy(occ_file+strlen(file), ".occ");
   int fd_occ = open(occ_file, O_RDONLY);
   if (fd_occ == -1) {
      fprintf(stderr, "error opening '%s': %s\n", occ_file, strerror(errno));
      return NULL;
   }
   free(occ_file);

   // Open SA file.
   char * sa_file = malloc(strlen(file)+5);
   strcpy(sa_file, file);
   strcpy(sa_file+strlen(file), ".sar");
   int fd_sa = open(sa_file, O_RDONLY);
   if (fd_sa == -1) {
      fprintf(stderr, "error opening '%s': %s\n", sa_file, strerror(errno));
      return NULL;
   }
   free(sa_file);

   // Open GEN file.
   char * gen_file = malloc(strlen(file)+5);
   strcpy(gen_file, file);
   strcpy(gen_file+strlen(file), ".gen");
   int fd_gen = open(gen_file, O_RDONLY);
   if (fd_gen == -1) {
      fprintf(stderr, "error opening '%s': %s\n", gen_file, strerror(errno));
      return NULL;
   }
   free(gen_file);

   // Open LCP file.
   char * lcp_file = malloc(strlen(file)+5);
   strcpy(lcp_file, file);
   strcpy(lcp_file+strlen(file), ".lcp");
   int fd_lcp = open(lcp_file, O_RDONLY);
   if (fd_lcp == -1) {
      fprintf(stderr, "error opening '%s': %s\n", lcp_file, strerror(errno));
      return NULL;
   }
   free(lcp_file);

   // Open LUT file.
   /*
   char * lut_file = malloc(strlen(file)+5);
   strcpy(lut_file, file);
   strcpy(lut_file+strlen(file), ".lut");
   int fd_lut = open(lut_file, O_RDONLY);
   if (fd_lut == -1) {
      fprintf(stderr, "error opening '%s': %s\n", lut_file, strerror(errno));
      return NULL;
   }
   free(lut_file);
   */
   // Load OCC index.
   files->occ_len = lseek(fd_occ, 0, SEEK_END);
   lseek(fd_occ, 0, SEEK_SET);
   files->occ_file = mmap(NULL, files->occ_len, PROT_READ, MMAP_FLAGS, fd_occ, 0);
   close(fd_occ);
   if (files->occ_file == NULL) {
      fprintf(stderr, "error mmaping .occ index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load SA index.
   files->sa_len = lseek(fd_sa, 0, SEEK_END);
   lseek(fd_sa, 0, SEEK_SET);
   files->sa_file = mmap(NULL, files->sa_len, PROT_READ, MMAP_FLAGS, fd_sa, 0);
   close(fd_sa);
   if (files->sa_file == NULL) {
      fprintf(stderr, "error mmaping .sar index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load GEN index.
   files->gen_len = lseek(fd_gen, 0, SEEK_END);
   lseek(fd_gen, 0, SEEK_SET);
   files->gen_file = mmap(NULL, files->gen_len, PROT_READ, MMAP_FLAGS, fd_gen, 0);
   close(fd_gen);
   if (files->gen_file == NULL) {
      fprintf(stderr, "error mmaping .gen index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load LCP index.
   files->lcp_len = lseek(fd_lcp, 0, SEEK_END);
   lseek(fd_lcp, 0, SEEK_SET);
   files->lcp_file = mmap(NULL, files->lcp_len, PROT_READ, MMAP_FLAGS, fd_lcp, 0);
   close(fd_lcp);
   if (files->lcp_file == NULL) {
      fprintf(stderr, "error mmaping .lcp index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load LUT index.
   /*
   idxsize = lseek(fd_lut, 0, SEEK_END);
   lseek(fd_lut, 0, SEEK_SET);
   files->lut_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_lut, 0);
   close(fd_lut);
   if (files->lut_file == NULL) {
      fprintf(stderr, "error mmaping .lut index file: %s.\n", strerror(errno));
      return NULL;
   }
   */
   return files;

}

index_t *
index_format
(
 idxfiles_t * files
)
{
   if (files == NULL) return NULL;
   // Alloc index struct.
   index_t * index = malloc(sizeof(index_t));
   // GEN.
   index->genome = (char *) files->gen_file;
   // OCC.
   index->occ_mark_int = *(uint64_t *) files->occ_file;
   index->c = ((uint64_t *) files->occ_file + 1);
   index->occ = index->c + NUM_BASES + 1;
   // Genome size.
   index->size = index->c[NUM_BASES]; // Count forward strand only.
   // SAR.
   index->sa_bits = *((uint64_t *) files->sa_file);
   index->sa      = ((uint64_t *) files->sa_file) + 1;
   // LUT.
   //   index->lut_kmer = *((uint64_t *) files->lut_file);
   //   index->lut = ((uint64_t *) files->lut_file) + 1;
   // LCP.
   // params
   index->lcp_mark_int = *((uint64_t *)files->lcp_file);
   index->lcp_min_depth =  *((uint64_t *)files->lcp_file + 1);
   uint64_t lcp_idx_size = *((uint64_t *)files->lcp_file + 2);//((2*index->size + LCP_WORD_SIZE - 1)/LCP_WORD_SIZE + index->lcp_mark_int - 1)/index->lcp_mark_int * (index->lcp_mark_int + 1) + 1;
   // index
   index->lcp_sample_idx = ((uint64_t *) files->lcp_file + 3);
   uint64_t ext_idx_size = *(index->lcp_sample_idx + lcp_idx_size);
   index->lcp_extend_idx = index->lcp_sample_idx + lcp_idx_size + 1;
   // alloc structures
   index->lcp_sample = malloc(sizeof(lcpdata_t));
   if (index->lcp_sample == NULL) return NULL;
   index->lcp_extend = malloc(sizeof(list64_t));
   if (index->lcp_extend == NULL) return NULL;
   // data
   index->lcp_sample->size = *(index->lcp_extend_idx + ext_idx_size)/2;
   index->lcp_sample->lcp = (lcpval_t *)(index->lcp_extend_idx + ext_idx_size + 1);
   uint64_t * lcpext_size = (uint64_t *)(index->lcp_sample->lcp + index->lcp_sample->size);
   index->lcp_extend->size = *lcpext_size;
   index->lcp_extend->val = (int64_t *)(lcpext_size + 1);
   //CHR.
   index->chr = files->chr;
   

   return index;
}

chr_t *
read_CHRindex
(
 char * filename
)
{
   // Files
   FILE * input = fopen(filename,"r");
   if (input == NULL) {
      fprintf(stderr, "error in 'read_CHRindex' (fopen): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   int chrcount = 0;
   int structsize = CHRSTR_SIZE;

   long *  start = malloc(structsize*sizeof(long));
   char ** names = malloc(structsize*sizeof(char*));

   // File read vars
   ssize_t rlen;
   size_t  sz = BUFFER_SIZE;
   char  * buffer = malloc(BUFFER_SIZE);
   int lineno = 0;

   // Read chromosome entries.
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      lineno++;
      char *name;
      int i = 0;
      while (buffer[i] != '\t' && buffer[i] != '\n') i++;
      if (buffer[i] == '\n') {
         fprintf(stderr, "illegal format in %s (line %d) - ignoring chromosome: %s\n", filename, lineno, buffer);
         continue;
      }
      
      if (buffer[rlen-1] == '\n') buffer[rlen-1] = 0;

      // Realloc stacks if necessary.
      if (chrcount >= structsize) {
         structsize *= 2;
         start = realloc(start, structsize*sizeof(long));
         names = realloc(names, structsize*sizeof(char*));
         if (start == NULL || names == NULL) {
            fprintf(stderr, "error in 'read_CHRindex' (realloc): %s\n", strerror(errno));
            exit(EXIT_FAILURE);
         }
      }

      // Save chromosome start position.
      buffer[i] = 0;
      start[chrcount] = atol(buffer);
      // Save chromosome name.
      name = buffer + i + 1;
      names[chrcount] = malloc(strlen(name)+1);
      strcpy(names[chrcount], name);
      // Inc.
      chrcount++;
   }

   // Return chromosome index structure.
   chr_t * chrindex = malloc(sizeof(chr_t));
   chrindex->nchr = chrcount;
   chrindex->start= start;
   chrindex->name = names;

   // Close files.
   fclose(input);

   // Free memory.
   free(buffer);

   return chrindex;
}


seqstack_t *
read_file
(
   FILE      * inputf
)
{
   // TODO: Set a maximum read size in MB, when this limit is reached return.
   // Once the read sequences are processed, read again and process, continuing
   // at the point of the file where the last 'read_file' returned.

   // Read first line of the file to guess format.
   // Store in global variable FORMAT.
   char c = fgetc(inputf);

   if (c == EOF) return NULL;
   if (ungetc(c, inputf) == EOF) {
      return NULL;
   }

   format_t format;

   switch(c) {
   case '>':
      format = FASTA;
      break;
   case '@':
      format = FASTQ;
      break;
   default:
      format = RAW;
   }

   ssize_t nread;
   size_t linesz = BUFFER_SIZE;
   size_t sublinesz = BUFFER_SIZE;
   size_t tempsz = BUFFER_SIZE;
   char *line = malloc(BUFFER_SIZE * sizeof(char));
   char *subline = malloc(BUFFER_SIZE * sizeof(char));
   char *temp = malloc(BUFFER_SIZE * sizeof(char));
   if (line == NULL) {
      return NULL;
   }
   seqstack_t *seqstack = new_seqstack(SEQSTACK_SIZE);
   if (seqstack == NULL) {
      return NULL;
   }

   char *seq  = NULL;
   char *tag  = NULL;
   char *q    = NULL;
   int lineno = 0;

   while ((nread = getline(&line, &linesz, inputf)) != -1) {
      lineno++;
      if (line[nread-1] == '\n') line[nread-1] = '\0';
      switch(format) {
      case RAW:
         seq = line;
         tag = line;
         break;
      case FASTA:
         if (line[0] == '>') {
            if((nread = getline(&subline, &sublinesz, inputf)) == -1) {
               fprintf(stderr, "incorrect FASTA format: line %d\n", lineno);
               continue;
            }
            lineno++;
            if (subline[nread-1] == '\n') subline[nread-1] = '\0';
            tag = line;
            seq = subline;
         } else continue;
         break;
      case FASTQ:
         if (line[0] == '@') {
            if((nread = getline(&subline, &sublinesz, inputf)) == -1) {
               fprintf(stderr, "incorrect FASTQ format: line %d\n", lineno);
               continue;
            }
            lineno++;
            if (subline[nread-1] == '\n') subline[nread-1] = '\0';
            tag = line;
            seq = subline;
            for (int i = 0; i < 2; i++) {
               if((nread = getline(&temp, &tempsz, inputf)) == -1) {
                  fprintf(stderr, "incorrect FASTQ format: line %d\n", lineno);
                  continue;
               }
               lineno++;
            }
            q = temp;
         } else continue;
         break;
      }

      /*size_t seqlen = strlen(seq);
      if (seqlen > MAXSEQLEN) {
         fprintf(stderr, "max sequence length exceeded (%d)\n", MAXSEQLEN);
         fprintf(stderr, "offending sequence:\n%s\n", seq);
         continue;
         }*/
      
      
      if (seq_push(&seqstack, tag, seq, q)) return NULL;
   }

   free(line);
   free(subline);
   free(temp);
   return seqstack;
}
