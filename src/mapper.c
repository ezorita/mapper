#include "mapper.h"
#include <assert.h>
#include <signal.h>
#include <execinfo.h>
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

int main(int argc, char *argv[])
{
   // Redirect signals.
   signal(SIGSEGV, SIGSEGV_handler);
   // Force verbose.
   int opt_verbose = 1;

   // Without params, print help.
   if (argc == 1) {
      say_map_usage();
      return EXIT_SUCCESS;
   }

   if (strcmp(argv[1],"index") == 0) {
      // Without params, print help.    
      if (argc == 2 || (argc == 3 && !strcmp(argv[2],"-h"))) {
         say_index_usage();
         return EXIT_SUCCESS;
      } else if (strcmp(argv[2], "add") == 0) {
         if (argc == 3) {
            say_add_usage();
            return EXIT_SUCCESS;
         }
         // Parse params.
         char * index_file;
         opt_add_t opt;
         int rval;
         // TODO: 
         // - Remove seed hash table params.
         if ((rval = parse_opt_add(argc, argv, &opt, &index_file))) {
            return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
         }
         // Load index.
         index_t * index = index_load_base(index_file);
         if (index == NULL)
            return EXIT_FAILURE;
         // Compute annotation.
         return index_add_annotation(opt.k, opt.d, opt.threads, index, index_file);
      } else if (strcmp(argv[2], "build") == 0) {
         if (argc == 3) {
            say_build_usage();
            return EXIT_SUCCESS;
         }
         // Parse params.
         char * genome_file;
         opt_add_t opt;
         int rval;
         if ((rval = parse_opt_build(argc, argv, &opt, &genome_file))) {
            return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
         }

         return write_index(genome_file, opt.k, opt.d, opt.threads);
      }
   } else if (strcmp(argv[1],"neigh") == 0) {
      if (argc != 5) {
         fprintf(stderr,"Find neighbors of a sequence:\n\n  mapper neigh [index] [dist] [sequence]\n");
         return EXIT_FAILURE;
      }
      // Load index.
      index_t * index = index_load_base(argv[2]);
      if (index == NULL) {
         fprintf(stderr, "error loading index.\n");
         exit(1);
      }
      int slen = strlen(argv[4]);
      // Sequence and exact repeats.
      fprintf(stdout, "%s in %s:\n", argv[4], argv[2]);
      bwpos_t repeat_pos;
      if(suffix_string(argv[4], slen, 0, &repeat_pos, index->bwt))
         fprintf(stdout, "[!] Sequence not found!\n");
      else
         fprintf(stdout, "Exact repeats: %ld\n", repeat_pos.ep - repeat_pos.sp + 1);
         
      // 1- neighbors.
      seedstack_t * stack = seedstack_new(10);
      seq_neighbors(argv[4], slen, 0, atoi(argv[3]), &stack, index->bwt->bwt_base, index);
      if (stack->pos == 0) {
         fprintf(stdout, "[!] No %d-neighbors found!\n", atoi(argv[3]));
         return EXIT_SUCCESS;
      }
      fprintf(stdout, "%d-neighbors of %s in %s:\n", atoi(argv[3]), argv[4], argv[2]);
      // Alloc buffer for neigh sequence.
      char * nseq = malloc(slen+1);
      nseq[slen] = 0;
      for (int i = 0; i < stack->pos; i++) {
         seed_t s = stack->seed[i];
         uint64_t sa = get_sa(s.ref_pos.sp, index->sar);
         memcpy(nseq, index->genome + sa, slen);
         fprintf(stdout, "%s\td=%d\tn=%ld\tsp=%ld\tsar[0]=%ld\n", nseq, s.errors, s.ref_pos.ep - s.ref_pos.sp + 1, s.ref_pos.sp, sa);
      }
      free(nseq);
   }

   // DEFAULT Options.
   /* ILLUMINA
   seedopt_t seedopt = { // MEM
      .min_len = 19,
      .max_len = 1000,
      .min_loci = 0,
      .max_loci = 50,
      .thr_loci = 20,
      .reseed_len = 28,
      .sensitive_mem = 0
   };
   */

   /*
   // NANOPORE 2D
   seedopt_t seedopt = { // MEM
      .min_len = 14,
      .max_len = 1000,
      .min_loci = 0,
      .max_loci = 500,
      .thr_loci = 100,
      .reseed_len = 10,
      .sensitive_mem = 0
   };

   filteropt_t filteropt = {
      .target_q = 20,
      .dist_accept = 10,
      .max_align_per_read = 200,
      .read_ref_ratio = 0.15,
      .align_accept_eexp = -100.0,
      .overlap_max_tolerance = 0.5,
      .align_seed_filter_thr = 0.5,
      .align_seed_filter_dif = 19,//38,
//      .align_filter_ident = 0.85,
      .align_filter_ident = 0.0,
      .align_filter_eexp = 0.0,
      .mapq_evalue_ratio = 0.5,
      .min_best_to_second = 0.75,
      .min_interval_size = 25,
      .min_interval_score = 30,
      .min_interval_hits = 30,
      .min_gap_id = 0.8
   };
   alignopt_t alignopt = {
      .bp_diagonal      = 0,
      .bp_thr           = 10,
      .bp_max_thr       = 50,
      .bp_resolution    = 1,
      .bp_period        = 5,
      .bp_repeats       = 4,
      .read_error       = 0.15,
      .rand_error       = 0.50,
      .width_ratio      = 0.15,
      .mismatch_penalty = 1
   };

   // Precompute logarithms.
   alignopt.logAe = log(alignopt.rand_error) - log(alignopt.read_error);
   alignopt.logBe = log(1-alignopt.rand_error) - log(1-alignopt.read_error);

   formatopt_t formatopt = {
      .print_first = 2, // 2 for first_only
      .mapq_thr = 0,
      .eval_thr = 0
   };
   param_map_t mapopt = {
      .seed = seedopt,
      .align = alignopt,
      .filter = filteropt,
      .format = formatopt,
      .threads = 1
   };

   // Parse input parameters.
   if (argc == 1) {
      say_map_usage();
      exit(EXIT_SUCCESS);
   }
   opt_map_t in_params = {0,0,0,0};
   char * qfile, *ifile;
   parse_opt_map(argc, argv, &in_params, &qfile, &ifile);
   if (in_params.threads) mapopt.threads = in_params.threads;
   

   // Parse query filename.
   FILE * queryfile = fopen(qfile, "r");
   if (queryfile == NULL) {
      fprintf(stderr, "error: could not open file: %s (%s).\n", argv[1],strerror(errno));
      return EXIT_FAILURE;
   }

   // Read index.
   if (opt_verbose) fprintf(stderr, "loading index...\n");
   index_t * index = index_load_base(ifile);
   if (index == NULL)
      return EXIT_FAILURE;

   // Read annotations and seed tables.
   annlist_t * ann = ann_index_read(ifile);
   if (ann == NULL)
      return EXIT_FAILURE;
   if (ann->count == 0) {
      fprintf(stderr, "[error] no annotations found.\n");
      return EXIT_FAILURE;
   }
   if (ann_index_load(ann, ifile))
      return EXIT_FAILURE;
   ann_print_index(ann);

   index->ann = ann;

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
   free(index);
   */
   return EXIT_SUCCESS;
}

int
mt_scheduler
(
 seqstack_t * seqs,
 param_map_t* options,
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
/*
   mtjob_t * job = (mtjob_t *) args;
   // Get params.
   seq_t       * seq   = job->seq;
   param_map_t * opt   = job->opt;
   index_t     * index = job->index;

   // Alloc data structures.
   matchlist_t * map_matches = matchlist_new(64);
   seedstack_t * seeds    = seedstack_new(SEEDSTACK_SIZE);
   seedstack_t * reseeds  = seedstack_new(SEEDSTACK_SIZE);
   seedstack_t * mems     = seedstack_new(SEEDSTACK_SIZE);
   vstack_t    * mem_intv = new_stack(SEEDSTACK_SIZE);
   pathstack_t * pstack   = pathstack_new(PATHSTACK_DEF_SIZE);
   //   int k = index->sht->sht[0].k;
   //   int d = index->sht->sht[0].d;
   //   int r = index->sht->sht[0].repeat_thr;

   // Parametrize
   int max_align = 20;
   int target_q = opt->filter.target_q;

   // Precompute log table.
   int cumlog_max = 0;
   double * cumlog = malloc(sizeof(double));
   *cumlog = 0;
   

   size_t mapped = 0;
   
   for (size_t i = 0; i < job->count; i++) {
      if (VERBOSE_DEBUG) fprintf(stdout, "%s\t%s\n", seq[i].tag, seq[i].seq);

      // Read length and alloc memory.
      size_t slen = strlen(seq[i].seq);
      uint8_t * query = malloc(slen);
      uint8_t * value = malloc(slen);
      for (int n = 0; n < slen; n++) value[n] = 255;

      // Extend cumlog table.
      if (slen > cumlog_max) {
         cumlog = realloc(cumlog, (slen+1)*sizeof(double));
         for (int c = cumlog_max+1; c <= slen; c++)
            cumlog[c] = cumlog[c-1] + log10(c);
         cumlog_max = slen;
      }

      // Translate query.
      for (int j = 0; j < slen; j++) query[j] = translate[(int)seq[i].seq[j]];
      // Reset matches.
      map_matches->pos = 0;
*/
      /**
      *** UNIQUE SEEDS
      **/
/*
      while (1) {
         // Do not seed intervals.
         if (map_matches->pos) {
            for (int j = 0; j < map_matches->pos; j++) {
               match_t m = map_matches->match[j];
               if (beg > m.read_e) continue;
               if (beg < m.read_s) break;
               beg = m.read_e + 1;
               break;
            }
         }
         // Find first unique seed.
         hit_t hit = find_uniq_seed(beg,slen,query,value,index);
         if (hit.errors) {
            // No seed found.
            break;
         } else {
            align_hits(seq[i].seq, &hit, 1, &map_matches, index, opt->filter, opt->align);
         }
         beg = hit.qrypos + 1;
      }

      // DEBUG.
      if (VERBOSE_DEBUG) {
         fprintf(stdout, "** path 1 **\n");
         for (int k = 0; k < map_matches->pos; k++) {
            match_t m = map_matches->match[k];
            fprintf(stdout, "match[%d]: read_s:%d, read_e:%d, ref_s:%ld, hits:%d, s_hits:%d (%d), ann:%d\n", k, m.read_s, m.read_e, m.ref_s, m.hits, m.s_hits, m.s_hits_cnt, m.ann_d);
         }
      }

      if (map_matches->pos) goto map_end;
*/

      /**
      *** UNIQUE MEMS
      **/

/* TO COMPILE
      if (VERBOSE_DEBUG) fprintf(stdout, "MEMS:\n");
      hit_t hit;
      int beg = 0;
      mems->pos = 0;
      seeds->pos = 0;
      reseeds->pos = 0;
      while (beg < slen) {
         seed_t mem;
         beg = next_mem(query,beg,slen-1,&mem,index);
         // Select mems longer than min_len.
         if (mem.ref_pos.depth < opt->seed.min_len) continue;
         // Check neighborhood.
         if (mem_unique(mem, index, &hit)) {
            if (VERBOSE_DEBUG) {
               fprintf(stdout, "unique mem: beg:%d, end:%d, loci:%ld\n", mem.qry_pos, mem.qry_pos+mem.ref_pos.depth, mem.ref_pos.ep - mem.ref_pos.sp + 1);
            }
            align_hits(seq[i].seq, &hit, 1, &map_matches, cumlog, index, opt->filter, opt->align);
         } else {
            if (VERBOSE_DEBUG) {
               fprintf(stdout, "beg:%d, end:%d, loci:%ld\n", mem.qry_pos, mem.qry_pos+mem.ref_pos.depth, mem.ref_pos.ep - mem.ref_pos.sp + 1);
            }
            if (mem.ref_pos.ep - mem.ref_pos.sp > opt->seed.max_loci)
               // Repeat seeds. Do not reseed.
               mem.ref_pos.ep = mem.ref_pos.sp + opt->seed.max_loci - 1;
            else if (mem.ref_pos.depth >= opt->seed.reseed_len)
               // Reseed mems longer that reseed_len.
               reseed_mem(query, mem, slen, opt->seed.min_len, index, &reseeds);
            // Store MEM.
            seedstack_push(mem, &mems);
         }
         // Update beg.
         for (int j = 0; j < map_matches->pos; j++) {
            match_t m = map_matches->match[j];
            if (beg > m.read_e) continue;
            if (beg < m.read_s - opt->seed.min_len) break;
            if ((m.mapq >= target_q || m.maxq < target_q) &&  m.gap_id >= opt->filter.min_gap_id)
               beg = m.read_e + 1;
            if (VERBOSE_DEBUG) {
               fprintf(stdout, "match[%d]: r_beg:%d, r_end:%d, hits:%d, a:%d, score:%d, mapQ=%d, maxQ=%d, gap_id:%.2f, next_beg:%d\n", j, m.read_s, m.read_e, m.hits, m.ann_d, m.score, m.mapq, m.maxq, m.gap_id, beg);
            }
         }
      }


      if (VERBOSE_DEBUG) {
         // Print reseeds.
         fprintf(stdout, "RESEEDS:\n");
         for (int j = 0; j < reseeds->pos; j++) {
            seed_t mem = reseeds->seed[j];
            fprintf(stdout, "beg:%d, end:%d, loci:%ld\n", mem.qry_pos, mem.qry_pos+mem.ref_pos.depth, mem.ref_pos.ep - mem.ref_pos.sp + 1);
         }

         // Print intervals.
         for (int j = 0; j < map_matches->pos; j++) {
            match_t m = map_matches->match[j];
            fprintf(stdout, "match[%d]: r_beg:%d, r_end:%d, hits:%d, a:%d, score:%d, mapQ=%d, maxQ=%d, gap_id:%.2f, next_beg:%d\n", j, m.read_s, m.read_e, m.hits, m.ann_d, m.score, m.mapq, m.maxq, m.gap_id, beg);
         }
      }
*/
      /**
      *** MEM/THRESHOLD ALIGNMENT
      **/

/* TO COMPILE
      // Recompute unmapped intervals.
      mem_intv->pos = 0;
      mem_intervals(map_matches, slen, target_q,opt->filter.min_gap_id, opt->filter.min_interval_size, &mem_intv);
      if (mem_intv->pos == 0) goto map_end;

      if (VERBOSE_DEBUG) {
         fprintf(stdout,"unmapped intervals:\n");
         for (int j = 0; j < mem_intv->pos; j+=2) {
            fprintf(stdout, "beg: %ld, end:%ld\n", mem_intv->val[j], mem_intv->val[j+1]);
         }
      }

      // Append reseeds.
      for (int k = 0; k < reseeds->pos; k++) {
         seed_t s = reseeds->seed[k];
         if (s.ref_pos.ep - s.ref_pos.sp > opt->seed.max_loci)
            s.ref_pos.ep = s.ref_pos.sp + opt->seed.max_loci - 1;
         seedstack_push(s, &mems);
      }

      // Remove mapped mems.
      seeds->pos = 0;
      for (int j = 0; j < mems->pos; j++) {
         seed_t m = mems->seed[j];
         for (int s = 0; s < mem_intv->pos; s+= 2) {
            if (mem_intv->val[s+1] < m.qry_pos) continue;
            if (mem_intv->val[s] >= m.qry_pos + m.ref_pos.depth) break;
            seedstack_push(m, &seeds);
            break;
         }
      }

      // Perform threshold seeding on unmapped intervals.
      for (int s = 0; s < mem_intv->pos; s+= 2) {
         int b = mem_intv->val[s];
         int e = mem_intv->val[s+1];
         int nsd = seeds->pos;
         seed_interv(query, b, e, opt->seed.min_len+1, opt->seed.thr_loci, index, &seeds);
         if (VERBOSE_DEBUG) {
            fprintf(stdout, "LAST:\n");
            for (int j = nsd; j < seeds->pos; j++) {
               seed_t s = seeds->seed[j];
               fprintf(stdout, "beg: %d, end: %d, loci: %ld\n", s.qry_pos, s.qry_pos + s.ref_pos.depth, s.ref_pos.ep - s.ref_pos.sp + 1);
            }
         }
      }

      // Remove duplicates.
      remove_dup_mem(seeds);

      // Chain seeds and align.
      if (seeds->pos) {
         size_t hit_cnt = 0;
         hit_t * hits = chain_seeds(seeds, slen, index, opt->filter.dist_accept, opt->filter.read_ref_ratio, &hit_cnt);
         //hit_t * hits = compute_hits(seeds, index, &hit_cnt);
         hit_cnt = min(opt->filter.max_align_per_read, hit_cnt);
         for (int j = 0; j < hit_cnt; j++)
            hits[j].s_errors = 0;
         // Sort and align up to max_align.
         if (hit_cnt) {
            align_hits(seq[i].seq, hits, hit_cnt, &map_matches, cumlog, index, opt->filter, opt->align);
         }
         free(hits);
      }
*/
      /**
      *** MISMATCHED MEMS
      **/

/* TO COMPILE
      // Recompute seed intervals.
      mem_intv->pos = 0;
      mem_intervals(map_matches, slen, target_q,opt->filter.min_gap_id, opt->filter.min_interval_size, &mem_intv);
      if (mem_intv->pos == 0) goto map_end;

      // Remove mapped mems.
      seeds->pos = 0;
      for (int j = 0; j < mems->pos; j++) {
         seed_t m = mems->seed[j];
         for (int s = 0; s < mem_intv->pos; s+= 2) {
            if (mem_intv->val[s+1] < m.qry_pos) continue;
            if (mem_intv->val[s] >= m.qry_pos + m.ref_pos.depth) break;
            seedstack_push(m, &seeds);
            break;
         }
      }

      // Compute 1-mismatch sequences.
      if (VERBOSE_DEBUG) fprintf(stdout, "Mismatch MEMS:\n");
      mems->pos = 0;
      for (int j = 0; j < seeds->pos; j++) {
         // Reset stacks.
         pstack->pos = 0;
         // Perform 1-mismatch (every k nt?) block search on seed.
         seed_t s = seeds->seed[j];
         blocksearch(query + s.qry_pos, s.ref_pos.depth, 1, index, &pstack);
         int64_t loci = 0;
         // Count loci and remove perfect seed.
         for (int h = 0; h < pstack->pos; h++) {
            // Remove perfect seed.
            if (pstack->path[h].score == 0)
               pstack->path[h] = pstack->path[--pstack->pos];
            // Count loci.
            loci += pstack->path[h].pos.sz;
         }
         if (VERBOSE_DEBUG) fprintf(stdout, "beg:%d, end:%d, seqs:%ld, loci:%ld\n", s.qry_pos, s.qry_pos + s.ref_pos.depth - 1, pstack->pos, loci);         
         // Continue if repeat_thr is exceeded.
         //         if (loci > r) {
         //            continue;
         //         } else if (loci > max_align) {
         if (loci > max_align) {
            // Sort by score, then loci.
            qsort(pstack->path, pstack->pos, sizeof(spath_t), compar_path_score);
            loci = 0;
            for (int h = 0; h < pstack->pos && loci < max_align; h++) {
               spath_t path = pstack->path[h];
               fmdpos_t fmdp = path.pos;
               loci += fmdp.sz;
               if (loci > max_align) fmdp.sz -= loci-max_align;
               bwpos_t pos = {.depth = fmdp.dp, .sp = fmdp.fp, .ep = fmdp.fp + fmdp.sz - 1};
               seedstack_push((seed_t) {.errors = path.score, .qry_pos = beg, .ref_pos = pos}, &mems);
            }
         } else {
            for (int h = 0; h < pstack->pos; h++) {
               spath_t path = pstack->path[h];
               fmdpos_t fmdp = path.pos;
               bwpos_t pos = {.depth = fmdp.dp, .sp = fmdp.fp, .ep = fmdp.fp + fmdp.sz - 1};
               seedstack_push((seed_t) {.errors = path.score, .qry_pos = beg, .ref_pos = pos}, &mems);
            }
         }
      }

      // Align 1-mismatch seeds.
      if (seeds->pos)
         align_seeds(seq[i].seq, seeds, &map_matches, cumlog, index, opt->filter, opt->align);

      goto map_end;

*/
      /**
      *** MEM SEEDING ON GAPS
      **/
/*
      seeds->pos = 0;
      for (int s = 0; s < mem_intv->pos; s += 2) {
         int beg = mem_intv->val[s];
         int end = mem_intv->val[s+1];
         mems->pos = 0;
         seed_mem(query, beg, end, index, &mems);
         // SMEM and SMEM reseeding.
         // TODO:
         // - Add repeat control and threshold number of alingments per seed.
         if (VERBOSE_DEBUG) fprintf(stdout, "SMEMS:\n");
         seeds->pos = 0;
         for (int j = 0; j < mems->pos; j++) {
            seed_t s = mems->seed[j];
            if (VERBOSE_DEBUG) fprintf(stdout, "beg: %d, end: %d, loci: %ld\n", s.qry_pos, s.qry_pos + s.ref_pos.depth, s.ref_pos.ep - s.ref_pos.sp + 1);
            // Simple thresholding. Will do this only on super-repeated seeds,
            // because the seed content will be used to chain in MEM mode...
            if (s.ref_pos.ep - s.ref_pos.sp > opt->seed.max_loci) {
               s.ref_pos.ep = s.ref_pos.sp + opt->seed.max_loci - 1;
               seedstack_push(s, &seeds);
               continue;
            }
            seedstack_push(s, &seeds);
            if (s.ref_pos.depth >= opt->seed.reseed_len) {
               reseeds->pos = 0;
               reseed_mem(query, s, slen, opt->seed.min_len, index, &reseeds);
               for (int k = 0; k < reseeds->pos; k++) {
                     seed_t s = reseeds->seed[k];
                     if (VERBOSE_DEBUG) fprintf(stdout, "'-> beg: %d, end: %d, loci: %ld\n", s.qry_pos, s.qry_pos + s.ref_pos.depth, s.ref_pos.ep - s.ref_pos.sp + 1);
                     if (s.ref_pos.ep - s.ref_pos.sp > opt->seed.max_loci)
                        s.ref_pos.ep = s.ref_pos.sp + opt->seed.max_loci - 1;
                     seedstack_push(s, &seeds);
               }
            }
         }
         
        int nsd = seeds->pos;

         // Threshold seeding.
         seed_interv(query, beg, end, opt->seed.min_len+1, opt->seed.thr_loci, index, &seeds);
         if (VERBOSE_DEBUG) {
            fprintf(stdout, "LAST:\n");
            for (int j = nsd; j < seeds->pos; j++) {
               seed_t s = seeds->seed[j];
               fprintf(stdout, "beg: %d, end: %d, loci: %ld\n", s.qry_pos, s.qry_pos + s.ref_pos.depth, s.ref_pos.ep - s.ref_pos.sp + 1);
            }
         }

      }

      // Chain seeds.
      if (seeds->pos) {
         size_t hit_cnt = 0;
         hit_t * hits = chain_seeds(seeds, slen, index, opt->filter.dist_accept, opt->filter.read_ref_ratio, &hit_cnt);
         hit_cnt = min(opt->filter.max_align_per_read, hit_cnt);
         for (int j = 0; j < hit_cnt; j++)
            hits[j].errors = 0;
         // Sort and align up to max_align.
         if (hit_cnt) {
            align_hits(seq[i].seq, hits, hit_cnt, &map_matches, index, opt->filter, opt->align);
         }
         free(hits);
      }

      if (opt->seed.sensitive_mem) {
         for (int s = 0; s < mem_intv->pos; s += 2) {
            int beg = mem_intv->val[s];
            int end = mem_intv->val[s+1];
            // Find MEMs with more than 2 hits.
            mems->pos = 0;
            seed_mem_bp(query, beg, end, 2, index, &mems);
            for (int j = 0; j < mems->pos; j++) {
               seed_t s = mems->seed[j];
               uint64_t loci = s.ref_pos.ep - s.ref_pos.sp + 1;
               // Save MEMS with loci <= r.
               if (loci > r) {
                  s.ref_pos.ep = s.ref_pos.sp + r - 1;
               }
               seedstack_push(s, &seeds);
            }
         }
      } else {
         seedstack_t * resd = seedstack_new(64);
         for (int s = 0; s < mem_intv->pos; s += 2) {
            int beg = mem_intv->val[s];
            int end = mem_intv->val[s+1];
            if (VERBOSE_DEBUG) fprintf(stdout, "MEM interval: beg:%d, end:%d\n", beg,end);
            // Find MEMS.
            mems->pos = 0;
            //            seed_mem(query, beg, end, index, &seeds);
            //            get_smem(&mems, seeds);
            seed_thr(query, beg, end, max_align, index, &mems);
            seeds->pos = 0;
            // DEBUG.
            for (int j = 0; j < mems->pos; j++) {
               seed_t s = mems->seed[j];
               if (s.ref_pos.depth < opt->seed.min_len) continue;
               uint64_t loci = s.ref_pos.ep - s.ref_pos.sp + 1;
               if (VERBOSE_DEBUG) fprintf(stdout, "MEM (%d/%ld): beg:%d, end:%d, loci:%ld.\n",j+1,mems->pos,s.qry_pos, s.qry_pos+s.ref_pos.depth-1, loci);
               if (loci > max_align) {
                  // Do maximum 'max_align' alignments per seed.
                  s.ref_pos.ep = s.ref_pos.sp + max_align - 1;
               }
               else if (loci == 1 && s.ref_pos.depth >= opt->seed.reseed_len) {
                  // Reseed.
                  seed_mem_bp(query, s.qry_pos, s.qry_pos + s.ref_pos.depth, loci + 1, index, &resd);
               }
               seedstack_push(s, &seeds);
            }
            mems->pos = 0;
            get_smem(&mems, resd);
            for (int j = 0; j < mems->pos; j++) {
               seed_t s = mems->seed[j];
               uint64_t loci = s.ref_pos.ep - s.ref_pos.sp + 1;
               if (VERBOSE_DEBUG) fprintf(stdout, "RESEED-MEM (%d/%ld): beg:%d, end:%d, loci:%ld.\n",j+1,resd->pos,s.qry_pos, s.qry_pos+s.ref_pos.depth-1, loci);
               if (loci > max_align) {
                  // Do maximum 'max_align' alignments per seed.
                  s.ref_pos.ep = s.ref_pos.sp + max_align - 1;
               }
               seedstack_push(s, &seeds);
            }
         }
         free(resd);
      }
*/

      /**
      *** NON-UNIQUE AND NOT FOUND
      **/

      // TODO:
      // - Adaptar aquesta part del codi per fer-lo funcionar amb MEMS.
      // - Aqui han quedat alguns espais buits dins el read o alguns alineaments
      //   que no tenen prou suport. El que fa falta fer es definir una estrategia
      //   per fer mismatched seeding en aquests intervals irregulars.

      // ### De moment ho deixo comentat pero aquesta estrategia va al soft final ###
      /*
      int v[2] = {2,0};
      // Align acceptable seeds (less than repeat_thr loci) then the not found ones.
      for (int l = 0; l < 2; l++) {
         // DEBUG.
         if (VERBOSE_DEBUG) fprintf(stdout, "** path %d **\n", v[l]);
         int beg = 0;
         while (beg <= (int)(slen - opt->seed.min_len)) {
            // Do not seed path 1 intervals or alignments with second best score.
            int next = 0;
            for (int j = 0; j < map_matches->pos; j++) {
               match_t m = map_matches->match[j];
               if (beg > m.read_e) continue;
               if (beg < m.read_s - opt->seed.min_len) break;
               if (m.s_hits_cnt > 0) {
                  beg = m.read_e + 1;
                  next = 1;
                  break;
               }
            }
            if (next) continue;
            // Check seed type            
            if (value[beg] == v[l]) {
               // Compute (k,d) matches.
               pstack->pos = 0;
               seeds->pos = 0;
               blocksearch(query + beg, k, d, index, &pstack);
               int64_t loci = 0;
               for (int h = 0; h < pstack->pos; h++)
                  loci += pstack->path[h].pos.sz;
               // Continue if repeat_thr is exceeded.
               if (loci > r) {
                  // DEBUG.
                  if (VERBOSE_DEBUG) fprintf(stdout, "[beg:%d] too many loci (%ld)\n", beg, loci);
                  beg++;
                  continue;
               } else if (loci > max_align) {
                  // Sort to prioritize low scores
                  qsort(pstack->path, pstack->pos, sizeof(spath_t), compar_path_score);
                  loci = 0;
                  for (int h = 0; h < pstack->pos && loci < max_align; h++) {
                     spath_t path = pstack->path[h];
                     fmdpos_t fmdp = path.pos;
                     loci += fmdp.sz;
                     if (loci > max_align) fmdp.sz -= loci-max_align;
                     bwpos_t pos = {.depth = fmdp.dp, .sp = fmdp.fp, .ep = fmdp.fp + fmdp.sz - 1};
                     seedstack_push((seed_t) {.errors = path.score, .qry_pos = beg, .ref_pos = pos}, &seeds);
                  }
               } else {
                  for (int h = 0; h < pstack->pos; h++) {
                     spath_t path = pstack->path[h];
                     fmdpos_t fmdp = path.pos;
                     bwpos_t pos = {.depth = fmdp.dp, .sp = fmdp.fp, .ep = fmdp.fp + fmdp.sz - 1};
                     seedstack_push((seed_t) {.errors = path.score, .qry_pos = beg, .ref_pos = pos}, &seeds);
                  }
               }
               if (seeds->pos) {
                  // DEBUG.
                  if (VERBOSE_DEBUG) fprintf(stdout, "[beg:%d] aligning %ld seeds\n", beg, seeds->pos);
                  // Align non-unique seeds.
                  align_seeds(seq[i].seq, seeds, &map_matches, cumlog, index, opt->filter, opt->align);
                  // DEBUG.
                  if (VERBOSE_DEBUG) {
                     for (int k = 0; k < map_matches->pos; k++) {
                        match_t m = map_matches->match[k];
                        fprintf(stdout, "match[%d]: read_s:%d, read_e:%d, ref_s:%ld, hits:%d, s_hits:%d (%d), ann:%d\n", k, m.read_s, m.read_e, m.ref_s, m.hits, m.s_hits, m.s_hits_cnt, m.ann_d);
                     }
                  }
               }
            } else if (value[beg] == 255) {
               // Check seed in seed table.
               hit_t hit = find_uniq_seed(beg,beg + k,query,value,index);
               if (hit.s_errors == 0) {
                  align_hits(seq[i].seq, &hit, 1, &map_matches, cumlog, index, opt->filter, opt->align);
               }
               continue;
            } 
            beg++;
         }
      }
      */
/* TO COMPILE
      // DEBUG.
      if (VERBOSE_DEBUG) {
         for (int k = 0; k < map_matches->pos; k++) {
            match_t m = map_matches->match[k];
            fprintf(stdout, "match[%d]: read_s:%d, read_e:%d, ref_s:%ld, hits:%d, s_hits:%d (%d), ann:%d\n", k, m.read_s, m.read_e, m.ref_s, m.hits, m.s_hits, m.s_hits_cnt, m.ann_d);
         }
      }
   map_end:
      // Score and output.
      //      map_score(map_matches,cumlog);
      mapped += print_and_free(seq[i], map_matches, index, opt->format);
      free(query);
      free(value);

   }

   // Free data structures.
   free(seeds);
   free(reseeds);
   free(mems);
   free(map_matches);
   
   // Report to scheduler.
   pthread_mutex_lock(job->mutex);
   job->control->mapped += mapped;
   job->control->active -= 1;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free args.
   free(args);
*/
   return NULL;
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
