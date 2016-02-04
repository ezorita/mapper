#include "mapper.h"
#include <assert.h>
#include <time.h>



int main(int argc, char *argv[])
{
   int opt_verbose = 1;
   size_t seqs_per_thread = 10000;

   // Without params, print help.
   if (argc == 1) {
      say_map_usage();
      return EXIT_SUCCESS;
   }

   if (strcmp(argv[1],"index") == 0) {
      // Without params, print help.    
      if (argc == 2) {
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
         if ((rval = parse_opt_add(argc, argv, &opt, &index_file))) {
            return (rval == -1 ? EXIT_FAILURE : EXIT_SUCCESS);
         }
         // Load index.
         index_t * index = index_load_base(index_file);
         if (index == NULL)
            return EXIT_FAILURE;
         // Compute annotation.
         return index_add_annotation(opt.k, opt.d, opt.repeat_thr, opt.mode,
                                     opt.threads, index, index_file);
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

         return write_index(genome_file, opt.k, opt.d, opt.repeat_thr, opt.threads);
      }
   }

   // DEFAULT Options.
   seedopt_t seedopt = { // MEM
      .min_len = 12,
      .max_len = 1000,
      .min_loci = 0,
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
//      .align_filter_ident = 0.85,
      .align_filter_ident = 0.0,
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

   shtlist_t * sht = sht_index_read(ifile);
   if (sht == NULL)
      return EXIT_FAILURE;
   if (sht->count == 0) {
      fprintf(stderr, "[error] no seed tables found.\n");
      return EXIT_FAILURE;
   }
   if (sht_index_load(sht, ifile))
      return EXIT_FAILURE;
   sht_print_index(sht);

   index->ann = ann;
   index->sht = sht;

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
   mtjob_t * job = (mtjob_t *) args;
   // Get params.
   seq_t       * seq   = job->seq;
   param_map_t * opt   = job->opt;
   index_t     * index = job->index;

   // Alloc data structures.
   matchlist_t * seed_matches = matchlist_new(opt->filter.max_align_per_read);
   matchlist_t * map_matches = matchlist_new(64);
   seedstack_t * seeds    = seedstack_new(SEEDSTACK_SIZE);
   seedstack_t * thrseeds = seedstack_new(SEEDSTACK_SIZE);
   seedstack_t * reeds    = seedstack_new(SEEDSTACK_SIZE);
   pathstack_t * pstack   = pathstack_new(PATHSTACK_DEF_SIZE);

   size_t mapped = 0;
   
   for (size_t i = 0; i < job->count; i++) {
      if (VERBOSE_DEBUG) fprintf(stdout, "seq: %s\n", seq[i].seq);
      size_t slen = strlen(seq[i].seq);
      // Reset seeds.
      seeds->pos = 0;
      reeds->pos = 0;
      thrseeds->pos = 0;

      // Translate query.
      uint8_t * query = malloc(slen);
      for (int j = 0; j < slen; j++) query[j] = translate[(int)seq[i].seq[j]];

      // Find unique seeds.
      //find_uniq_seeds(slen,25,query,index,&seeds,&thrseeds);
      /**/
      int zero_cnt = 0;
      seed_t sd = find_uniq_seed(0,slen,25,query,index,&zero_cnt);
      if (!sd.bulk) seedstack_push(sd, &seeds);
      /**/
      // Align.
      map_matches->pos = 0;
      align_seeds(seq[i].seq, seeds, reeds, &map_matches, index, opt->filter, opt->align);
      if (map_matches->pos == 0 || map_matches->match[0].interval != 1) {
         int loci, beg = 0;
         do {
            loci = 0;
            if((beg = find_thr_seed(beg, slen, 25, query, index)) == -1) break;
            pstack->pos = 0;
            blocksearch(query + beg, 25, 1, index, &pstack);
            for (int h = 0; h < pstack->pos; h++)
               loci += pstack->path[h].pos.sz;
            beg += 1;
         } while (loci > 20 || loci == 0);

         if (beg >= 0 && loci <= 20) {
            seeds->pos = 0;
            for (int h = 0; h < pstack->pos; h++) {
               fmdpos_t fmdp = pstack->path[h].pos;
               bwpos_t pos = {.depth = fmdp.dp, .sp = fmdp.fp, .ep = fmdp.fp + fmdp.sz - 1};
               seedstack_push((seed_t) {.bulk = 0, .qry_pos = beg-1, .ref_pos = pos}, &seeds);
            }
            align_seeds(seq[i].seq, seeds, reeds, &map_matches, index, opt->filter, opt->align);
         }
         /**/
      }
      // Reseed?
      if (map_matches->pos == 0 && zero_cnt) {

      } 
      else if (map_matches->match[0].interval == 2 && map_matches->match[0].score > 0) {

      }
      mapped += print_and_free(seq[i], map_matches, index, opt->format);
      free(query);
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
