#include "annotate.h"
#include <string.h>
#include <time.h>

// Annotation array.
// Information structure:
// 8 bits per genomic locus.
//   bit 7 \
//   bit 6 |
//   bit 5 | Num of neighbors, in dec (!bit3) or log2 (bit3).
//   bit 4 /
//   bit 3 - Log2/dec flag for bits 4-7.
//   bit 2 - Neighbor alignment bit.
//   bit 1 \ 
//   bit 0 / Distance of closest neighbor (00: 1, 01: 2, 10: 3, 11: 4).

// Inline functions.
inline void send_progress_to_scheduler (int64_t, int64_t*, annjob_t*);

void
job_ranges_rec
(
 fmdpos_t   pos,
 int        tau,
 int        n_cnt,
 int        depth,
 int        max_depth,
 int      * jobcnt,
 annjob_t * jobs,
 index_t  * index
)
{
   // Too many 'N'.
   if (n_cnt > tau) return;
   // If reached suffix depth, make and store job.
   if (depth == max_depth) {
      annjob_t job = {
         .beg = pos.fw,
         .end = pos.fw + pos.sz
      };
      jobs[(*jobcnt)++] = job;
      return;
   }
   
   // Otherwise continue exploring the suffix trie.
   fmdpos_t newpos[NUM_BASES];
   extend_fw_all(pos, newpos, index->bwt);
   for (int i = 0; i < NUM_BASES; i++)
      job_ranges_rec(newpos[i], tau, n_cnt + (i == UNKNOWN_BASE), depth + 1, max_depth, jobcnt, jobs, index);

   return;
}


annotation_t
annotate
(
 int        kmer,
 int        tau,
 index_t  * index,
 int        threads
)
{
   // Constant parameters.
   int job_to_thread_ratio = 5;
   int           word_size = max(3,tau) + 3;

   // Alloc variables.
   pthread_mutex_t * mutex = malloc(sizeof(pthread_mutex_t));
   pthread_cond_t  * monitor = malloc(sizeof(pthread_cond_t));
   uint64_t computed = 0;
   int32_t  done = 0;

   // Temporary info array.
   // Has 'index->size' elements with the following structure:
   //   bytes 0-1: number of different neighbors (16bit).
   //   byte    2: distance of closest neighbor.
   //   byte  3-?: M bytes encoding the mutated positions to find the closest neighbors.
   //   
   //   M is max(3,tau). If the number of mutations is > M, the alignment is not stored.

   uint8_t * info = calloc(index->size, word_size);

   // Initialization.
   pthread_mutex_init(mutex,NULL);
   pthread_cond_init(monitor,NULL);

   // Divide job ranges.
   // For consistency, the best is to split the jobs in ranges of prefixes.
   // For instance, for 6 threads, it may be a good idea to create one job
   // for each prefix of 2nt, i.e. jobs = [1: 'AA', 2:'AC', 3:'AG', 4:'AT',
   // 5:'AN', ..., 25: 'NN']
   // So increase the prefix depth until the number of prefixes (jobs) is at least
   // K(=5?) times the number of threads.
   int prefix_depth = 1;
   int num_jobs = NUM_BASES;
   while (num_jobs < job_to_thread_ratio*threads) {
      prefix_depth++;
      num_jobs *= NUM_BASES;
   }

   // Create jobs based on prefixes of depth 'prefix_depth'.
   annjob_t * jobs = malloc(num_jobs*sizeof(annjob_t));
   num_jobs = 0;
   // Fill specific job info.
   job_ranges_rec(index->bwt->fmd_base, tau, 0, 0, prefix_depth, &num_jobs, jobs, index);
   // Fill constant job info.
   for (int i = 0; i < num_jobs; i++) {
      jobs[i]->kmer     = (uint32_t) kmer;
      jobs[i]->tau      = (uint16_t) tau;
      jobs[i]->wsize    = (uint16_t) word_size;
      jobs[i]->computed = &computed;
      jobs[i]->done     = &done;
      jobs[i]->mutex    = mutex;
      jobs[i]->monitor  = monitor;
      jobs[i]->index    = index;
      jobs[i]->info     = info;
   }

   clock_t clk = clock();
   // Run threads.
   for (int i = 0; i < min(threads, num_jobs); i++) {
      pthread_t thread;
      pthread_create(&thread,NULL,annotate_mt,jobs+i);
      pthread_detach(thread);
      //      annotate_mt(job[i]);
   }
   // Sleep and wait for thread signals.
   int runcount = threads;
   pthread_mutex_lock(mutex);
   while (done < num_jobs) {
      // Run new threads.
      for (int i = runcount; i < min(done + threads, num_jobs); i++) {
         pthread_t thread;
         pthread_create(&thread,NULL,annotate_mt,jobs+i);
         pthread_detach(thread);
      }
      // Update submitted jobs.
      runcount = done + threads;
      // Report progress.
      fprintf(stderr, "[proc] annotating... %.2f%%\r", computed*100.0/index->bwt->fmd_base.sz);
      // Wait for signal.
      pthread_cond_wait(monitor,mutex);
   }
   pthread_mutex_unlock(mutex);

   //   fprintf(stderr, "annotating... %ld/%ld\n", *computed, index->size);
   fprintf(stderr, "[proc] annotating... done [%.3fs]\n", ((clock()-clk)*1.0/CLOCKS_PER_SEC)/threads);
   // Free variables.
   free(jobs);
   free(mutex);
   free(monitor);
   return (annotation_t){.bitfield = repeat_bf};
}

int
next_seq
(
 int32_t    qlen,
 int32_t    tau,
 int64_t    max_sa,
 int64_t  * sa_ptr,
 uint8_t  * query,
 fmdpos_t * pos,
 index_t  * index
)
{
   int n_cnt;
   int trail;
   uint8_t * tmp = malloc(qlen*sizeof(uint8_t));
   memcpy(tmp, query, qlen*sizeof(uint8_t));

   // Iterate until a valid sequence is found (num(N) <= tau).
   do {
      // Reset variables.
      n_cnt = 0;
      trail = 0;
      // Get sequence from genome.
      char * seq = index->genome + get_sa(*sa_ptr, index->sar);
      // Iterate over qlen nucleotides.
      for (int i = 0; i < qlen; i++) {
         int nt = translate[(int)query[i]];
         n_cnt += nt == UNKNOWN_BASE;
         // Extend BWT search.
         if (trail == i && nt == tmp[i]) 
            trail++;
         else
            pos[i+1] = extend_fw(nt, pos[i], index->bwt);

         // Update tmp query.
         tmp[i] = nt;
      }
      // Update sa_ptr (point the start of the next suffix).
      *sa_ptr += pos[qlen].sz;
      
      // Return -1 if we're already out of segment.
      if (pos[qlen].fp >= max_sa) {
         free(tmp);
         return -1;
      }

   } while (n_cnt > tau);

   // Recompute trail wrt last query.
   trail = 0;
   while (tmp[trail] == query[trail]) trail++;

   // Copy tmp to query.
   memcpy(query, tmp, qlen*sizeof(uint8_t));

   // Free memory.
   free(tmp);

   // Return trail.
   return trail;
}

void
store_hits
(
 annjob_t  * job,
 pstree_t  * stack_tree
)
{
   // TODO.
}

inline
void
send_progress_to_scheduler
(
 int64_t     sa_pos,
 int64_t   * computed,
 annjob_t  * job
)
{
   // Send progress every 1M processed sequences.
   if (__builtin_expect((computed>>20) > (sa_pos>>20),0)) {
      pthread_mutex_lock(job->mutex);
      *(job->computed) += sa_pos - *computed;
      pthread_cond_signal(job->monitor);
      pthread_mutex_unlock(job->mutex);
      *computed = sa_pos;
   }
   return;
}

void *
annotate_mt
(
 void * argp
 )
{
   // Read job.
   annjob_t * job    = (annjob_t *)argp;
   index_t  * index  = job->index;
   int64_t    sa_ptr = job->beg;
   int        kmer   = job->kmer;
   int        tau    = job->tau;

   // Translated query.
   uint8_t * query = malloc(kmer*sizeof(uint8_t));
   memset(query, NUM_BASES, kmer*sizeof(uint8_t));

   // BWT path.
   fmdpos_t * path = malloc((kmer+1)*sizeof(fmdpos_t));
   path[0] = index->bwt->fmd_base;

   // Hit stack.
   pstree_t * stack_tree = alloc_stack_tree(tau);

   // Aux vars.
   int64_t computed = sa_ptr;

   // Iterate over all kmers.
   while (sa_ptr < job->end) {
      // Report progress to scheduler.
      send_progress_to_scheduler(sa_ptr, &computed, job);
      // This updates query and path, and makes sa_ptr point to the position of the next sequence in the SA.
      int trail = next_seq(kmer, tau, job->end, &sa_ptr, query, path, index);
      // Query the sequence.
      blocksearch_trail_sc(query, path, kmer, tau, trail, index, stack_tree);
      // Update annotation info.
      store_hits(job, stack_tree);
   }

   // Report job end to scheduler.
   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += job->end - computed;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free memory.
   free_stack_tree(stack_tree);
   free(query);
   free(path);

   return NULL;

}

/*
   // aux vars.
   uint64_t sa_pos = job->beg;
   uint64_t computed = job->beg;
   uint8_t * query = malloc(kmer);
   uint8_t * lastq = malloc(kmer);
   uint64_t kmers = 0;
   uint64_t unique = 0;

   int w_bits = 0;
   while ((tau >> w_bits) > 0) w_bits++;
   // Two extra bits to store log10(num of different neighbors).
   w_bits += 2;

   // Force first trail to be 0.
   memset(lastq, 5, kmer);

   // Iterate over all kmers.
   while (sa_pos < job->end) {

      // Get query values.
      uint64_t locus = get_sa(sa_pos, index->sar);
      char * seq = index->genome + locus;
      // Translate query sequence.
      for (int i = 0; i < kmer; i++) query[i] = translate[(int)seq[i]];
      // Find a non-reverse duplicate sequence without 'N'.
      while (contains_n(query, kmer) || reverse_duplicate(query, kmer)) {
         // Find next kmer.
         sa_pos++;
         if (sa_pos >= job->end) goto free_and_return;
         locus = get_sa(sa_pos, index->sar);
         seq = index->genome + locus;
         // Translate query sequence.
         for (int i = 0; i < kmer; i++) query[i] = translate[(int)seq[i]];
      }
      // Increase kmer count.
      kmers++;
      // Compute trail and update lastq.
      int trail = 0;
      while (query[trail] == lastq[trail] && trail < kmer) trail++;
      memcpy(lastq, query, kmer);

      // Alloc count buffers.
      int64_t * loci_count = calloc(tau+1,sizeof(int64_t));
      int64_t * hit_count = calloc(tau+1,sizeof(int64_t));
      // Search sequence.
      blocksearch_trail(query, kmer, tau, trail, index, stack_tree);
      // Count hits.
      pathstack_t * hits = stack_tree->stack;
      fmdpos_t   self = {.fp = 0, .rp = 0, .sz = 1};
      for (int i = 0; i < hits->pos; i++) {
         spath_t p = hits->path[i];
         loci_count[p.score] += p.pos.sz;
         hit_count[p.score]++;
         if (p.score == 0)
            self = p.pos;
      }
      loci_count[0] -= 1;
      
      // Store annotation.
      // Annotation format.
      // - Exact match: 0b000..0
      // - No match (d<=tau): 0b110..0
      // - Otherwise:
      //   bit[n-1] bit[n-2] = log10(neighbor diversity)
      //   bit[n-3].. bit[0] = matching tau.
      if (!loci_count[0]) {
         // Set locus in forward strand.
         if (locus >= index->size/2) locus = index->size - 1 - locus - kmer;
         // Words.
         size_t bit = w_bits*locus;
         size_t word = bit >> 3;
         uint8_t shift = bit & 7;
         uint16_t uw = 0, lw = 0, w = 0;
         // Find match distance.
         // Find match distance.
         int dist = tau+1;
         for (int d = 1; d <= tau; d++) {
            if (loci_count[d]) {
               dist = d;
               break;
            }
         }
         // No match (write 0b110..0).
         if (dist > tau) {
            w = 3 << (w_bits-2);
         } else {
            w = dist;
            if      (hit_count[dist] < 5)   {}
            else if (hit_count[dist] < 50)  w |= 1 << (w_bits-2);
            else if (hit_count[dist] < 500) w |= 2 << (w_bits-2);
            else                            w |= 3 << (w_bits-2);
         }
         lw = (w << shift) & 0x00FF;
         uw = (w >> (8-shift)) & 0x00FF;
         // Mutex write.
         // TODO:
         // Since each thread is workin on a different
         // region of the SA, this can be done safely
         // without mutex lock.
         pthread_mutex_lock(job->mutex);
         job->repeat_bf[word]   |= lw;
         job->repeat_bf[word+1] |= uw;
         pthread_mutex_unlock(job->mutex);
      }
      free(loci_count);
      free(hit_count);
      // Update position.
      sa_pos = self.fp + self.sz;
   }

   free_and_return:
   // Report to scheduler.
   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += job->end - computed;
   *(job->kmers) += kmers;
   *(job->unique) += unique;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free memory.
   free_stack_tree(stack_tree);
   free(query);
   free(lastq);

   return NULL;
}
*/
