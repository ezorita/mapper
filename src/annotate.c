#include "annotate.h"
#include <string.h>
#include <time.h>


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
         .beg = pos.fp,
         .end = pos.fp + pos.sz
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

void
locus_annotation
(
 uint8_t * info,
 uint8_t * ann,
 int       alnsize
)
{
   // Annotation data structure:
   // 8 bits per genomic locus.
   //   bit 7 - Locus alignment info.
   //   bit 6 - Locus alignment flag.
   //   bit 5 *
   //   bit 4 - Distance to closest neighbor (00: 1, 01: 2, 10: 3, 11: 4).
   //   bit 3 *
   //   bit 2 *
   //   bit 1 * 
   //   bit 0 - Neighbor count
   uint16_t  cnt = *((uint16_t *)info);
   uint8_t   tau = *(info + sizeof(uint16_t));
   uint8_t * aln = info + sizeof(uint16_t) + sizeof(uint8_t);

   if (cnt == 0 || cnt == NO_INFO) return;

   // Count.
   if      (cnt <= 10)               *ann |= (uint8_t)cnt;
   else if (cnt > 10 && cnt <= 20)   *ann |= 0x0B;
   else if (cnt > 20 && cnt <= 50)   *ann |= 0x0C;
   else if (cnt > 50 && cnt <= 100)  *ann |= 0x0D;
   else if (cnt > 100 && cnt <= 500) *ann |= 0x0E;
   else if (cnt > 500)               *ann |= 0x0F;
   
   // Tau.
   *ann |= ((tau-1) & 0x03) << 4;
   
   // Alignment info.
   if (*aln != 255) {
      // Set alignment flag in locus.
      *ann |= (uint8_t)1 << 6;
      // Store alignment info.
      for (int i = 0; i < alnsize && aln[i]; i++)
         // Alignment values are 1-based.
         *(ann + aln[i]-1) |= (uint8_t)1 << 7;
   }

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

   uint8_t * tmp_info = calloc(index->size, word_size*sizeof(uint8_t));

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
      jobs[i].kmer     = (uint32_t) kmer;
      jobs[i].tau      = (uint16_t) tau;
      jobs[i].wsize    = (uint16_t) word_size;
      jobs[i].computed = &computed;
      jobs[i].done     = &done;
      jobs[i].mutex    = mutex;
      jobs[i].monitor  = monitor;
      jobs[i].index    = index;
      jobs[i].info     = tmp_info;
   }

   clock_t clk = clock();
   // Run threads.
   for (int i = 0; i < min(threads, num_jobs); i++) {
      pthread_t thread;
      pthread_create(&thread,NULL,annotate_mt,jobs+i);
      pthread_detach(thread);
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
      fprintf(stderr, "\r[proc] computing neighbors... %.2f%%", computed*100.0/index->bwt->fmd_base.sz);
      // Wait for signal.
      pthread_cond_wait(monitor,mutex);
   }
   pthread_mutex_unlock(mutex);

   // Report time.
   double elapsed = ((clock()-clk)*1.0/CLOCKS_PER_SEC);
   fprintf(stderr, "\r[proc] computing neighbors... done [%.3fs with %d threads (%.3fs/thread)]\n", elapsed, threads, elapsed/threads);
   fprintf(stderr, "[proc] compressing annotation...\n");
   clk = clock();

   // Alloc annotation.
   annotation_t ann;
   ann.kmer = kmer;
   ann.tau  = tau;
   ann.info = calloc(index->size/2+1, sizeof(uint8_t));
   ann.size = index->size/2+1;

   // Compute and store annotation from info.
   uint64_t i = 0;
   uint8_t * info = tmp_info;
   while (1) {
      // info structure:
      //   bytes     0-1: number of different neighbors (16bit).
      //   byte        2: distance to closest neighbor.
      //   byte  3-wsize: bytes encoding the mutated positions to find the closest neighbors.

      while (i < index->size && (*(uint16_t *)info == NO_INFO || *(uint16_t *)info == 0)) {
         i++;
         info += word_size;
      }
      if (i >= index->size) break;

      // Compute k-mer range in SA.
      uint64_t size = 1;
      while (i+size < index->size && *(uint16_t *)(info + size*word_size) == 0)
         size++;
      // Get SA values.
      size_t * range = malloc(size*sizeof(size_t));
      get_sa_range(i, size, range, index->sar);
      for (uint64_t j = 0; j < size; j++) {
         // Store annotation in forward strand.

         if (range[j] >= index->size/2) {
            // Convert reverse to fw strand.
            range[j] = index->size - range[j] - kmer - 1;
            // Reverse positions of bits set to '1'.
            uint8_t * aln = info + (sizeof(uint16_t) + sizeof(uint8_t));
            for (int k = 0; k < word_size - 3 && aln[k]; k++)
               aln[k] = kmer + 1 - aln[k];
         }
         locus_annotation(info, ann.info + range[j], word_size-3);
      }
      free(range);

      // Update pointer.
      i    += size;
      info += size*word_size;
   }
   
   // Free variables.
   free(jobs);
   free(mutex);
   free(monitor);
   free(tmp_info);
   return ann;
}

int
next_seq
(
 int32_t    qlen,
 int32_t    tau,
 int64_t    sa_ptr,
 int64_t  * next_sa,
 int32_t  * last_fragment,
 uint8_t  * query_1,
 uint8_t  * query_2,
 fmdpos_t * pos,
 index_t  * index
)
{
   int n_cnt = 0;
   int trail_1 = 0, trail_2 = 0;
   uint8_t * tmp = malloc(qlen*sizeof(uint8_t));
   // Get sa_ptr sequence from genome.
   char * seq = index->genome + get_sa(sa_ptr, index->sar);
   // Iterate over qlen nucleotides.
   for (int i = 0; i < qlen; i++) {
      int nt = translate[(int)seq[i]];
      // Report sequences containing '$'.
      if (nt == 5) {
         pos[qlen].sz = 0;
         break;
      }
      n_cnt += nt == UNKNOWN_BASE;
      // Track trail.
      if (trail_1 == i && nt == query_1[i])
         trail_1++;
      if (trail_2 == i && nt == query_2[i])
         trail_2++;
      // Extend BWT search.
      pos[i+1] = extend_fw(nt, pos[i], index->bwt);

      // Update tmp query.
      tmp[i] = nt;
   }

   // Check whether sequence is valid (all in same strand).
   if (pos[qlen].sz == 0) {
      *next_sa = sa_ptr + 1;
      free(tmp);
      return -1;
   }

   // Update sa_ptr (point the start of the next suffix).
   *next_sa = sa_ptr + pos[qlen].sz;

   // Check N count.
   if (n_cnt > tau) {
      free(tmp);
      return -1;
   }

   // Copy tmp to query.
   if (pos[qlen].fp >= pos[qlen].rp) {
      *last_fragment = 1;
      memcpy(query_2, tmp, qlen*sizeof(uint8_t));
   } else {
      *last_fragment = 0;
      memcpy(query_1, tmp, qlen*sizeof(uint8_t));
   }

   // Free memory.
   free(tmp);

   // Return trail.
   return *last_fragment ? trail_2 : trail_1;
}

void
merge_alignments
(
 uint8_t * a,
 uint8_t * b,
 int       len
)
{
   uint8_t * tmp = calloc(2*len, sizeof(uint8_t));

   // Merge values, halt if a or b equals 0.
   int i = 0, j = 0, k = 0;
   while (i < len && j < len && k <= len) {
      if (a[i] == 0 || b[j] == 0) break;
      if (a[i] == b[j]) {
         tmp[k++] = a[i];
         i++; j++;
      }
      else if (a[i] < b[j])
         tmp[k++] = a[i++];
      else if (a[i] > b[j])
         tmp[k++] = b[j++];
   }

   // Fill in remaining values.
   while (i < len && k <= len && a[i])
      tmp[k++] = a[i++];
   while (j < len && k <= len && b[j])
      tmp[k++] = b[j++];
   while (k < len)
      tmp[k++] = 0;

   // If we found too many align positions, flag 0xFF.
   if (k > len) 
      memset(a, 0xFF, len*sizeof(uint8_t));
   else 
      memcpy(a, tmp, len*sizeof(uint8_t));

   // Return.
   free(tmp);
   return;
}

void
compute_aln_positions
(
 uint64_t * bits,
 uint8_t  * pos,
 int        nbits,
 int        npos,
 int        reverse
)
{
   uint8_t * tmp = calloc(npos+1, sizeof(uint8_t));
   int a = 0;

   // Store positions of bits set to '1'.
   for (int i = 0; i < nbits && a <= npos; i++) {
      int w = i / ALIGN_WORD_SIZE;
      int b = i % ALIGN_WORD_SIZE;
      if ((bits[w] >> b) & (uint64_t)1)
         tmp[a++] = (reverse ? nbits - i : i + 1);
   }
   // If we found too many align positions, flag 0xFF.
   if (a > npos) {
      memset(pos, 0xFF, npos*sizeof(uint8_t));
   } else {
      memcpy(pos, tmp, npos*sizeof(uint8_t));
   }
   
   free(tmp);
   return;
}

void
store_hits
(
 annjob_t     * job,
 pathstack_t  * stack,
 fmdpos_t       query
)
// Stores new information discovered in the last query.
// Information is always stored/updated in the lexicographically
// smallest between the queried/hit sequence and its rev complement.
{
   int hits     = 0;
   int tau      = job->tau + 1;
   int aln_size = job->wsize - 3;

   // Flag lexicographically greater as NO_INFO.
   if (query.fp > query.rp) {
      pthread_mutex_lock(job->mutex);
      *(uint16_t *)(job->info + query.fp * job->wsize) = NO_INFO;
      pthread_mutex_unlock(job->mutex);
   }

   // No matches (self hit only).
   if (stack->pos < 2) {
      pthread_mutex_lock(job->mutex);
      // If there is no previous info, set smaller as NO_INFO.
      if (*(uint16_t *)(job->info + min(query.fp,query.rp) * (size_t)(job->wsize)) == 0)
         *(uint16_t *)(job->info + min(query.fp,query.rp) * (size_t)(job->wsize)) = NO_INFO;
      pthread_mutex_unlock(job->mutex);
      return;
   }

   // Aggregate query alignment.
   uint64_t * qalign = calloc(ALIGN_WORDS, sizeof(uint64_t));

   // Store hits info in their SA start position.
   for (int i = 0; i < stack->pos; i++) {
      spath_t path = stack->path[i];
      // Remove self hit.
      if (path.pos.fp == query.fp) continue;

      // Select lexicographically smallest between the neighbor and its revcomp.
      int64_t nptr = (path.pos.rp < path.pos.fp ? path.pos.rp : path.pos.fp);
      int     nrev = path.pos.rp < path.pos.fp;

      // info structure:
      //   bytes     0-1: number of different neighbors (16bit).
      //   byte        2: distance to closest neighbor.
      //   byte  3-wsize: bytes encoding the mutated positions to find the closest neighbors.
      uint8_t  * ninfo = job->info + nptr * job->wsize;
      uint16_t * ncnt  = (uint16_t *) ninfo;
      uint8_t  * ntau  = ninfo + sizeof(uint16_t);
      uint8_t  * naln  = ntau + 1;

      // Remote updates (neighbor).
      pthread_mutex_lock(job->mutex);
      // Update neighbor.
      if (*ntau == path.score) {
         // Add one to neighbor count.
         *ncnt = (uint16_t) min(0xFFFE, ((uint32_t) *ncnt) + 1);
         // Check if neighbor has free align slots.
         if (*naln != 0xFF) {
            // Compute mismatch positions.
            uint8_t hit_aln[aln_size];
            compute_aln_positions(path.align, hit_aln, job->kmer, aln_size, nrev);
            // Merge mismatch positions and store in neighbor memory (naln).
            merge_alignments(naln, hit_aln, aln_size);
         }
      }
      // Reset neighbor.
      else if (*ncnt == 0 || *ncnt == NO_INFO || *ntau > path.score) {
         // Set count to 1 and neighbor's tau to path.score.
         *ncnt = 1;
         *ntau = path.score;
         // Compute mismatch positions.
         uint8_t hit_aln[aln_size];
         compute_aln_positions(path.align, hit_aln, job->kmer, aln_size, nrev);
         // Store mismatch positions.
         memcpy(naln, hit_aln, aln_size);
      }
      pthread_mutex_unlock(job->mutex);

      // Local updates (query).
      // Update query info.
      if (path.score == tau) {
         // Update query align positions.
         for (int w = 0; w < ALIGN_WORDS; w++)
            qalign[w] |= path.align[w];
         // Update hit count.
         hits++;
      } 
      // Found score smaller than tau, update.
      else if (path.score < tau) {
         memcpy(qalign, path.align, ALIGN_WORDS*sizeof(uint64_t));
         tau = path.score;
         hits = 1;
      }
   }

   // Select lexicographically smallest between the query and its revcomp.
   int64_t qptr = (query.rp < query.fp ? query.rp : query.fp);
   int     qrev = query.rp < query.fp;

   // Query info.
   uint8_t  * info = job->info + qptr * job->wsize;
   uint16_t * qcnt = (uint16_t *) info;
   uint8_t  * qtau = info + sizeof(uint16_t);
   uint8_t  * qaln = qtau + 1;

   // If query is empty or tau < qtau, reset query.
   pthread_mutex_lock(job->mutex);
   if (*qcnt == 0 || *qcnt == NO_INFO || *qtau > tau) {
      // Set tau and neighbor count.
      *qtau = (uint8_t)  tau;
      *qcnt = (uint16_t) min(NO_INFO - 1, hits);
      // Compute and directly store mismatch positions for query sequence.
      compute_aln_positions(qalign, qaln, job->kmer, aln_size, qrev);
   }
   // Update query.
   else if (*qtau == tau) {
      // Update hit count.
      *qcnt = (uint16_t) min(NO_INFO-1, ((uint32_t)hits)+((uint32_t)*qcnt));
      // Check if query has free align slots.
      if (*qaln != 0xFF) {
         // Compute mismatch positions for query sequence.
         uint8_t tmp_aln[aln_size];
         compute_aln_positions(qalign, tmp_aln, job->kmer, aln_size, qrev);
         // Update alignments in query by merging.
         merge_alignments(qaln, tmp_aln, aln_size);
      }
   }
   pthread_mutex_unlock(job->mutex);

   free(qalign);
}

void
send_progress_to_scheduler
(
 int64_t     sa_pos,
 int64_t   * computed,
 annjob_t  * job
)
{
   // Send progress every 1M processed sequences.
   if (__builtin_expect(( (*computed) >> 20) > (sa_pos>>20) , 0 )) {
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
   uint8_t * query_1 = malloc(kmer*sizeof(uint8_t));
   uint8_t * query_2 = malloc(kmer*sizeof(uint8_t));
   memset(query_1, NUM_BASES, kmer*sizeof(uint8_t));
   memset(query_2, NUM_BASES, kmer*sizeof(uint8_t));
   // BWT path.
   fmdpos_t * path = malloc((kmer+1)*sizeof(fmdpos_t));
   path[0] = index->bwt->fmd_base;

   // Hit stack.
   pstree_t * stack_tree_1 = alloc_stack_tree(tau);
   pstree_t * stack_tree_2 = alloc_stack_tree(tau);

   // Aux vars.
   int64_t computed = sa_ptr;
   int64_t next_sa;
   int32_t last_fragment;

   // Iterate over all kmers.
   while (sa_ptr < job->end) {
      // Report progress to scheduler.
      send_progress_to_scheduler(sa_ptr, &computed, job);
      // This updates query and path, and makes sa_ptr point to the position of the next sequence in the SA.
      int trail = next_seq(kmer, tau, sa_ptr, &next_sa, &last_fragment, query_1, query_2, path, index);
      // Check valid query sequence.
      if (trail < 0) {
         // Flag 'sa_ptr' as NO_INFO.
         *(uint16_t *)(job->info + sa_ptr * job->wsize) = NO_INFO;
      } else {
         // Query the sequence.
         if (last_fragment)
            blocksc_trail(query_2, path, kmer, tau, trail, index, stack_tree_2);
         else
            blocksc_trail(query_1, path, kmer, tau, trail, index, stack_tree_1);
         // Update annotation info.
         store_hits(job, (last_fragment ? stack_tree_2->stack : stack_tree_1->stack), path[kmer]);
      }
      sa_ptr = next_sa;
   }

   // Report job end to scheduler.
   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += job->end - computed;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free memory.
   free_stack_tree(stack_tree_1);
   free_stack_tree(stack_tree_2);
   free(query_1);
   free(query_2);
   free(path);

   return NULL;

}
