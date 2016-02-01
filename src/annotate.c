#include "annotate.h"
#include <string.h>
#include <time.h>

int
reverse_duplicate
(
 uint8_t * query,
 int       kmer
 )
/*
** This function takes a seq from the genome and compares
** it with its reverse complement. Returns 1 (duplicate) if
** the reverse complement is lexicographically smaller.
** The aim of this is to avoid computing the same thing twice.
** Since we have an index with forward and reverse strands,
** the forward and reverse complement of each kmer must yield
** exactly the same number of hits. Therefore we will only
** compute and store the lexicographically smallest of each
** forward/reverse pair.
*/
{
   for (int i = 0, j = kmer-1; i < kmer; i++, j--) {
      if (3-query[j] == query[i]) continue;
      return 3-query[j] < query[i];
   }
   return 0;
}

int
contains_n
(
 uint8_t * query,
 int       kmer
)
{
   for (int i = 0; i < kmer; i++)
      if (query[i] > 3) return 1;
   return 0;
}

int
cmp_seq
(
 char * seq_a,
 char * seq_b,
 int    k
)
{
   for (int i = 0; i < k; i++)
      if (translate[(int)seq_a[i]] != translate[(int)seq_b[i]]) return 0;
   return 1;
}


annotation_t
annotate
(
 int        kmer,
 int        tau,
 int        repeat_thr,
 index_t  * index,
 int        threads
 )
{

   // Alloc variables.
   pthread_mutex_t * mutex = malloc(sizeof(pthread_mutex_t));
   pthread_mutex_t * ht_mutex = malloc(sizeof(pthread_mutex_t));
   pthread_cond_t  * monitor = malloc(sizeof(pthread_cond_t));
   uint64_t computed = 0;
   uint64_t kmers = 0;
   uint64_t collision = 0;
   int32_t  done = 0;
   uint8_t  * repeat_bf = malloc(((index->size >> 4) + 1) * sizeof(uint8_t));
   // Hash table.
   int bits = 0;
   while ((index->size >> bits) > 0) bits++;
   // Alloc Hash table with same size as genome.
   htable_t * ht = htable_new(bits+2);
   // Verbose structure size.
   fprintf(stderr, "structure size: [hash table] %.2f MB\t[bit field] %.2f MB\n", ((sizeof(htable_t) + ((uint64_t)1<<bits))*1.0/1024)/1024, (((index->size >> 4) + 1)*1.0/1024)/1024);

   // Initialization.
   pthread_mutex_init(mutex,NULL);
   pthread_mutex_init(ht_mutex,NULL);
   pthread_cond_init(monitor,NULL);

   // Find 'NXXX..' interval.
   bwpos_t n_pos = {.sp = 0, .ep = index->size-1, .depth = 0};
   uint64_t eff_size = index->size;
   suffix_extend(translate['N'], n_pos, &n_pos, index);
   if (n_pos.ep >= n_pos.sp) eff_size = n_pos.sp;
   // Compute thread jobs.
   int      njobs = 10*threads;
   uint64_t seqs  = eff_size/njobs;
   annjob_t ** job = malloc(njobs*sizeof(annjob_t*));
   for (int i = 0; i < njobs; i++) {
      // Allocate job.
      job[i] = malloc(sizeof(annjob_t));
      // Compute sequence offset in BWT.
      uint64_t offset = i*seqs;
      // Set job values.
      job[i]->beg = offset;
      job[i]->mutex = mutex;
      job[i]->ht_mutex = ht_mutex;
      job[i]->monitor = monitor;
      job[i]->computed = &computed;
      job[i]->done = &done;
      job[i]->index = index;
      job[i]->kmer = kmer;
      job[i]->tau = tau;
      job[i]->kmers = &kmers;
      job[i]->collision = &collision;
      job[i]->repeat_thr = repeat_thr;
      job[i]->repeat_bf = repeat_bf;
      job[i]->htable = ht;
   }
   // Set job end indices.
   for (int i = 0; i < njobs-1; i++) job[i]->end = job[i+1]->beg;
   job[njobs-1]->end = eff_size;

   // Do not break intervals between jobs.
   for (int i = 1; i < njobs; i++) {
      while (cmp_seq(index->genome + get_sa(job[i-1]->end - 1,index->sa,index->sa_bits), index->genome + get_sa(job[i]->beg, index->sa, index->sa_bits), kmer))
         job[i]->beg += 1;
   }

   // Run threads.
   for (int i = 0; i < threads; i++) {
      pthread_t thread;
      pthread_create(&thread,NULL,annotate_mt,job[i]);
      pthread_detach(thread);
      //      annotate_mt(job[i]);
   }
   // Sleep and wait for thread signals.
   int runcount = threads;
   pthread_mutex_lock(mutex);
   while (done < njobs) {
      // Run new threads.
      for (int i = runcount; i < min(done + threads, njobs); i++) {
         pthread_t thread;
         pthread_create(&thread,NULL,annotate_mt,job[i]);
         pthread_detach(thread);
      }
      // Update submitted jobs.
      runcount = done + threads;
      // Report progress.
      fprintf(stderr, "annotating... %.2f%%\r", computed*100.0/eff_size);
      // Wait for signal.
      pthread_cond_wait(monitor,mutex);
   }
   pthread_mutex_unlock(mutex);

   //   fprintf(stderr, "annotating... %ld/%ld\n", *computed, index->size);
   fprintf(stderr, "annotating... %.2f%%\n", computed*100.0/eff_size);
   fprintf(stderr, "kmers: %ld different kmers out of %ld processed loci (%.2f%%).\n",kmers,index->size/2,kmers*200.0/index->size);
   fprintf(stderr, "hash table collisions: %ld (%.2f%%)\n", collision, collision*100.0/kmers);

   // Free variables.
   for (int i = 0; i < threads; i++) { free(job[i]); }
   free(job);
   free(mutex);
   free(ht_mutex);
   free(monitor);
   return (annotation_t){.htable = ht, .bitfield = repeat_bf};
}

void *
annotate_mt
(
 void * argp
 )
{
   annjob_t * job = (annjob_t *)argp;
   int kmer = job->kmer;
   int tau  = job->tau;
   int thr  = job->repeat_thr;
   index_t * index = job->index;
   // hit stack.
   pstree_t * stack_tree = alloc_stack_tree(tau);
   // aux vars.
   uint64_t sa_pos = job->beg;
   uint64_t computed = job->beg;
   uint8_t * query = malloc(kmer);
   uint8_t * lastq = malloc(kmer);
   uint64_t kmers = 0;
   uint64_t ht_collision = 0;

   // Force first trail to be 0.
   memset(lastq, 5, kmer);

   // Iterate over all kmers.
   while (sa_pos < job->end) {
      // Report progress to scheduler.
      if ((computed>>20) > (sa_pos>>20)) {
         pthread_mutex_lock(job->mutex);
         *(job->computed) += sa_pos - computed;
         pthread_cond_signal(job->monitor);
         pthread_mutex_unlock(job->mutex);
         computed = sa_pos;
      }

      // Get query values.
      uint64_t locus = get_sa(sa_pos, index->sa, index->sa_bits);
      char * seq = index->genome + locus;
      // Translate query sequence.
      for (int i = 0; i < kmer; i++) query[i] = translate[(int)seq[i]];
      // Find a non-reverse duplicate sequence without 'N'.
      while (contains_n(query, kmer) || reverse_duplicate(query, kmer)) {
         // Find next kmer.
         sa_pos++;
         if (sa_pos >= job->end) goto free_and_return;
         locus = get_sa(sa_pos, index->sa, index->sa_bits);
         seq = index->genome + locus;
         // Translate query sequence.
         for (int i = 0; i < kmer; i++) query[i] = translate[(int)seq[i]];
      }
      // Increase kmer count.
      kmers++;
      // Compute query key.
      uint64_t key = XXH64(query, kmer, 0);
      // Compute trail and update lastq.
      int trail = 0;
      while (query[trail] == lastq[trail] && trail < kmer) trail++;
      memcpy(lastq, query, kmer);
      // Search sequence.
      blocksearch_trail(query, kmer, tau, trail, index, stack_tree);
      // Count hits.
      uint64_t hit_count = 0;
      pathstack_t * hits = stack_tree->stack;
      fmdpos_t self = {.fp = 0, .rp = 0, .sz = 1};
      for (int i = 0; i < hits->pos; i++) {
         spath_t p = hits->path[i];
         if (p.score == 0)
            self = p.pos;
         hit_count += p.pos.sz;
      }
      if (hit_count > 1 && hit_count <= thr) hit_count = 2;
      else if(hit_count > thr) hit_count = 3;
      // Store count.
      pthread_mutex_lock(job->ht_mutex);
      // Set unique bit field.
      if (hit_count == 1) {
         if (locus >= index->size/2) locus = index->size - 1 - locus - kmer;
         job->repeat_bf[locus>>3] |= 1 << (locus&7);
      }
      // Store hit count in hash table.
      int ht_val = htable_get(job->htable, key);
      if (!ht_val || ht_val > hit_count)
         htable_set(job->htable, key, hit_count);
      // DEBUG.
      ht_collision += (ht_val && ht_val != hit_count);
      // DEBUG END.
      pthread_mutex_unlock(job->ht_mutex);
      // Update position.
      sa_pos = self.fp + self.sz;
   }

   free_and_return:
   // Report to scheduler.
   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += job->end - computed;
   *(job->kmers) += kmers;
   *(job->collision) += ht_collision;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free memory.
   free_stack_tree(stack_tree);

   return NULL;
}


/*
htable_t *
htable_new
(
 int    index_bits,
 int    probe_bits,
 size_t elm_size
 )
{
   htable_t * htable = calloc(sizeof(htable_t) + elm_size * (1 << index_bits), 1);
   if (htable == NULL) {
      return NULL;
   }
   // Compute bitmasks.
   htable->index_bits = index_bits;
   htable->index_bitmask = 0xFFFFFFFFFFFFFFFF >> (64 - index_bits);
   htable->probe_bits = probe_bits;
   htable->probe_bitmask = (0xFFFFFFFFFFFFFFFF >> (64 - probe_bits)) << index_bits;
   htable->probe_datamask = 0xFF >> (8 - probe_bits);
   // Set element size and reset count.
   htable->count = 0;
   htable->elm_size = elm_size;
   return htable;
}

int
htable_insert
(
 uint8_t         * key,
 uint32_t          keylen,
 uint8_t           data,
 uint32_t          data_pos,
 htable_t        * table,
 pthread_mutex_t * mutex
 )
{
   if (data_pos >= table->elm_size) return -1;
   if (table->count >= table->index_bitmask) return -1;
   // Generate key Hash.
   uint64_t k = XXH64(key, keylen, 0);
   // Split index and probe bits.
   uint64_t idx = k & table->index_bitmask;
   uint8_t chk = ((k & table->probe_bitmask) >> table->index_bits) & 0xFF;
   // chk = 0 reserved for free cells.
   if (chk == 0) chk = 1;
   // Find spot in table. Simple linear probing.
   uint8_t * beg_addr = table->c + idx*table->elm_size + data_pos;
   uint8_t * max_addr = table->c + table->index_bitmask * table->elm_size;
   uint8_t * cell = beg_addr;
   // This part of the code is concurrent, mutex lock required.
   pthread_mutex_lock(mutex);
   while (*cell != 0 && ((*cell & table->probe_datamask) != chk)) {
      // Next cell (Linear probing).
      cell += table->elm_size;
      // Loop until cell is free or has the same probe.
      if (cell > max_addr) cell = table->c + data_pos;
   }
   // Store data in cell.
   *cell = data;
   // End of concurrent code.
   pthread_mutex_unlock(mutex);

   return 0;
}
*/
