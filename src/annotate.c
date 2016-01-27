#include "annotate.h"
#include <string.h>
#include <time.h>

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

int
reverse_duplicate
(
 char * seq,
 int    kmer
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
      if (translate_rc[(int)seq[j]] == translate[(int)seq[i]]) continue;
      return translate_rc[(int)seq[j]] < translate[(int)seq[i]];
   }
   return 0;
}

int
contains_n
(
 char * seq,
 int   kmer
)
{
   for (int i = 0; i < kmer; i++)
      if (translate[(int)seq[i]] > 3) return 1;
   return 0;
}


int
annotate
(
 int        kmer,
 int        tau,
 index_t  * index,
 int        threads
 )
{

   // Alloc variables.
   pthread_mutex_t * mutex = malloc(sizeof(pthread_mutex_t));
   pthread_cond_t  * monitor = malloc(sizeof(pthread_cond_t));
   uint64_t * computed = malloc(sizeof(uint64_t));
   uint64_t * kmers = malloc(sizeof(uint64_t));
   uint8_t * counts = malloc(index->size * sizeof(uint8_t));
   int * done = malloc(sizeof(int));

   // Initialization.
   *computed = 0;
   *done = 0;
   *kmers = 0;
   pthread_mutex_init(mutex,NULL);
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
      job[i]->monitor = monitor;
      job[i]->computed = computed;
      job[i]->done = done;
      job[i]->index = index;
      job[i]->kmer = kmer;
      job[i]->tau = tau;
      job[i]->kmers = kmers;
      job[i]->counts = counts;
   }
   // Set job end indices.
   for (int i = 0; i < njobs-1; i++) job[i]->end = job[i+1]->beg;
   job[njobs-1]->end = eff_size;


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
   while (*done < njobs) {
      // Run new threads.
      for (int i = runcount; i < min(*done + threads, njobs); i++) {
         pthread_t thread;
         pthread_create(&thread,NULL,annotate_mt,job[i]);
         pthread_detach(thread);
      }
      // Update submitted jobs.
      runcount = *done + threads;
      // Report progress.
      fprintf(stderr, "annotating... %.2f%%\r", *computed*100.0/index->size);
      // Wait for signal.
      pthread_cond_wait(monitor,mutex);
   }
   pthread_mutex_unlock(mutex);

   //   fprintf(stderr, "annotating... %ld/%ld\n", *computed, index->size);
   fprintf(stderr, "annotating... %.2f%%\n", *computed*100.0/index->size);
   fprintf(stderr, "kmers: %ld unique kmers out of %ld genomic loci.\n",*kmers,index->size-1);

   for (int i = 0; i < threads; i++) { free(job[i]); }
   free(job);
   free(mutex);
   free(monitor);
   free(computed);
   free(done);
   free(kmers);
   return 0;
}

void *
annotate_mt
(
 void * argp
 )
{
   uint64_t computed = 0, last = 0;
   annjob_t * job = (annjob_t *)argp;
   int kmer = job->kmer;
   int tau  = job->tau;
   index_t * index = job->index;
   // hit stack.
   pathstack_t * hits = pathstack_new(PATHSTACK_DEF_SIZE);
   // aux vars.
   uint64_t sa_pos = job->beg;
   uint8_t * query = malloc(kmer);
   uint64_t kmers = 0;

   // Iterate over all kmers.
   while (sa_pos < job->end) {
      // Increase kmer count.
      kmers++;
      
      // Report progress to scheduler.
      if ((computed>>20) > (last>>20)) {
         pthread_mutex_lock(job->mutex);
         *(job->computed) += computed - last;
         pthread_cond_signal(job->monitor);
         pthread_mutex_unlock(job->mutex);
         last = computed;
      }

      // Get query values.
      uint64_t locus = get_sa(sa_pos, index->sa, index->sa_bits);
      char * seq = index->genome + locus;
      uint64_t offset = 0;
      // Find a non-reverse duplicate sequence without 'N'.
      while (contains_n(seq, kmer) || reverse_duplicate(seq, kmer)) {
         if (locus >= index->size - kmer) offset++;
         else {
            bwpos_t p = {0,0,index->size-1};
            for (int i = kmer-1; i >=0; i--)
               suffix_extend(translate[(int)seq[i]], p, &p, index);
            offset += p.ep - p.sp + 1;
         }
         if (sa_pos + offset >= job->end) break;
         locus = get_sa(sa_pos + offset, index->sa, index->sa_bits);
         seq = index->genome + locus;
      }
      if (sa_pos + offset >= job->end) break;
      // Translate query sequence.
      for (int i = 0; i < kmer; i++) query[i] = translate[(int)seq[i]];
      // Reset hits.
      hits->pos = 0;
      // Search sequence.
      blocksearch(query, kmer, tau, index, &hits);
      // Count hits.
      uint64_t hit_count = 0;
      fmdpos_t self = {.fp = sa_pos, .rp = 0, .sz = 1};
      for (int i = 0; i < hits->pos; i++) {
         spath_t p = hits->path[i];
         if (p.score == 0)
            self = p.pos;
         hit_count += p.pos.sz;
      }
      // Store count.
      job->counts[self.fp] = min(255,hit_count);
      // Update computed values.
      computed += self.fp + self.sz - sa_pos;
      // Update position.
      sa_pos = self.fp + self.sz;
   }

   // Report to scheduler.
   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += computed - last;
   *(job->kmers) += kmers;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);

   // Free memory.
   free(hits);

   return NULL;
}
