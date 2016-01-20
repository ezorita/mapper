#include "annotate.h"
#include <string.h>

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
ppush
(
 pstack_t ** stackp,
 pebble_t    pebble
 )
{
   pstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      size_t newsize = 2*stack->size;
      stack = *stackp = realloc(stack, sizeof(pstack_t) + newsize*sizeof(pebble_t));
      if (stack == NULL) return -1;
      stack->size = newsize;
   }
   stack->pebble[stack->pos++] = pebble;

   return 0;
}

pstack_t *
pstack_new
(
 size_t size
 )
{
   if (size < 1) size = 1;
   pstack_t * stack = malloc(sizeof(pstack_t) + size*sizeof(pebble_t));
   if (stack == NULL) return NULL;
   stack->pos = 0;
   stack->size = size;
   return stack;
}

int
seq_trail
(
 char * seq_a,
 char * seq_b,
 int    kmer
 )
{
   int lcp = 0;
   while (translate[(int)seq_a[lcp]] == translate[(int)seq_b[lcp]] && lcp < kmer) lcp++;
   return lcp;
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
   for (int i = 0; i < kmer; i++) {
      if (seq[i] == 'n' || seq[i] == 'N') return 1;
   }
   return 0;
}


int
find_next
(
 bwpos_t  * pos,
 uint8_t  * seq,
 char     * locus,
 int        kmer,
 bwpos_t  * pos_cache,
 index_t  * index
 )
{
   int lcp = 0;
   int revdup = 0;
   do {
      uint64_t ptr = pos->ep;
      uint64_t sa_beg;
      // Get next valid interval: SA value must be < genome_size - kmer.
      // TODO: Also check that the kmer does not fall in the border between
      // forward and reverse strands of the genome, i.e. the SA value must not
      // fall in the region (genome->size/2-kmer , genome->size/2].
      do {
         ptr++;
         if (ptr >= index->size) {
            pos->sp = index->size;
            pos->ep = index->size;
            return 0;
         }
         sa_beg = get_sa(ptr,index->sa,index->sa_bits);
      } while (sa_beg >= index->size - kmer);

      // Keep pointer to next kmer (genome sequence).
      char * next_locus = index->genome+sa_beg;

      // Compute LCP between current and next intervals.
      lcp = seq_trail(locus, next_locus, kmer);

      // Compute next interval range using pos cache.
      // ftr: Forward to reverse strand conversion table.
      int ftr[5] = {3,2,1,0,4};
      // Start at LCP.
      for (int i = lcp; i < kmer ; i++) {
         // Read the sequence backwards.
         int nt = translate[(int)next_locus[i]];
         // Store the reverse complement in seq.
         seq[i] = ftr[nt];
         // Update cache with backward search of new nucleotides (using reverse complement).
         suffix_extend(seq[i],pos_cache[i],pos_cache+i+1,index);
      }

      // Compute the interval of the new kmer using the revcomp interval size.
      *pos = (bwpos_t){.sp = ptr, .ep = ptr + pos_cache[kmer].ep - pos_cache[kmer].sp, .depth = kmer};
      revdup = reverse_duplicate(next_locus, kmer) | contains_n(next_locus, kmer);
   } while (revdup);
   // Repeat the whole process if the kmer is a duplicate.
   
   return lcp;
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
      // Find a non-duplicated sequence without 'N'.
      int duplicate = 0;
      uint8_t * seq = malloc(kmer);
      bwpos_t * cache = malloc((kmer+1) * sizeof(bwpos_t));
      do {
         // Extend kmer from suffix array.
         uint64_t base;
         do {
            base = get_sa(++offset, index->sa, index->sa_bits);
         } while (base >= index->size - kmer);

         // Alloc interval cache.
         cache[0] = (bwpos_t){0,0,index->size-1};;
         int ftr[5] = {3,2,1,0,4};
         for (int j = 0; j < kmer; j++) {
            // Backward search the reverse complement to reuse the cached intervals.
            int nt = translate[(int)index->genome[base+j]];
            seq[j] = ftr[nt];
            suffix_extend(seq[j],cache[j],cache+j+1,index);
         }
         duplicate = reverse_duplicate(index->genome+base, kmer) | contains_n(index->genome+base, kmer);
      } while (duplicate && offset < eff_size);

      // Break if eff_size is covered.
      if (offset >= eff_size) {
         njobs = i;
         break;
      }

      // Compute pos using the forward start pointer and revcomp interval size.
      bwpos_t pos = {.sp = offset, .ep = offset + cache[kmer].ep - cache[kmer].sp, .depth = kmer};
      
      job[i]->beg = pos;
      job[i]->seq = seq;
      job[i]->mutex = mutex;
      job[i]->monitor = monitor;
      job[i]->computed = computed;
      job[i]->done = done;
      job[i]->cache = cache;
      job[i]->index = index;
      job[i]->kmer = kmer;
      job[i]->tau = tau;
      job[i]->kmers = kmers;
      job[i]->counts = counts;
   }
   // Set job end indices.
   for (int i = 0; i < njobs-1; i++) job[i]->end = job[i+1]->beg.sp;
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

   for (int i = 0; i < threads; i++) { free(job[i]->seq); free(job[i]->cache); free(job[i]); }
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
   // pebble stack.
   pstack_t ** pebbles = malloc(kmer*sizeof(pstack_t*));
   for (int i = 0; i < kmer; i++) pebbles[i] = pstack_new(64);
   pstack_t * hits = pstack_new(64);
   // hit stack.

#if COMPUTE_INDELS == 1
   pebble_t root = (pebble_t){(bwpos_t){0,0,index->size-1},0,0xFFFFFFFF};
   int8_t * r = root.row + MAXTAU + 1;
   for (int i = 1; i <= MAXTAU + 1; i++) r[i] = r[-i] = i;
#else
   pebble_t root = (pebble_t){(bwpos_t){0,0,index->size-1},0,0};
#endif
   ppush(pebbles, root); 
   // args.
   arg_t args = {job->seq, kmer, 0, tau, index, pebbles, &hits};
   // aux vars.
   bwpos_t current, pos = job->beg;
   uint8_t * seq = malloc(kmer);
   char * locus = index->genome + get_sa(pos.sp, index->sa, index->sa_bits);
   uint64_t kmers = 0;
   int start = 0;

   // Iterate over all kmers.
   while (pos.sp < job->end) {
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

      // Get next seq.
      current = pos;
      memcpy(seq, args.query, kmer);
      // DEBUG.
      args.trail = find_next(&pos, seq, locus, kmer, job->cache, index);
      // DISABLE POUCET SEARCH.
      /*
      args.trail = 0;
      find_next(&pos, seq, locus, kmer, job->cache, index);
      */
      // Reset hits and pebbles.
      for (int j = start + 1; j <= args.trail; j++) pebbles[j]->pos = 0;

      // Iterate over pebbles.
      hits->pos = 0;
      for (int p = 0 ; p < pebbles[start]->pos ; p++) {
         // Next pebble.
         pebble_t pebble = pebbles[start]->pebble[p];
         // Poucet search.
#if COMPUTE_INDELS == 1
         poucet(pebble, args);         
#else
         poucet_mismatch(pebble, args);
#endif
      }
      uint64_t current_loci = current.ep - current.sp + 1;
      uint64_t hit_count = 0;
      for (int h = 0; h < hits->pos; h++) {
         job->counts[hits->pebble[h].pos.sp] += current_loci;
         hit_count += hits->pebble[h].pos.ep - hits->pebble[h].pos.sp + 1;
      }
      // Add the remaining hit count to current loci.
      job->counts[current.sp] = (uint8_t)min(255,job->counts[current.sp] + hit_count + current_loci);

      // Update computed values.
      computed += pos.sp - current.sp;
      // update query sequence.
      memcpy(args.query, seq, kmer);
      // update bw interval and start depth.
      start = args.trail;
   }

   // Report to scheduler.
   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += computed - last;
   *(job->kmers) += kmers;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);
   
   free(seq);
   free(hits);
   for (int i = 0; i < kmer; i++) free(pebbles[i]);

   return NULL;
}

#if COMPUTE_INDELS == 1
int
poucet
(
 pebble_t  pebble,
 arg_t     arg
 )
{
   int8_t * prow = pebble.row + MAXTAU + 1;
   int depth = pebble.depth;
   int tau = arg.tau;

   // Penalty for match/mismatch and insertion/deletion resepectively.
   int mmatch;
   int shift;
   int hits = 0;

   // Part of the cache that is shared between all the children.
   char r[2*MAXTAU+3];
   for (int i = 0; i < MAXTAU + 2; i++)
      r[i] = r[2*MAXTAU+2-i] = MAXTAU + 1 - i;
   char * row = r + MAXTAU + 1;

   // Upper arm of the L (need the path).
   for (int a = tau ; a > 0 ; a--) {
      mmatch = prow[a] + (((pebble.path >> (3*a)) & 0x07) != arg.query[depth]);
      shift = min(prow[a-1], row[a+1]) + 1;
      row[a] = min(mmatch, shift);
   }

   // Compute new index intervals.
   fmdpos_t newpos[NUM_BASES];
   extend_fw_all(pebble.pos, newpos, arg.index);

   for (int nt = 0 ; nt < NUM_BASES ; nt++) {
      // Check whether child 'i' exists.
      if (newpos[nt].sz < 1) {
         // DEBUG.
         //fprintf(stdout,"\n");
         continue;
      }

      // Horizontal arm of the L (need previous characters).
      for (int i = tau ; i > 0 ; i--) {
         mmatch = prow[-i] + (depth-i < 0 ? 1 : (nt != arg.query[depth-i]));
         shift = min(prow[1-i], row[-i-1]) + 1;
         row[-i] = min(mmatch, shift);
      }

      // Center cell (need both arms to be computed).
      mmatch = prow[0] + (nt != arg.query[depth]);
      shift = min(row[-1], row[1]) + 1;
      row[0] = min(mmatch, shift);

      int score = tau+1;
      for (int i = -tau ; i <= tau; i++)
         if (row[i] <= score)
            score = row[i];
      // DEBUG.
      //fprintf(stdout,"\tscore:%d\t%s\n",score,(depth+1==arg.kmer && score <= tau ? "***" : ""));

      // Stop searching if 'tau' is exceeded.
      if (score > tau) continue;

      // Reached height of the trie: it's a hit!
      // Add the interval size.
      if (depth + 1 == arg.kmer) {
         hits += newpos[nt].sz;
         continue;
      }

      // Update pebble.
      pebble_t newpebble = {
         .pos  = newpos[nt],
         .depth = depth + 1,
         .path = (pebble.path << 3) | nt
      };
      memcpy(&(newpebble.row[0]), r, 2*MAXTAU+3);

      
      // Cache nodes in pebbles when trailing.
      if (depth < arg.trail)
         ppush(arg.pebbles + depth + 1, newpebble);

      // Dash path if mismatches exhausted.
      if (depth >= arg.trail && score == tau) {
         hits += dash(newpebble, arg);
         continue;
      }

      // Recursive call.
      hits += poucet(newpebble, arg);
   }

   return hits;
}


int
dash
(
 pebble_t    pebble,
 const arg_t arg
 )
{
   int hits = 0;
   int8_t * r = pebble.row +MAXTAU+1;
   // DEBUG.
   //   fprintf(stdout,"dash\n");
   for (int j = -arg.tau; j <= arg.tau; j++) {
      /*
      // DEBUG.
      fprintf(stdout,"r[%d]:%d\n",j,r[j]);
      */
      if (r[j] == arg.tau) {
         int i = pebble.depth;
         fmdpos_t pos = pebble.pos;
         if (j > 0) {
            int mism = 0;
            for (int k = j - 1; k >= 0 && i < arg.kmer ; k--) {
               if (((pebble.path >> (3*k)) & 0x07) != arg.query[i++]) {
                  mism = 1;
                  break;
               }
            }
            if (mism) continue;
         } else {
            i += j;
         }
         for (; i < arg.kmer; i++) {
            int nt = arg.query[i];
            pos = extend_fw(nt, pos, arg.index);
            /*
            // DEBUG.
            fprintf(stdout,"Q:%c,fp:%ld,sz:%ld\n",pc_nt[nt],pos.fp,pos.sz);
            */
            if (pos.sz < 1) break;
         }
         /*
         // DEBUG.
         fprintf(stdout, "fp:%ld,sz:%ld\t%s\t%s\n",pos.fp,pos.sz,pc_cur,pos.sz > 0 ? "***" : "");
         */

         if (pos.sz > 0) hits += pos.sz;
      }
   }

   return hits;
}


#else


uint64_t
poucet_mismatch
(
 pebble_t  pebble,
 arg_t     arg
 )
{
   int depth = pebble.depth;
   int tau = arg.tau;

   // Compute new index intervals.
   bwpos_t newpos[NUM_BASES];
   suffix_extend_all(pebble.pos, newpos, arg.index);

   // Since we query the reverse complement of the kmers, we start by the last
   // (alphabetically), hence we don't need to recompute the right parts of the
   // tree since these were queried and reported by the previous kmers.
   for (int nt = (pebble.score ? 0 : arg.query[depth]) ; nt < NUM_BASES; nt++) {
   //   for (int nt = 0 ; nt <= (pebble.score ? 3 : arg.query[depth]) ; nt++) {
   //   for (int nt = 0 ; nt < NUM_BASES; nt++) {
      int8_t score = pebble.score;
      // Check whether child 'i' exists.
      if (newpos[nt].ep < newpos[nt].sp)
         continue;

      // Center cell.
      score += (nt != arg.query[depth]);

      // Stop searching if 'tau' is exceeded.
      if (score > tau) continue;

      // Update pebble.
      pebble_t newpebble = {
         .pos  = newpos[nt],
         .depth = depth + 1,
         .score = score
      };

      // Reached height of the trie: it's a hit!
      // Just add the interval size.
      if (depth + 1 == arg.kmer) {
         if (score) ppush(arg.hits, newpebble);
         continue;
      }
     
      // Cache nodes in pebbles when trailing.
      if (depth < arg.trail)
         ppush(arg.pebbles + depth + 1, newpebble);

      // Dash path if mismatches exhausted.
      if (depth >= arg.trail && score == tau) {
         dash_mismatch(newpebble, arg);
         continue;
      }

      // Recursive call.
      poucet_mismatch(newpebble, arg);
   }

   return 0;
}

uint64_t
dash_mismatch
(
 pebble_t    pebble,
 const arg_t arg
 )
{
   bwpos_t pos = pebble.pos;
   for (int i = pebble.depth; i < arg.kmer; i++) {
      int nt = arg.query[i];
      suffix_extend(nt,pos,&pos,arg.index);
      if (pos.ep < pos.sp) break;
   }
   if (pos.sp <= pos.ep) 
      ppush(arg.hits,(pebble_t){.pos = pos, .depth = arg.kmer, .score = arg.tau});

   return 0;
}
#endif

