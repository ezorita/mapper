#include "annotate.h"

/*
// DEBUG.
char pc_cur[11];
char * pc_nt = "ACGTN";
*/
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
find_next
(
 fmdpos_t current,
 fmdpos_t * next,
 uint8_t * next_seq,
 int       kmer,
 index_t * index
)
{
   uint64_t ptr = current.fp + current.sz - 1;
   uint64_t sa_beg;
   do {
      ptr++;
      if (ptr >= index->size) {
         next->fp = index->size;
         next->sz = 1;
         return 0;
      }
      sa_beg = get_sa(ptr,index->sa,index->sa_bits);
   } while (sa_beg >= index->size - kmer);
   int lcp = seq_trail(index->genome+get_sa(current.fp,index->sa,index->sa_bits), index->genome+sa_beg, kmer);
   // Find range.
   *next = (fmdpos_t) {0,0,index->size};
   for (int i = 0; i < kmer; i++) {
      next_seq[i] = translate[(int)index->genome[sa_beg+i]];
      *next = extend_fw(next_seq[i],*next,index);
   }
   return lcp;
}


int
annotate
(
 int        kmer,
 int        tau,
 uint8_t  * counts,
 index_t  * index,
 int        threads
)
{
   // Compute thread jobs.
   uint64_t seqs = index->size/threads;
   annjob_t ** job = malloc(threads*sizeof(annjob_t*));
   pthread_mutex_t * mutex = malloc(sizeof(pthread_mutex_t));
   pthread_cond_t  * monitor = malloc(sizeof(pthread_cond_t));
   uint64_t * computed = malloc(sizeof(uint64_t));
   int * done = malloc(sizeof(int));
   *computed = 0;
   *done = 0;
   pthread_mutex_init(mutex,NULL);
   pthread_cond_init(monitor,NULL);
   for (int i = 0; i < threads; i++) {
      // Create job.
      job[i] = malloc(sizeof(annjob_t));
      uint64_t offset = i*seqs + 1;
      // Extend kmer from suffix array.
      fmdpos_t pos = {0,0,index->size};
      uint64_t base = index->size;
      do {
         base = get_sa(offset++, index->sa, index->sa_bits);
      } while (base >= index->size - kmer);
      uint8_t * seq = malloc(kmer);
      for (int j = 0; j < kmer; j++) {
         seq[j] = translate[(int)index->genome[base+j]];
         pos = extend_fw(seq[j],pos,index);
      }
      job[i]->beg = pos;
      job[i]->seq = seq;
      job[i]->mutex = mutex;
      job[i]->monitor = monitor;
      job[i]->computed = computed;
      job[i]->done = done;
      job[i]->index = index;
      job[i]->kmer = kmer;
      job[i]->tau = tau;
      job[i]->counts = counts;
   }
   for (int i = 0; i < threads-1; i++) job[i]->end = job[i+1]->beg.fp;
   job[threads-1]->end = index->size;

   for (int i = 0; i < threads; i++) {
      pthread_t thread;
      pthread_create(&thread,NULL,annotate_mt,job[i]);
      pthread_detach(thread);
      //      annotate_mt(job[i]);
  
   }

   pthread_mutex_lock(mutex);
   while (*done < threads) {
      //      fprintf(stderr, "annotating... %ld/%ld\r", *computed, index->size);
      fprintf(stderr, "annotating... %.2f%%\r", *computed*100.0/index->size);
      pthread_cond_wait(monitor,mutex);
   }
   pthread_mutex_unlock(mutex);

   //   fprintf(stderr, "annotating... %ld/%ld\n", *computed, index->size);
   fprintf(stderr, "annotating... 100.00%%\n");

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
   pebble_t root = (pebble_t){(fmdpos_t){0,0,index->size},0,0xFFFFFFFF};
   int8_t * r = root.row + MAXTAU + 1;
   for (int i = 1; i <= MAXTAU + 1; i++) r[i] = r[-i] = i;
   ppush(pebbles, root); 
   // args.
   arg_t args = {job->seq, kmer, 0, tau, index, pebbles};
   // aux vars.
   fmdpos_t current = job->beg, next;
   uint8_t * next_seq = malloc(kmer);
   int start = 0;
   // Iterate over all kmers.
   while (current.fp != job->end) {
      if ((computed>>20) > (last>>20)) {
         pthread_mutex_lock(job->mutex);
         *(job->computed) += computed - last;
         pthread_cond_signal(job->monitor);
         pthread_mutex_unlock(job->mutex);
         last = computed;
      }
      // Get next seq.
      args.trail = find_next(current, &next, next_seq, kmer, index);
      int hits = 0;
      // Reset pebbles.
      for (int j = start + 1; j <= args.trail; j++) pebbles[j]->pos = 0;

      /*
      // DEBUG.
      for (int j = 0; j < kmer; j++) pc_cur[j] = pc_nt[args.query[j]];
      pc_cur[kmer] = 0;
      fprintf(stdout, "sequence fp:%ld,sz:%ld: %s\n", current.fp, current.sz, pc_cur);
      */

      // Iterate over pebbles.
      for (int p = 0 ; p < pebbles[start]->pos ; p++) {
         // Next pebble.
         pebble_t pebble = pebbles[start]->pebble[p];
         /*
         // DEBUG.
         int score = MAXTAU;
         for (int i = 0; i < 2*MAXTAU+3; i++) if (pebble.row[i] < score) score = pebble.row[i];
         memcpy(pc_cur,index->genome+get_sa(pebble.pos.fp, index->sa, index->sa_bits), pebble.depth);
         pc_cur[pebble.depth] = 0;
         fprintf(stdout, "pebble[score:%d], fp:%ld. sz:%ld\n", score, pebble.pos.fp, pebble.pos.sz);
         */
         // Poucet search.
         hits += poucet_mismatch(pebble, args);
         //hits += poucet(pebble, args);
      }
      // write result.
      uint64_t * sa_values = malloc(current.sz*sizeof(uint64_t));
      get_sa_range(current.fp, current.sz, index->sa, index->sa_bits, sa_values);
      uint8_t value = (uint8_t) min(255,hits);
      for (uint64_t i = 0; i < current.sz; i++) {
         job->counts[sa_values[i]] = value;
      } 
      computed += current.sz;
      free(sa_values);
      // swap query sequences.
      uint8_t * tmp = args.query;
      args.query = next_seq;
      next_seq = tmp;
      // update bw interval and start depth.
      current = next;
      start = args.trail;
   }

   pthread_mutex_lock(job->mutex);
   *(job->done) += 1;
   *(job->computed) += computed - last;
   pthread_cond_signal(job->monitor);
   pthread_mutex_unlock(job->mutex);
   
   free(args.query);
   free(next_seq);
   for (int i = 0; i < kmer; i++) free(pebbles[i]);

   return NULL;
}

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
      /*
      // DEBUG.
      pc_cur[depth] = pc_nt[nt];
      pc_cur[depth+1] = 0;
      fprintf(stdout, "fp:%ld,sz:%ld,h:%d\t%s",newpos[nt].fp,newpos[nt].sz,hits,pc_cur);
      */
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
      // Just add the interval size.
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
poucet_mismatch
(
 pebble_t  pebble,
 arg_t     arg
)
{
   int8_t prow = pebble.row[MAXTAU + 1];
   int depth = pebble.depth;
   int tau = arg.tau;

   // Penalty for match/mismatch and insertion/deletion resepectively.
   int hits = 0;

   // Compute new index intervals.
   fmdpos_t newpos[NUM_BASES];
   extend_fw_all(pebble.pos, newpos, arg.index);

   for (int nt = 0 ; nt < NUM_BASES ; nt++) {
      // Check whether child 'i' exists.
      if (newpos[nt].sz < 1) continue;

      // Center cell.
      int score = prow + (nt != arg.query[depth]);

      // Stop searching if 'tau' is exceeded.
      if (score > tau) continue;

      // Reached height of the trie: it's a hit!
      // Just add the interval size.
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
      newpebble.row[MAXTAU+1] = score;

      
      // Cache nodes in pebbles when trailing.
      if (depth < arg.trail)
         ppush(arg.pebbles + depth + 1, newpebble);

      // Dash path if mismatches exhausted.
      if (depth >= arg.trail && score == tau) {
         hits += dash_mismatch(newpebble, arg);
         continue;
      }

      // Recursive call.
      hits += poucet_mismatch(newpebble, arg);
   }

   return hits;
}

int
dash_mismatch
(
 pebble_t    pebble,
 const arg_t arg
)
{
   int hits = 0;
   fmdpos_t pos = pebble.pos;
   for (int i = pebble.depth; i < arg.kmer; i++) {
      int nt = arg.query[i];
      pos = extend_fw(nt,pos,arg.index);
      if (pos.sz < 1) break;
   }
   if (pos.sz > 0) hits += pos.sz;
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
