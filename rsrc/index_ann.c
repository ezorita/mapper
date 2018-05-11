#include "index_ann.h"

// Define interface structures.

struct ann_t {
   size_t      mmap_len;
   void      * mmap_ptr;
   int64_t     kmer;
   int64_t     tau;
   int64_t     size;
   uint8_t   * info;
};


// Private structures.

typedef struct annjob_t annjob_t;

struct annjob_t {
   uint32_t          kmer;
   uint16_t          tau;
   uint16_t          wsize;
   uint64_t          beg;
   uint64_t          end;
   uint64_t        * computed;
   int             * done;
   pthread_mutex_t * mutex;
   pthread_cond_t  * monitor;
   bwt_t           * bwt;
   sar_t           * sar;
   uint8_t         * info;
};


// Private function headers.

void  * ann_build_mt    (void *);
void    jobs_by_prefix  (bwtquery_t * q, int tau, int n_cnt, int depth, int max_depth, int * jobcnt, annjob_t * jobs);
void    neigh_push      (uint8_t * info, uint8_t * ann, int alnsize, int reverse_kmer);
int     neigh_next      (int32_t, int32_t, int64_t, int64_t *, int32_t *, uint8_t *, uint8_t *, bwtquery_t **, sar_t *);
void    aln_merge       (uint8_t * a, uint8_t * b, int len);
void    aln_positions   (uint64_t * bits, uint8_t * pos, int nbits, int npos, int reverse);
void    hits_push       (annjob_t * job, pathstack_t * stack, bwtquery_t * q);
void    update_progress (int64_t sa_pos, int64_t * computed, annjob_t * job);

// Private macros.
#define max(x,y)  ((x) > (y) ? (x) : (y))
#define min(x,y)  ((x) < (y) ? (x) : (y))


/*
** Interface functions source.
*/

ann_t *
ann_build
(
 int        kmer,
 int        tau,
 bwt_t    * bwt,
 sar_t    * sar,
 int        threads
)
{
   if (kmer < 2 || tau >= 4 || tau >= kmer || tau < 1 || threads < 1)
      return NULL;

   if (bwt == NULL || sar == NULL)
      return NULL;

   // Constant parameters.
   int job_to_thread_ratio = 5;
   int           word_size = max(3,tau) + 3;

   // Progress vars.
   uint64_t computed = 0;
   int32_t  done = 0;

   // Memory that will be alloc'ed.
   pthread_mutex_t * mutex    = NULL;
   pthread_cond_t  * monitor  = NULL;
   uint8_t         * tmp_info = NULL;
   annjob_t        * jobs     = NULL;
   bwtquery_t      * q        = NULL;
   ann_t           * ann      = NULL;
   int64_t         * range    = NULL;

   // Alloc variables.
   mutex = malloc(sizeof(pthread_mutex_t));
   monitor = malloc(sizeof(pthread_cond_t));
   
   if (mutex == NULL || monitor == NULL)
      goto failure_return;

   // Index information.
   txt_t   * txt      = bwt_get_text(bwt);
   int64_t   tlen     = txt_length(txt);
   int       num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Temporary info array.
   // Has 'tlen' elements with the following structure:
   //   bytes 0-1: number of different neighbors (16bit).
   //   byte    2: distance of closest neighbor.
   //   byte  3-?: M bytes encoding the mutated positions to find the closest neighbors.
   //   
   //   M is max(3,tau). If the number of mutations is > M, the alignment is not stored.

   tmp_info = calloc(tlen, word_size*sizeof(uint8_t));
   if (tmp_info == NULL)
      goto failure_return;

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
   int num_jobs = num_symb;
   while (num_jobs < job_to_thread_ratio*threads) {
      prefix_depth++;
      num_jobs *= num_symb;
   }

   // Create jobs based on prefixes of depth 'prefix_depth'.
   jobs = malloc(num_jobs*sizeof(annjob_t));
   if (jobs == NULL) 
      goto failure_return;

   num_jobs = 0;
   // Fill specific job info.
   q = bwt_new_query(bwt);
   if (q == NULL)
      goto failure_return;

   jobs_by_prefix(q, tau, 0, 0, prefix_depth, &num_jobs, jobs);

   // Fill job info.
   for (int i = 0; i < num_jobs; i++) {
      jobs[i].kmer     = (uint32_t) kmer;
      jobs[i].tau      = (uint16_t) tau;
      jobs[i].wsize    = (uint16_t) word_size;
      jobs[i].computed = &computed;
      jobs[i].done     = &done;
      jobs[i].mutex    = mutex;
      jobs[i].monitor  = monitor;
      jobs[i].bwt      = bwt;
      jobs[i].sar      = sar;
      jobs[i].info     = tmp_info;
   }

   clock_t clk = clock();
   // Run threads.
   for (int i = 0; i < min(threads, num_jobs); i++) {
      pthread_t thread;
      pthread_create(&thread, NULL, ann_build_mt, jobs+i);
      pthread_detach(thread);
   }

   // Sleep and wait for thread signals.
   int runcount = threads;
   pthread_mutex_lock(mutex);
   while (done < num_jobs) {
      // Run new threads.
      for (int i = runcount; i < min(done + threads, num_jobs); i++) {
         pthread_t thread;
         pthread_create(&thread, NULL, ann_build_mt, jobs+i);
         pthread_detach(thread);
      }
      // Update submitted jobs.
      runcount = done + threads;
      // Report progress.
      fprintf(stderr, "\r[proc] computing neighbors... %.2f%%", computed*100.0/tlen);
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
   ann = malloc(sizeof(ann_t));
   if (ann == NULL)
      goto failure_return;

   ann->kmer = kmer;
   ann->tau  = tau;
   ann->size = tlen/2;
   ann->info = calloc(ann->size, sizeof(uint8_t));

   
   if (ann->info == NULL)
      goto failure_return;

   // Compute and store annotation from info.
   uint64_t i = 0;
   uint8_t * info = tmp_info;
   while (1) {
      // info structure:
      //   bytes     0-1: number of different neighbors (16bit).
      //   byte        2: distance to closest neighbor.
      //   byte  3-wsize: bytes encoding the mutated positions to find the closest neighbors.

      while (i < tlen && (*(uint16_t *)info == ANN_NO_INFO || *(uint16_t *)info == 0)) {
         i++;
         info += word_size;
      }
      if (i >= tlen) break;

      // Compute k-mer range in SA.
      uint64_t size = 1;
      while (i+size < tlen && *(uint16_t *)(info + size*word_size) == 0)
         size++;
      // Get SA values.
      range = malloc(size*sizeof(int64_t));
      if (range == NULL)
         goto failure_return;

      sar_get_range(i, size, range, sar);
      
      // Store info on all SA occurrences.
      for (uint64_t j = 0; j < size; j++) {
         // Store annotation in forward strand.
         if (range[j] >= tlen/2)
            // Convert reverse to fw strand.
            neigh_push(info, ann->info + tlen - range[j] - kmer, word_size-3, kmer);
         else
            neigh_push(info, ann->info + range[j], word_size-3, 0);
      }

      // Free temp memory.
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
   free(q);
   return ann;

 failure_return:
   free(range);
   free(jobs);
   free(mutex);
   free(monitor);
   free(tmp_info);
   free(q);
   ann_free(ann);
   return NULL;
}

void
ann_free
(
  ann_t * ann
)
{
   if (ann != NULL) {
      if (ann->mmap_ptr == NULL) {
         free(ann->info);
      } else {
         munmap(ann->mmap_ptr, ann->mmap_len);
      }
      free(ann);
   }

   return;
}


locinfo_t *
ann_query
(
  int64_t    pos,
  ann_t    * ann
)
{
   if (ann == NULL)
      return NULL;

   if (pos < 0 || pos >= ann->size*2)
      return NULL;

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

   int strand = 0;
   if (pos > ann->size) {
      strand = 1;
      pos    = (ann->size-1)*2 - pos;
   }
   
   uint8_t   info    = ann->info[pos];
   int32_t   aln_cnt = 0;
   int32_t * aln_pos = malloc(ann->kmer*sizeof(int32_t));
   if ((info >> 6) & 1) {
      for (int i = 0; i < ann->kmer; i++) {
         if ((ann->info[pos+i] >> 7) & 1) {
            aln_pos[aln_cnt++] = (strand ? ann->kmer - 1 - i : i);
         }
      }
   }
   
   locinfo_t * locinfo = malloc(sizeof(locinfo_t) + aln_cnt*sizeof(int32_t));
   if (locinfo == NULL)
      return NULL;

   int cnt = (info & 0x0F);
   
   locinfo->dist      = (cnt ? ((info >> 4) & 3) + 1 : 0);
   locinfo->align_cnt = aln_cnt;
   memcpy(&(locinfo->align_pos[0]), aln_pos, aln_cnt*sizeof(int32_t));

   // Convert neigh_cnt.
   if      (cnt <= 10)   locinfo->neigh_cnt = cnt;
   else if (cnt == 0x0B) locinfo->neigh_cnt = 15;
   else if (cnt == 0x0C) locinfo->neigh_cnt = 40;
   else if (cnt == 0x0D) locinfo->neigh_cnt = 75;
   else if (cnt == 0x0E) locinfo->neigh_cnt = 300;
   else if (cnt == 0x0F) locinfo->neigh_cnt = 1000;

   free(aln_pos);

   return locinfo;
}


int
ann_file_write
(
  char   * filename,
  ann_t  * ann
)
{
   // Check arguments.
   if (filename == NULL || ann == NULL)
      return -1;
   
   // Open file.
   int fd = creat(filename, 0644);
   if (fd == -1)
      return -1;

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = ANN_FILE_MAGICNO;

   // Write magic.
   if (write(fd, &magic, sizeof(uint64_t)) == -1)
      goto close_and_error;
   
   // Write kmer.
   if (write(fd, (int64_t *)&(ann->kmer), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write tau.
   if (write(fd, (int64_t *)&(ann->tau), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write size.
   if (write(fd, (int64_t *)&(ann->size), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write info array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint8_t *)ann->info + e_cnt, (ann->size - e_cnt)*sizeof(uint8_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < ann->size);

   close(fd);
   return 0;
   
 close_and_error:
   close(fd);
   return -1;
}


ann_t *
ann_file_read
(
  char * filename
)
{
   // Check arguments.
   if (filename == NULL)
      return NULL;

   // Open file.
   int fd = open(filename, O_RDONLY);
   if (fd == -1)
      return NULL;

   // Alloc memory.
   ann_t * ann = malloc(sizeof(ann_t));
   if (ann == NULL)
      goto free_and_return;

   // Set NULL pointers.
   ann->info = NULL;

   // Get file len and mmap file.
   struct stat sb;
   fstat(fd, &sb);
   if (sb.st_size < 48)
      goto free_and_return;

   int64_t * data = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
   if (data == NULL)
      goto free_and_return;
   
   ann->mmap_len = sb.st_size;
   ann->mmap_ptr = (void *) data;

   // Read file.
   // Read magic number.
   uint64_t magic = data[0];
   if (magic != ANN_FILE_MAGICNO)
      goto free_and_return;

   ann->kmer = data[1];
   ann->tau  = data[2];
   ann->size = data[3];
   ann->info = (uint8_t *) (data + 4);

   close(fd);

   return ann;

 free_and_return:
   close(fd);
   ann_free(ann);
   return NULL;
}

/*
** Private functions source.
*/

void *
ann_build_mt
(
 void * argp
 )
{
   // Read job.
   annjob_t * job      = (annjob_t *)argp;
   bwt_t    * bwt      = job->bwt;
   sar_t    * sar      = job->sar;
   int64_t    sa_ptr   = job->beg;
   int        kmer     = job->kmer;
   int        tau      = job->tau;
   int        num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));

   // Translated query.
   uint8_t * query_1 = malloc(kmer*sizeof(uint8_t));
   uint8_t * query_2 = malloc(kmer*sizeof(uint8_t));
   memset(query_1, num_symb, kmer*sizeof(uint8_t));
   memset(query_2, num_symb, kmer*sizeof(uint8_t));
   // BWT path.
   bwtquery_t ** path = malloc((kmer+1)*sizeof(bwtquery_t *));
   for (int i = 0; i <= kmer; i++) {
      path[i] = bwt_new_query(bwt);
   }

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
      update_progress(sa_ptr, &computed, job);
      // This updates query and path, and makes sa_ptr point to the position of the next sequence in the SA.
      int trail = neigh_next(kmer, tau, sa_ptr, &next_sa, &last_fragment, query_1, query_2, path, sar);
      // Check valid query sequence.
      if (trail < 0) {
         // Flag 'sa_ptr' as ANN_NO_INFO.
         *(uint16_t *)(job->info + sa_ptr * job->wsize) = ANN_NO_INFO;
      } else {
         // Query the sequence.
         if (last_fragment)
            blocksc_trail(query_2, path, kmer, tau, trail, stack_tree_2);
         else
            blocksc_trail(query_1, path, kmer, tau, trail, stack_tree_1);
         // Update annotation info.
         hits_push(job, (last_fragment ? stack_tree_2->stack : stack_tree_1->stack), path[kmer]);
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
   for (int i = 0; i <= kmer; i++) {
      free(path[i]);
   }


   return NULL;

}


void
jobs_by_prefix
(
 bwtquery_t  * q,
 int           tau,
 int           n_cnt,
 int           depth,
 int           max_depth,
 int         * jobcnt,
 annjob_t    * jobs
)
{
   // Too many 'N'.
   if (n_cnt > tau) 
      return;
   // If reached suffix depth, make and store job.
   if (depth == max_depth) {
      annjob_t job = {
         .beg = bwt_start(q),
         .end = bwt_start(q) + bwt_size(q)
      };
      jobs[(*jobcnt)++] = job;
      return;
   }
   
   // Otherwise continue exploring the suffix trie.
   bwt_t       * bwt = bwt_get_bwt(q);
   bwtquery_t ** qv  = bwt_new_vec(bwt);
   
   bwt_query_all(BWT_QUERY_SUFFIX, q, qv);

   // Iterate over all symbols.
   int num_symb = sym_count(txt_get_symbols(bwt_get_text(bwt)));
   for (int i = 0; i < num_symb; i++) {
      jobs_by_prefix(qv[i], tau, n_cnt + (i == UNKNOWN_BASE), depth + 1, max_depth, jobcnt, jobs);
   }
   bwt_free_vec(qv);

   return;
}


void
neigh_push
(
 uint8_t * info,
 uint8_t * ann,
 int       alnsize,
 int       reverse_kmer
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

   // Return if kmer has perfect match or no info.
   if (cnt == 0 || cnt == ANN_NO_INFO) return;

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
      // Copy alignment data.
      aln = malloc(alnsize);
      memcpy(aln, info + sizeof(uint16_t) + sizeof(uint8_t), alnsize);

      // Update alignments from reverse strand.
      if (reverse_kmer)
         for (int k = 0; k < alnsize && aln[k]; k++)
            aln[k] = reverse_kmer + 1 - aln[k];

      // Set alignment flag in locus.
      *ann |= (uint8_t)1 << 6;
      // Store alignment info.
      for (int i = 0; i < alnsize && aln[i]; i++)
         // Alignment values are 1-based.
         *(ann + aln[i]-1) |= (uint8_t)1 << 7;

      free(aln);
   }

   return;
}


int
neigh_next
(
 int32_t       qlen,
 int32_t       tau,
 int64_t       sa_ptr,
 int64_t     * next_sa,
 int32_t     * last_fragment,
 uint8_t     * query_1,
 uint8_t     * query_2,
 bwtquery_t ** q,
 sar_t       * sar
)
{
   int n_cnt = 0;
   int trail_1 = 0, trail_2 = 0, truncated = 0;

   // Index data.
   bwt_t * bwt = bwt_get_bwt(q[qlen]);
   txt_t * txt = bwt_get_text(bwt);
   int     num_symb = sym_count(txt_get_symbols(txt));

   // Get suffix array.
   int64_t txt_pos = sar_get(sa_ptr, sar);

   // Increase next_sa by one (error default).
   *next_sa = sa_ptr + 1;

   // Check text limits.
   if (txt_pos + qlen > txt_length(txt)) {
      return -1;
   }

   // Get symbols.
   uint8_t * seq = txt_sym_range(txt_pos, qlen, txt);
   
   // Report sequences containing '$'.
   for (int i = 0; i < qlen; i++) {
      if (seq[i] >= num_symb) {
         return -1;
      }
   }

   // Iterate over qlen nucleotides.
   uint8_t * tmp = malloc(qlen*sizeof(uint8_t));
   for (int i = 0; i < qlen && !truncated; i++) {
      n_cnt += seq[i] == UNKNOWN_BASE;
      // Track trail.
      if (trail_1 == i && seq[i] == query_1[i])
         trail_1++;
      if (trail_2 == i && seq[i] == query_2[i])
         trail_2++;
      // Extend BWT search.
      bwt_query(seq[i], BWT_QUERY_SUFFIX, q[i], q[i+1]);

      // Update tmp query.
      tmp[i] = seq[i];
   }

   // Check whether sequence is valid (all in same strand).
   if (bwt_size(q[qlen]) == 0 || truncated) {
      free(tmp);
      return -1;
   }

   // Update sa_ptr (point the start of the next suffix).
   *next_sa = sa_ptr + bwt_size(q[qlen]);

   // Check N count.
   if (n_cnt > tau) {
      free(tmp);
      return -1;
   }

   // Copy tmp to query.
   if (bwt_start(q[qlen]) >= bwt_rcstart(q[qlen])) {
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
aln_merge
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
aln_positions
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
hits_push
(
 annjob_t     * job,
 pathstack_t  * stack,
 bwtquery_t   * q
)
// Stores new information discovered in the last query.
// Information is always stored/updated in the lexicographically
// smallest between the queried/hit sequence and its rev complement.
{
   int hits     = 0;
   int tau      = job->tau + 1;
   int aln_size = job->wsize - 3;

   int64_t fp = bwt_start(q);
   int64_t rp = bwt_rcstart(q);

   // Flag lexicographically greater as ANN_NO_INFO.
   if (fp > rp) {
      pthread_mutex_lock(job->mutex);
      *(uint16_t *)(job->info + fp * job->wsize) = ANN_NO_INFO;
      pthread_mutex_unlock(job->mutex);
   }

   // No matches (self hit only).
   if (stack->pos < 2) {
      pthread_mutex_lock(job->mutex);
      // If there is no previous info, set smaller as ANN_NO_INFO.
      if (*(uint16_t *)(job->info + min(fp, rp) * (size_t)(job->wsize)) == 0)
         *(uint16_t *)(job->info + min(fp, rp) * (size_t)(job->wsize)) = ANN_NO_INFO;
      pthread_mutex_unlock(job->mutex);
      return;
   }

   // Aggregate query alignment.
   uint64_t * qalign = calloc(ALIGN_WORDS, sizeof(uint64_t));

   // Store hits info in their SA start position.
   for (int i = 0; i < stack->pos; i++) {
      spath_t path = stack->path[i];
      int64_t path_fp = bwt_start(path.bwtq);
      int64_t path_rp = bwt_rcstart(path.bwtq);
      // Remove self hit.
      if (path_fp == fp) continue;

      // Select lexicographically smallest between the neighbor and its revcomp.
      int64_t nptr = (path_rp < path_fp ? path_rp : path_fp);
      int     nrev = path_rp < path_fp;

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
            aln_positions(path.align, hit_aln, job->kmer, aln_size, nrev);
            // Merge mismatch positions and store in neighbor memory (naln).
            aln_merge(naln, hit_aln, aln_size);
         }
      }
      // Reset neighbor.
      else if (*ncnt == 0 || *ncnt == ANN_NO_INFO || *ntau > path.score) {
         // Set count to 1 and neighbor's tau to path.score.
         *ncnt = 1;
         *ntau = path.score;
         // Compute mismatch positions.
         uint8_t hit_aln[aln_size];
         aln_positions(path.align, hit_aln, job->kmer, aln_size, nrev);
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
   int64_t qptr = (rp < fp ? rp : fp);
   int     qrev = rp < fp;

   // Query info.
   uint8_t  * info = job->info + qptr * job->wsize;
   uint16_t * qcnt = (uint16_t *) info;
   uint8_t  * qtau = info + sizeof(uint16_t);
   uint8_t  * qaln = qtau + 1;

   // If query is empty or tau < qtau, reset query.
   pthread_mutex_lock(job->mutex);
   if (*qcnt == 0 || *qcnt == ANN_NO_INFO || *qtau > tau) {
      // Set tau and neighbor count.
      *qtau = (uint8_t)  tau;
      *qcnt = (uint16_t) min(ANN_NO_INFO - 1, hits);
      // Compute and directly store mismatch positions for query sequence.
      aln_positions(qalign, qaln, job->kmer, aln_size, qrev);
   }
   // Update query.
   else if (*qtau == tau) {
      // Update hit count.
      *qcnt = (uint16_t) min(ANN_NO_INFO - 1, ((uint32_t)hits)+((uint32_t)*qcnt));
      // Check if query has free align slots.
      if (*qaln != 0xFF) {
         // Compute mismatch positions for query sequence.
         uint8_t tmp_aln[aln_size];
         aln_positions(qalign, tmp_aln, job->kmer, aln_size, qrev);
         // Update alignments in query by merging.
         aln_merge(qaln, tmp_aln, aln_size);
      }
   }
   pthread_mutex_unlock(job->mutex);

   free(qalign);
}

void
update_progress
(
 int64_t     sa_pos,
 int64_t   * computed,
 annjob_t  * job
)
{
   // Send progress every 1M processed sequences.
   if (__builtin_expect(( (*computed) >> 20) > (sa_pos>>20) , 0)) {
      pthread_mutex_lock(job->mutex);
      *(job->computed) += sa_pos - *computed;
      pthread_cond_signal(job->monitor);
      pthread_mutex_unlock(job->mutex);
      *computed = sa_pos;
   }
   return;
}
