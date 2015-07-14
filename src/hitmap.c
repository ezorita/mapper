#include "hitmap.h"

int
hitmap
(
 index_t    * index,
 seqstack_t * seqs,
 hmargs_t     hmargs
 )
{
   int maxtau = 0;
   for (int i = 0; i < hmargs.search_rounds; i++) if (hmargs.tau[i] > maxtau) maxtau = hmargs.tau[i];

   // Define number of seqs that will be mixed in each poucet search.
   int seq_block = hm_min(seqs->pos, hmargs.sequence_blocks);

   pstack_t **pebbles = NULL, **hits = NULL;
   vstack_t **hitmaps = NULL;
   matchlist_t * seeds = NULL;
   // Initialize pebble and hit stacks.
   pebbles = malloc((MAX_TRAIL+1)*sizeof(pstack_t *));
   hits    = malloc((maxtau+1) * sizeof(pstack_t *));
   for (int i = 0 ; i <= MAX_TRAIL ; i++) pebbles[i] = new_pstack(PBSTACK_SIZE);
   for (int i = 0 ; i <= maxtau ; i++)    hits[i]    = new_pstack(HITSTACK_SIZE);

   // Push root node.
   pebble_t root = {.sp = 0, .ep = index->gsize-1};
   for (int i = 0; i < MAXTAU + 2; i++)
      root.nwrow[i] = root.nwrow[2*MAXTAU+2-i] = MAXTAU + 1 - i;

   ppush(pebbles, root);

   // Allocate hitmaps.
   hitmaps = malloc(seq_block * sizeof(vstack_t *));
   for (int i = 0 ; i < seq_block ; i++) hitmaps[i] = new_stack(HITMAP_SIZE);

   // Allocate loci stack.
   seeds = matchlist_new(hmargs.max_align_per_read);

   // Allocate match buffer.
   matchlist_t ** seqmatches;
   seqmatches = malloc(seq_block * sizeof(matchlist_t *));
   for (int i = 0; i < seq_block; i++) seqmatches[i] = matchlist_new(HITBUF_SIZE);

   // Count number of mapped sequences.
   uint64_t mapped = 0;

   // Iterate over sequence blocks.
   for (int s = 0 ; s < seqs->pos; s += seq_block) {
      // Pre-process sequence block: chunk in subsequences and sort.
      int         numseqs = (seq_block < seqs->pos - s ? seq_block : seqs->pos - s);
      sublist_t * subseqs = process_subseq(seqs->seq+s, numseqs, hmargs.kmer_size[0], hmargs.kmer_offset[0], hmargs.qthr[0], hitmaps);

      // Verbose processed blocks.
      if(hmargs.verbose) fprintf(stderr, "processing reads from %d to %d\n", s+1, s+numseqs);

      // Iterate over increasing mismatch tolerance.
      for (int a = 0; a < hmargs.search_rounds; a++) {
         int tau = hmargs.tau[a];
         int kmer = hmargs.kmer_size[a];
         // Verbose sorting.
         if(hmargs.verbose) fprintf(stderr, "[k:%d,d:%d,o:%d,t:%d,q:%d] %d seeds\nsorting  ...", kmer, tau, hmargs.kmer_offset[a], hmargs.seed_thr[a], hmargs.qthr[a],subseqs->pos);
         clock_t tstart = clock();
         // Sort subsequences.
         mergesort_mt(subseqs->sub, subseqs->pos, sizeof(sub_t), kmer, 1, compar_seqsort);
         // Verbose sort time.
         if(hmargs.verbose) fprintf(stderr, "  [%.3fms]\n", (clock()-tstart)*1000.0/CLOCKS_PER_SEC);

         // Clear hitmaps.
         for (int i = 0 ; i < numseqs ; i++) hitmaps[i]->pos = 0;
         // Search subsequences.
         if (tau > 0)
            poucet_search(subseqs, pebbles, hits, index, tau, kmer, kmer, hmargs.seed_max_loci, hmargs.seed_abs_max_loci, hmargs.verbose);
         else 
            bw_search(subseqs, kmer, hmargs.seed_thr[a], index, hmargs);

         // Reset the subsequence list. It will be filled up with the unmatched regions.
         subseqs->pos = 0;

         // Sort all hitmaps
         long recompute = 0, matched = 0, totalnt = 0, num_aligns = 0;
         // Reset verbose variables.
         double progress = 0.0;
         tstart = clock();
         clock_t himap_time = 0;
         clock_t align_time = 0;
         // Process reads.
         for (int i = 0 ; i < numseqs ; i++) {
            // Verbose progress.
            if (hmargs.verbose && i*100.0 / numseqs - progress > 0.1) {
               progress = i*100.0/numseqs;
               fprintf(stderr, "aligning ...  %.1f%%\r", progress);
            }

            // Sequence length.
            int slen = strlen(seqs->seq[i].seq);
            totalnt += slen;

            // Process hitmap.
            clock_t tmp = clock();
            hitmap_analysis(hitmaps[i], seeds, kmer, slen, hmargs);
            /*
            // DEBUG.
            for (int i = 0; i < seeds->pos; i++) {
            match_t * match = seeds->match[i];
            long     g_start = index->gsize - match->ref_e;
            long       g_end = index->gsize - match->ref_s - 1;
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
            fprintf(stdout, "%d\t%s:%ld-%ld\t%d-%d\n", match->hits, chr->name[chrnum],g_start-chr->start[chrnum], g_end-chr->start[chrnum], match->read_s, match->read_e);
            }
            fprintf(stdout, "SEED LIST END\n");
            */

            himap_time += clock() - tmp;
            num_aligns += seeds->pos;
            // Smith-Waterman alignment.
            tmp = clock();
            int newdata = align_seeds(seqs->seq[i], seeds, seqmatches + i, index, hmargs);
            align_time += clock() - tmp;

            if (newdata) {
               // Fuse matches.
               //               fuse_matches(seqmatches+i, slen, hmargs);

               // Find sequence repeats.
               find_repeats(seqmatches[i], hmargs.repeat_min_overlap);
            }
            // Assemble read.
            matchlist_t * intervals = combine_matches(seqmatches[i], hmargs.overlap_tolerance);

            if (a < hmargs.search_rounds - 1) {
               // Feedback gaps between intervals and low identity matches.
               recompute += feedback_gaps(i, seqs->seq[i], intervals, &subseqs, hitmaps + i, hmargs, hmargs.kmer_size[a+1], hmargs.kmer_offset[a+1], hmargs.qthr[a+1]);
            }
            else {
               // Print results if a == maxtau.
               // Fill read gaps allowing greater overlap.
               fill_gaps(&intervals, seqmatches[i], slen, 1 - hmargs.overlap_tolerance, hmargs.overlap_max_tolerance);

               // Print intervals.
               if (intervals->pos) {
                  mapped++;
                  print_intervals(seqs->seq[s+i].tag, intervals, index, hmargs.repeat_print_num);
               } else fprintf(stdout, "%s\t-\n", seqs->seq[s+i].tag);
            }

            // Compute matched nucleotides.
            for (int i = 0; i < intervals->pos; i++) {
               match_t * current = intervals->match[i];
               long next_s;
               if (i < intervals->pos - 1) next_s = intervals->match[i+1]->read_s;
               else next_s = slen;
               matched += hm_min(next_s, current->read_e + 1) - current->read_s;
            }

            free(intervals);
         }
         if (hmargs.verbose) {
            fprintf(stderr, "aligning ...  [%.3fms]\n", (clock()-tstart)*1000.0/CLOCKS_PER_SEC);
            //            fprintf(stderr, "hitmap analysis: [%.3fms], %ld alignments in [%.3fms]\n", himap_time*1000.0/CLOCKS_PER_SEC, num_aligns, align_time*1000.0/CLOCKS_PER_SEC);
            fprintf(stderr, "matched  ...  %.1f%%\n", matched*100.0/totalnt);
         }
      }
   }
   if (hmargs.verbose) fprintf(stderr, "stats: %.2f%% mapped reads.\n", mapped*100.0/seqs->pos);

   return 0;
}


int
bw_search
(
 sublist_t  * subseqs,
 int          kmer_size,
 int          seed_loci,
 index_t    * index,
 hmargs_t     hmargs
)
{
   static const char translate[256] = {[0 ... 255] = 3,
                          ['a'] = 0, ['c'] = 1, ['g'] = 2, ['n'] = 3, ['t'] = 4,
                          ['A'] = 0, ['C'] = 1, ['G'] = 2, ['N'] = 3, ['T'] = 4 };

   int start = 0;
   clock_t tstart = clock();
   double  progress = 0.0;
   bwpos_t * cache = malloc((kmer_size+1) * sizeof(bwpos_t));
   // Insert root.
   cache[0] = (bwpos_t){0,index->gsize-1};

   for (int i = 0; i < subseqs->pos; i++) {
      // Verbose progress.
      if (hmargs.verbose && i*100.0/subseqs->pos - progress > 0.1) {
         progress = i*100.0/subseqs->pos;
         fprintf(stderr, "mapping  ...  %.1f%%\r", progress);
      }
      // Subseq.
      sub_t query = subseqs->sub[i];
      int   trail = 0;

      // Compute trail depth.
      if (i < subseqs->pos - 1) {
         sub_t next = subseqs->sub[i+1];
         while (query.seq[trail] == next.seq[trail] && trail < kmer_size) trail++;
      }
   
      // Query index.
      uint64_t    sp = cache[start].sp;
      uint64_t    ep = cache[start].ep;
      uint64_t  * c  = index->c;
      uint64_t ** occs = index->occ;
      int nt;
      int n = start;
      /* VARIANT 1: SEED THRESHOLDING
      while (n < kmer_size && ep >= seed_loci + sp) {
         nt = translate[(uint8_t)query.seq[n++]];
         sp = c[nt] + compute_occ(sp-1, occs[nt]);
         ep = c[nt] + compute_occ(ep  , occs[nt]) - 1;
         cache[n] = (bwpos_t){sp,ep};
      }
       */
      /* VARIANT 2: SEED RESSURECTION */
      uint64_t tsp = 0, tep = seed_loci;
      while (n < kmer_size && ep >= sp) {
         nt = translate[(uint8_t)query.seq[n++]];
         tsp = sp;
         tep = ep;
         sp = c[nt] + compute_occ(sp-1, occs[nt]);
         ep = c[nt] + compute_occ(ep  , occs[nt]) - 1;
         cache[n] = (bwpos_t){sp,ep};
      }
      if (ep < sp && tep - tsp < seed_loci) {
         sp = tsp;
         ep = tep;
         n--;
      }
      /* */

      if (n < trail) trail = n;

      if (ep >= sp && ep - sp < hmargs.seed_abs_max_loci) {
         vstack_t * hmap = *(query.hitmap);
         long n_hits = ep - sp + 1;
         if (hmap->pos + n_hits >= hmap->size) {
            long newsize = hmap->size + n_hits + 1;
            hmap = *(query.hitmap) = realloc(hmap, sizeof(vstack_t) + newsize * sizeof(long));
            if (hmap == NULL) {
               fprintf(stderr, "error in 'map_hits' (realloc): %s\n", strerror(errno));
               return -1;
            }
            hmap->size = newsize;
         }

         // Copy hits.
         long idstamp = 2*(query.seqid & KMERID_MASK);
         for (long k = sp; k <= ep ; k++) {
            // Delay offset to match the beginning of the sequence. Offset contains the alignment length.
            hmap->val[hmap->pos++] = ((long)(index->pos[k] + n - 1) << KMERID_BITS) | (idstamp + (n_hits > hmargs.seed_max_loci));
         }
      }

      // Update start.
      start = trail;
   }

   if (hmargs.verbose) fprintf(stderr, "mapping  ...  [%.3fms]\n", (clock()-tstart)*1000.0/CLOCKS_PER_SEC);
   return 0;
}


int
poucet_search
(
 sublist_t * subseqs,
 pstack_t ** pebbles,
 pstack_t ** hits,
 index_t   * index,
 int         tau,
 int         kmer_size,
 int         max_trail,
 int         max_loci_per_hit,
 int         abs_max_loci_per_hit,
 int         verbose
 )
{
   // Translate alphabet.
   char translate[256] = {[0 ... 255] = 3,
                          ['a'] = 0, ['c'] = 1, ['g'] = 2, ['n'] = 3, ['t'] = 4,
                          ['A'] = 0, ['C'] = 1, ['G'] = 2, ['N'] = 3, ['T'] = 4 };

   // Poucet variables.
   int start = 0;
   // Poucet algorithm over the k-mers.
   clock_t tstart = clock();
   double  progress = 0.0;
   for (int i = 0; i < subseqs->pos; i++) {
      if (verbose && i*100.0/subseqs->pos - progress > 0.1) {
         progress = i*100.0/subseqs->pos;
         fprintf(stderr, "mapping  ...  %.1f%%\r", progress);
      }
      // Extract subseq.
      sub_t query = subseqs->sub[i];
      int    qlen = kmer_size;
      int   trail = 0;
      if (i < subseqs->pos - 1) {
         // Compute trail depth.
         int k = 1;
         while (i+k < subseqs->pos && strcmp(query.seq, subseqs->sub[i+k].seq) == 0) k++;
         if (i+k < subseqs->pos) {
            sub_t next = subseqs->sub[i+k];
            while (query.seq[trail] == next.seq[trail]) trail++;
         }
      }
      trail = min(max_trail, trail);

      // Reset hits.
      for (int j = 0; j <= tau; j++) hits[j]->pos = 0;

      // Reset the pebbles that will be overwritten.
      for (int j = start+1 ; j <= trail ; j++) pebbles[j]->pos = 0;

      // Translate the query string. The first 'char' is kept to store
      // the length of the query, which shifts the array by 1 position.
      char translated[qlen+2];
      translated[qlen+1] = EOS;
      for (int j = 0 ; j < qlen ; j++) {
         translated[j+1] = translate[(int) query.seq[j]];
      }

      // Set the search options.
      struct arg_t arg = {
         .query    = translated,
         .tau      = tau,
         .trail    = trail,
         .qlen     = qlen,
         .index    = index,
         .pebbles  = pebbles,
         .hits     = hits
      };

      // Run recursive search from cached pebbles.
      char path[qlen+tau+1];

      for (int p = 0 ; p < pebbles[start]->pos ; p++) {
         // Next pebble.
         pebble_t pebble = pebbles[start]->pebble[p];
         // Compute current alignment from alignment trie.
         /*
         int wingsz;
         if (tau > 0) {
             trie_getrow(*trie, pebble.rowid >> PATH_SCORE_BITS, pebble.rowid & PATH_SCORE_MASK, &wingsz, nwrow);
         } else {
            wingsz = 0;
            *nwrow = 0;
         }
         */
         char * nwrow = pebble.nwrow + MAXTAU + 1;
         // Recover the current path.
         long gpos = index->pos[pebble.sp];
         for (int j = 0; j < start; j++) {
            path[start-j] = translate[(int)index->genome[gpos+j]];
         }
         poucet(pebble.sp, pebble.ep, nwrow, start + 1, path, &arg);
      }
      // Map hits.
      map_hits(hits, query.hitmap, index, tau, query.seqid, max_loci_per_hit, abs_max_loci_per_hit);

      start = trail;
   }
   if (verbose) fprintf(stderr, "mapping  ...  [%.3fms]\n", (clock()-tstart)*1000.0/CLOCKS_PER_SEC);
   return 0;
}



int
hitmap_analysis
(
 vstack_t    * hitmap,
 matchlist_t * matchlist,
 int           kmer_size,
 int           readlen,
 hmargs_t      hmargs
 )
{
   if (matchlist->size < 1) return -1;

   // Merge-sort loci.
   mergesort_mt(hitmap->val, hitmap->pos, sizeof(long), 0, 1, compar_long);

   long minv = 0;
   int  min  = 0;
   int  accept_d = hmargs.dist_accept;
   double rg_ratio = hmargs.read_ref_ratio;

   // Reset matchlist.
   matchlist->pos = 0;

   if (hitmap->pos == 0) return 0;

   // Find clusters in hitmap.
   long ref = 0;
   long * loc_list = malloc(hitmap->pos * sizeof(long));

   while (ref >= 0) {
      // Reference position.
      long refloc = hitmap->val[ref] >> KMERID_BITS;
      long refkmer = (hitmap->val[ref] & KMERID_MASK) >> 2;
      long refdir = (hitmap->val[ref] >> 1) & 1;
      long refend = refloc + refkmer * rg_ratio;
      int  non_significant = hitmap->val[ref] & 1;

      // Insert kmer in list.
      loc_list[0] = refloc;

      // Mark subseq as used.
      hitmap->val[ref] *= -1;

      // Start streak.
      int streak = 1;
      int sgnf_streak = non_significant;

      // Explore hitmap.
      long i = ref;
      //      long refdbg = ref;
      long streakloc = refloc;
      long streakkmer = refkmer;
      ref = -1;
      while (i < hitmap->pos - 1) {
         // Increase pointer.
         i++;
         
         // Check if this subseq has been used before.
         if (hitmap->val[i] < 0) continue;

         // Get next subsequence.
         long loc = hitmap->val[i] >> KMERID_BITS;
         long kmer = (hitmap->val[i] & KMERID_MASK) >> 2;
         long dir = (hitmap->val[i] >> 1) & 1;
         int  sgnf = hitmap->val[i] & 1;

         // Distances.
         long r_dist = refkmer - kmer;
         long g_dist = loc - refloc; // Recall that the genome is stored backwards.

         // DEBUG.
         //         if (dir == refdir) fprintf(stdout, "[%ld<->%ld]:\tend=%ld(refkmer=%ld, readlen=%d)\tr_dist=%ld\tg_dist=%ld\tfactor=%f\tkmers=[%ld:%c,%ld:%c]\tloc=[%ld,%ld]\n", refdbg, i, refend, refkmer, readlen, r_dist, g_dist, ((float)g_dist)/r_dist, refkmer, (refdir ? '+' : '-'), kmer, (dir ? '+' : '-'), refloc, loc);

         // Check if the compared sequence is too far away.
         if (loc > refend) {
            if (ref == -1) ref = i;
            break;
         }
      
         // Check if the sequences are placed similarly both in the read and the genome.
         if (dir == refdir && r_dist > 0 && ( g_dist < (rg_ratio * r_dist) || (g_dist < accept_d && r_dist < accept_d))) {
            // Save the last index of the streak.
            streakkmer = kmer;
            streakloc = loc;
            loc_list[streak++] = loc;
            if (sgnf) sgnf_streak++;
            // Mark the subsequence as used.
            hitmap->val[i] *= -1;
         }
         else {
            if (ref < 0) ref = i;
         }
      }

      int hits = kmer_size;
      for (int i = 0; i < streak - 1; i++) hits += hm_min(loc_list[i+1] - loc_list[i], kmer_size);

      // Save streak.
      // Non-significant streaks will not be saved.
      if (hits > minv && (sgnf_streak < streak)) {
         // Allocate match.
         match_t * match = malloc(sizeof(match_t));
         match->ref_e  = refloc - kmer_size;
         match->read_e = refkmer + kmer_size;
         match->dir    = refdir;
         match->ref_s  = streakloc; //- kmer_size;
         match->read_s = streakkmer;
         match->hits   = hits;
         match->flags  = 0;
         match->score  = -1;
         //         match->repeats = matchlist_new(REPEATS_SIZE);

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
         }
         else {
            match_t * match_min = matchlist->match[min];
            free(match_min);
            matchlist->match[min] = match;
         }
               
         // Find the minimum that will be substituted next.
         if (matchlist->pos == matchlist->size) {
            minv = matchlist->match[0]->hits;
            for (int j = 1 ; j < matchlist->pos; j++) {
               if (minv > matchlist->match[j]->hits) {
                  min  = j;
                  minv = matchlist->match[j]->hits;
               }
            }
         }
      }
   }

   free(loc_list);

   return 0;
}


int
map_hits
(
 pstack_t ** hits,
 vstack_t ** hitmap,
 index_t   * index,
 int         tau,
 long        id,
 int         max_loci,
 int         abs_max_loci
 )
{
   vstack_t * hmap    = *hitmap;
   long       idstamp = 2*(id & KMERID_MASK);

   // Push hits to hitmap.
   for (int a = 0 ; a <= tau ; a++) {
      for (int h = 0; h < hits[a]->pos; h++) {
         pebble_t hit    = hits[a]->pebble[h];
         long     n_hits = hit.ep - hit.sp + 1;
         // The offset between the genome start and end points in the alignment has
         // been stored in rowid during the poucet search.
         long     offset = hit.nwrow[0];
         // Filter out too abundant sequences. (To speed up the sorting).
         if (n_hits <= abs_max_loci) {
            // Realloc hitmap if needed.
            if (hmap->pos + n_hits >= hmap->size) {
               long newsize = hmap->size + n_hits + 1;
               hmap = *hitmap = realloc(hmap, sizeof(vstack_t) + newsize * sizeof(long));
               if (hmap == NULL) {
                  fprintf(stderr, "error in 'map_hits' (realloc): %s\n", strerror(errno));
                  return -1;
               }
               hmap->size = newsize;
            }

            // Copy hits.
            for (long i = hit.sp; i <= hit.ep ; i++) {
               // Delay offset to match the beginning of the sequence. Offset contains the alignment length.
               hmap->val[hmap->pos++] = ((long)(index->pos[i] + offset) << KMERID_BITS) | (idstamp + (n_hits > max_loci));
            }
         }
      }
   }
   return 0;
}

int
align_seeds
(
 seq_t          seq,
 matchlist_t  * seeds,
 matchlist_t ** seqmatches,
 index_t      * index,
 hmargs_t       hmargs
 )
{
   if (seeds->pos == 0) return 0;
   int significant = 0;
   int slen = strlen(seq.seq);
   matchlist_t * matches = *seqmatches;

   // Sort by seeded nucleotides.
   mergesort_mt(seeds->match, seeds->pos, sizeof(match_t *), 0, 1, compar_seedhits);

   for(long k = 0; k < seeds->pos; k++) {
      // Get next seed.
      match_t * seed = seeds->match[k];

      // Get query sequence.
      char * read = (seed->dir ? seq.rseq : seq.seq);

      // Compute alignment limits.
      long r_min  = seed->read_e - seed->read_s + 1;
      long l_min  = 0;
      long r_qstart = seed->read_s + 1;
      long l_qstart = hm_max(r_qstart - 1,0);
      long r_rstart = seed->ref_s - 1;
      long l_rstart = r_rstart + 1;
     
      // Extend previous alignment.
      double last_eexp = INFINITY;
      match_t * extend = NULL;
      int extend_score = 0;
      int cancel_align = 0;
      for (int i = 0; i < matches->pos; i++) {
         match_t * match = matches->match[i];
         if (match->dir != seed->dir || match->e_exp > last_eexp) continue;
         long read_s = (seed->dir ? slen - 1 - match->read_e : match->read_s);
         long read_e = (seed->dir ? slen - 1 - match->read_s : match->read_e);
         long ref_s = match->ref_s;
         long ref_e = match->ref_e;
         // Check whether the seed is contiguous.
         long r_distr = seed->read_e - read_e;
         long g_distr = ref_e - seed->ref_e;
         long r_distl = read_s - seed->read_s;
         long g_distl = seed->ref_s - ref_s;
         long d_accept = hmargs.dist_accept;
         double rg_ratio = hmargs.read_ref_ratio;
         int ext_r = g_distr > -d_accept && r_distr > -d_accept && ((g_distr < (r_distr * rg_ratio)) || (g_distr < d_accept && r_distr < d_accept));
         int ext_l = g_distl > -d_accept && r_distl > -d_accept && ((g_distl < (r_distl * rg_ratio)) || (g_distl < d_accept && r_distl < d_accept));
         if (ext_r || ext_l) {
            r_min  = (ext_r ? r_distr : 0);
            l_min  = (ext_l ? r_distl : 0);
            r_qstart = read_e + 1;
            r_rstart = ref_e - 1;
            l_qstart = read_s - 1;
            l_rstart = ref_s + 1;
            last_eexp = match->e_exp;
            extend = match;
            extend_score = match->score;
         } else if (match->e_exp < hmargs.align_accept_eexp) {
            // Check overlap
            int span = hm_min(seed->read_e - seed->read_s, match->read_e - match->read_s);
            int overlap = hm_max(0, hm_min(seed->read_e, match->read_e) - hm_max(seed->read_s, match->read_s));
            overlap = overlap > span*hmargs.overlap_max_tolerance;
            int seed_ratio = seed->hits < match->hits*hmargs.align_seed_filter_thr;
            // If overlap is higher than the maximum overlap tolerance, cancel the alignment.
            if (overlap && seed_ratio) {
               cancel_align = 1;
               break;
            }
         } else {
            if (r_distr < 0) {r_distr = -r_distr; g_distr = -g_distr;}
            if (r_distl < 0) {r_distl = -r_distl; g_distl = -g_distl;}
            int same_align = (g_distr > -d_accept && (g_distr < (r_distr*rg_ratio) || g_distr < d_accept)) || (g_distl > -d_accept && (g_distl < (r_distl*rg_ratio) || g_distl < d_accept));
            if (same_align) {
               cancel_align = 1;
               break;
            }
         }
      }
      long r_qlen = slen - r_qstart;
      long l_qlen = l_qstart + 1;
      if (cancel_align || (!r_qlen && !l_qlen)) continue;
      long r_rlen = align_min((long)(r_qlen * (1 + hmargs.align.width_ratio)), index->gsize - r_rstart);
      long l_rlen = align_min((long)(l_qlen * (1 + hmargs.align.width_ratio)), r_qstart);
      char * r_qry = read + r_qstart;
      char * l_qry = read + l_qstart;
      char * r_ref = index->genome + r_rstart;
      char * l_ref = index->genome + l_rstart;
      
      // Align forward and backward starting from (read_s, ref_s).
      long read_s, read_e, ref_s, ref_e;
      path_t align_r = (path_t){0,0,0}, align_l = (path_t){0,0,0};
      if(r_qlen) {
         align_r = dbf_align(r_qlen, r_qry, r_rlen, r_ref, r_min, ALIGN_FORWARD, ALIGN_BACKWARD, hmargs.align);
         read_e = r_qstart + align_r.row;
         ref_e = r_rstart - align_r.col;
      }
      else {
         read_e = r_qstart - 1;
         ref_e  = r_rstart + 1;
      }
      if(l_qlen) {
         align_l = dbf_align(l_qlen, l_qry, l_rlen, l_ref, l_min, ALIGN_BACKWARD, ALIGN_FORWARD, hmargs.align);
         read_s =  l_qstart - align_l.row;
         ref_s = l_rstart + align_l.col;
      }
      else {
         read_s = l_qstart + 1;
         ref_s  = l_rstart - 1;
      }
      
      // Compute significance.
      long score = extend_score + align_l.score + align_r.score;
      double ident = 1.0 - (score*1.0)/(hm_max(read_e-read_s+1,ref_e-ref_s+1));
      double e_exp = INFINITY;
      if (ident > hmargs.align_filter_ident) 
         e_exp = e_value(ref_s - ref_e + 1, extend_score + align_l.score + align_r.score, index->gsize);

      // If significant, store.
      if (e_exp < hmargs.align_filter_eexp) {
         match_t * hit;
         if (extend != NULL) {
            hit = extend;
            hit->hits += seed->hits;
         }
         else {
            hit = malloc(sizeof(match_t));
            hit->repeats = matchlist_new(REPEATS_SIZE);
            hit->flags  = 0;
            hit->hits   = seed->hits;
         }

         // Fill/Update hit.
         hit->score  = score;
         hit->ident  = ident;
         hit->ref_e  = ref_e;
         hit->ref_s  = ref_s;
         hit->read_e = (seed->dir ? slen - 1 - read_s : read_e);
         hit->read_s = (seed->dir ? slen - 1 - read_e : read_s);
         hit->dir    = seed->dir;
         hit->e_exp  = e_exp;

         // Add to significant matchlist.
         significant = 1;
         if (extend == NULL) {
            matchlist_add(seqmatches, hit);
            matches = *seqmatches;
         }
      }
            
      free(seed);
   }

   return significant;
}

int
feedback_gaps
(
 int            seqnum,
 seq_t          seq,
 matchlist_t  * intervals,
 sublist_t   ** subseqsp,
 vstack_t    ** hitmap,
 hmargs_t       hmargs,
 int            next_kmer_size,
 int            next_kmer_offset,
 char           next_qthr
)
// Find sequence gaps and feedback them to a sublist_t.
// Both unmapped gaps and intervals with identity below
// FEEDBACK_ID_THR will be mapped again after increasing tau.
{
   int slen = strlen(seq.seq);
   int recompute = 0;
   int gap_start = 0;
   sublist_t * subseqs = *subseqsp;

   for (long k = 0; k <= intervals->pos; k++) {
      int gap_end;
      int next_start = 0;
      if (k == intervals->pos) {
         gap_end = slen;
      } else {
         match_t * match = intervals->match[k];
         // If maximum identity within repeats is not enough, feedback the whole interval.
         if (match->repeats->pos > 1) {
            mergesort_mt(match->repeats->match, match->repeats->pos, sizeof(match_t *), 0, 1, compar_matcheexp);
            match = match->repeats->match[0];
         }
         if (match->e_exp > hmargs.feedback_eexp_thr) continue;
         // Check the gap length otherwise.
         gap_end = match->read_s;
         next_start = match->read_e;
      }

      if (gap_end - gap_start >= hmargs.feedback_gap_minlen) {
         int last_rec = gap_start;
         for (int j = gap_start, cnt = 0; j < gap_end; j++) {
            if (seq.q[j] < next_qthr) {
               cnt = 0;
               continue;
            } else {
               cnt++;
               if (cnt < next_kmer_size) continue;
               int idx = j + 1 - next_kmer_size;

               // Add nucleotides to recompute count.
               recompute += hm_min(j - last_rec + 1, next_kmer_size);
               last_rec = j+1;
               // Realloc sublist.
               if (subseqs->pos + 1 >= subseqs->size) {
                  int newsize = subseqs->size * 2;
                  subseqs = *subseqsp = realloc(subseqs, sizeof(sublist_t) + newsize * sizeof(sub_t));
                  if (subseqs == NULL) return -1;
                  subseqs->size = newsize;
               }
               sub_t sseq = (sub_t) {
                  .seqid = (seqnum << KMERID_BITS | (idx*2 & KMERID_MASK)),
                  .seq   = seq.seq + idx,
                  .hitmap = hitmap 
               };
               subseqs->sub[subseqs->pos++] = sseq;
               if (seq.rseq != NULL) {
                  sseq.seqid = (seqnum << KMERID_BITS | (((slen- idx - next_kmer_size)*2+1) & KMERID_MASK));
                  sseq.seq = seq.rseq + slen - idx - next_kmer_size;
                  subseqs->sub[subseqs->pos++] = sseq;
               }
               cnt = hm_max(next_kmer_size - next_kmer_offset, 0);
            }
         }
      }
      gap_start = next_start;
   }
   return recompute;
}

void
print_intervals
(
 char        * tagname,
 matchlist_t * intervals,
 index_t     * index,
 int           max_repeats
)
{
   //  int cnt = 0;
   chr_t * chr = index->chr;
   for (long k = 0; k < intervals->pos; k++) {
      match_t * match = intervals->match[k];
      if (match->repeats->pos > 1) {
         mergesort_mt(match->repeats->match, match->repeats->pos, sizeof(match_t *), 0, 1, compar_matcheexp);
         //         fprintf(stdout, "[interval %d]\t*%d* significant loci\n", ++cnt, match->repeats->pos);
         for (int j = 0; j < hm_min(max_repeats, match->repeats->pos); j++) {
            match_t * rmatch = match->repeats->match[j];
            int          dir = rmatch->dir;
            long     g_start = index->gsize - rmatch->ref_s - 1;
            long       g_end = index->gsize - rmatch->ref_e - 1;
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
            //            fprintf(stdout, "\t\t(%d,%d)\t%s:%ld-%ld:%c\t%.2f\t(%.0f%%)\t%c%c\n",
            fprintf(stdout, "%s \t%d\t%d\t%s:%ld-%ld:%c\t%.2f\t%.1f%%\t%c%c\t*\n",
                    tagname,
                    rmatch->read_s+1, rmatch->read_e+1,
                    chr->name[chrnum],
                    g_start - chr->start[chrnum]+1,
                    g_end - chr->start[chrnum]+1,
                    dir ? '-' : '+',
                    rmatch->e_exp,
                    rmatch->ident*100.0,
                    (match->flags & WARNING_OVERLAP ? 'o' : '-'),
                    (match->flags & FLAG_FUSED ? 'f' : '-'));

         }

         //         if (match->repeats->pos > max_repeats) fprintf(stdout, "\t\t...\n");
      }
      else {
         int      dir = match->dir;
         long     g_start = index->gsize - match->ref_s - 1;
         long       g_end = index->gsize - match->ref_e - 1;
         int   chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
         // Print results.
         //         fprintf(stdout, "[interval %d]\t(%d,%d)\t%s:%ld-%ld:%c\t%.2f\t(%.2f%%)\t%c%c\n",
         fprintf(stdout, "%s \t%d\t%d\t%s:%ld-%ld:%c\t%.2f\t%.1f%%\t%c%c\n",
                 //                 ++cnt,
                 tagname,
                 match->read_s+1, match->read_e+1,
                 chr->name[chrnum],
                 g_start - chr->start[chrnum]+1,
                 g_end - chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 match->e_exp,
                 match->ident*100.0,
                 (match->flags & WARNING_OVERLAP ? 'o' : '-'),
                 (match->flags & FLAG_FUSED ? 'f' : '-'));
      }
   }
}


sublist_t *
process_subseq
(
 seq_t     * seqs,
 int         numseqs,
 int         k,
 int         offset,
 char        qthr,
 vstack_t ** hitmaps
 )
/*
** Note:
** - numseqs must be smaller than 1 << KMERID_BITS.
*/
{
   if (numseqs < 1) return NULL;
   int    rev   = (seqs[0].rseq != NULL);
   char * seq   = seqs[0].seq;
   int    lsize = (strlen(seq)-k+1)*(1+rev)*numseqs;
   sublist_t * list = malloc(sizeof(sublist_t) + lsize * sizeof(sub_t));
   list->size = lsize;
   list->pos  = 0;
   if (list == NULL) return NULL;

   for (long i = 0; i < numseqs; i++) {
      seq_t s  = seqs[i];
      int slen = strlen(s.seq);
      
      // Convert to uppercase.
      for (int k = 0; k < slen; k++) if (s.seq[k] > 90) s.seq[k] -= 32;
      if (rev) for (int k = 0; k < slen; k++) if (s.rseq[k] > 90) s.rseq[k] -= 32;

      for (int j = 0, cnt = 0; j < slen; j++) {
         // Quality control.
         if (s.q[j] < qthr || s.seq[j] == 'N') {
            cnt = 0;
            continue;
         }
         else {
            cnt++;
            if (cnt < k) continue;
            int idx = j - k + 1;
            // Realloc full list.
            if (list->pos + rev >= list->size) {
               int newsize = list->size * 2;
               list = realloc(list, sizeof(sublist_t) + newsize * sizeof(sub_t));
               if (list == NULL) return NULL;
               list->size = newsize;
            }
            // Generate subseq and save to list.
            sub_t subseq = {
               .seqid = (i << KMERID_BITS) | (idx*2 & KMERID_MASK),
               .seq = s.seq + idx,
               .hitmap = hitmaps + i
            };
            list->sub[list->pos++] = subseq;
            // Insert the reverse complement as well.
            if (rev) {
               subseq.seqid = (i << KMERID_BITS) | (((slen - idx - k)*2+1) & KMERID_MASK);
               subseq.seq = s.rseq + slen - idx - k;
               list->sub[list->pos++] = subseq;
            }
            j = hm_max(j, idx + offset - 1);
            cnt = hm_max(k - offset,0);
         }
      }
   }
   return list;
}

void
fuse_matches
(
 matchlist_t ** listp,
 int            slen,
 hmargs_t       hmargs
)
{
   matchlist_t * list = *listp;

   // Sort candidate matches by genome start position.
   mergesort_mt(list->match, list->pos, sizeof(match_t *), 0, 1, compar_refstart);

   int current_pos = list->pos;
   matchlist_t * to_combine = matchlist_new(current_pos);

   for (int i = 0; i < current_pos - 1; i++) {
      match_t * match = list->match[i];
      if (match == NULL) continue;
      long max_distance = (slen - match->read_s)*hmargs.read_ref_ratio;
      to_combine->pos = 0;

      for (int j = i + 1; j < current_pos; j++) {
         match_t * cmpar = list->match[j];
         if (cmpar == NULL) continue;
         if (cmpar->dir != match->dir) continue;
         // Compute distances.
         long ref_d  = cmpar->ref_s - match->ref_s;
         long read_d = cmpar->read_s - match->read_s;
         if (match->dir) read_d = -read_d;
         // Break if we're already outside of the read scope.
         if (ref_d > max_distance) break;
         // Fusion coverage. At least fuse_min_spanratio of the fusion candidates
         // must have been aligned. (This is to avoid short random sequences to
         // be fused together spanning big regions of the read).
         long align_span = cmpar->ref_e - cmpar->ref_s + match->ref_e - match->ref_s;
         long read_span  = cmpar->ref_e - match->ref_s;
         double span_ratio = align_span * 1.0 / read_span;
         if (span_ratio < hmargs.fuse_min_spanratio) break;
         // If distances are comparable, fuse matches.
         if (ref_d <= read_d * hmargs.read_ref_ratio || (read_d < hmargs.dist_accept && ref_d < hmargs.dist_accept)) {
            matchlist_add(&to_combine, cmpar);
            list->match[j] = NULL;
         }
      }
      
      // Fuse matches.
      if (to_combine->pos > 0) {
         int totalnt = 0;
         double cmid = 0.0;
         // Matches are sorted by starting point.
         // Compute a 'crappy' mean identity. (There may be gaps that were never aligned)
         match_t * newmatch = malloc(sizeof(match_t));
         newmatch->repeats = matchlist_new(REPEATS_SIZE);

         // Copy from reference.
         newmatch->read_s = match->read_s;
         newmatch->read_e = match->read_e;
         newmatch->ref_s  = match->ref_s;
         newmatch->ref_e  = match->ref_e;
         newmatch->dir    = match->dir;
         newmatch->flags  = FLAG_FUSED;
         totalnt = match->read_e - match->read_s;
         cmid += totalnt * match->ident;

         // Expand combination.
         for (int j = 0; j < to_combine->pos; j++) {
            match_t * comb = to_combine->match[j];
            if (comb->ref_s < newmatch->ref_s) newmatch->ref_s = comb->ref_s;
            if (comb->ref_e > newmatch->ref_e) newmatch->ref_e = comb->ref_e;
            if (comb->read_e > newmatch->read_e) newmatch->read_e = comb->read_e;
            if (comb->read_s < newmatch->read_s) newmatch->read_s = comb->read_s;
            int readspan = comb->read_e - comb->read_s;
            totalnt += readspan;
            cmid += readspan * comb->ident;
         }
         newmatch->ident = cmid / totalnt;
         newmatch->score = (newmatch->read_e - newmatch->read_s) * (1 - newmatch->ident);
         matchlist_add(listp, newmatch);
         list = *listp;
         list->match[i] = NULL;
      }
   }

   // Rearrange NULL gaps.
   int gap = -1;
   int elm = 0;
   for (int i = 0; i < list->pos; i++) {
      if (list->match[i] == NULL) {
         if (gap < 0) gap = i;
      } else {
         elm++;
         if (gap >= 0) {
            list->match[gap] = list->match[i];
            list->match[i] = NULL;
            gap++;
         }
      }
   }
   
   list->pos = elm;
}

int
find_repeats
(
 matchlist_t * list,
 double        overlap
)
{
   // Sort candidate matches by genome start position.
   mergesort_mt(list->match, list->pos, sizeof(match_t *), 0, 1, compar_readstart);

   // Delete previous repeat record.
   for (int i = 0; i < list->pos; i++) list->match[i]->repeats->pos = 0;

   // Iterate over all the matches and compare overlaps.
   for (int i = 0; i < list->pos; i++) {
      match_t * ref = list->match[i];
      // Add a feedback link. This will be helpful when sorting.
      matchlist_add(&(ref->repeats), ref);
      // Compute the position of the last nucleotide.
      int last_nt = ref->read_s + (int)((ref->read_e - ref->read_s)*(1 - overlap));
      // Iterate over the other matches.
      for (int j = i+1; j < list->pos; j++) {
         match_t * cmp = list->match[j];
         if (cmp->read_s > last_nt) break;
         // Compute combined sequence span and overlap.
         int span = hm_max(cmp->read_e, ref->read_e) - hm_min(cmp->read_s, ref->read_s);
         int ovlp = hm_min(cmp->read_e, ref->read_e) - hm_max(cmp->read_s, ref->read_s);
         // Check whether the two compared matches share at least REPEAT_OVERLAP.
         if (ovlp*1.0/span >= overlap) {
            if(matchlist_add(&(ref->repeats), cmp)) return -1;
            if(matchlist_add(&(cmp->repeats), ref)) return -1;
         }
      }
   }
   return 0;
}

matchlist_t *
combine_matches
(
 matchlist_t * list,
 double        overlap_tolerance
)
{
   // Sort candidate matches by size.
   mergesort_mt(list->match, list->pos, sizeof(match_t *), 0, 1, compar_matcheexp);

   // Allocate intervals.
   matchlist_t * interv = matchlist_new(list->pos);

   // Fill read with maximum significance.
   for (int i = 0; i < list->pos; i++) {
      match_t * match = list->match[i];
      int append = 1;
      for (int j = 0; j < interv->pos; j++) {
         match_t * cmp = interv->match[j];
         int      span = hm_min(match->read_e - match->read_s, cmp->read_e - cmp->read_s);
         if (match->read_s <= cmp->read_s && cmp->read_e <= match->read_e) {
            append = 0;
            break;
         }
         else if (match->read_s < cmp->read_s) {
            int overlap = hm_max(0, match->read_e - cmp->read_s);
            if (overlap*1.0/span > overlap_tolerance) {
               append = 0;
               break;
            }
         } else if (cmp->read_e < match->read_e) {
            int overlap = hm_max(0, cmp->read_e - match->read_s);
            if (overlap*1.0/span > overlap_tolerance) {
               append = 0;
               break;
            }
         } else {
            append = 0;
            break;
         }
      }
      if (append) matchlist_add(&interv, list->match[i]);
   }

   // Sort intervals by start position.
   mergesort_mt(interv->match, interv->pos, sizeof(match_t *), 0, 1, compar_readstart);

   return interv;
}

int
fill_gaps
(
 matchlist_t ** intervp,
 matchlist_t *  matches,
 int            seq_len,
 double         gap_coverage,
 double         max_overlap
 )
// Intervals must be sorted by read start position.
// Matches that fill any gap with minimum coverage of (1-overlap_tolerance) will be included in the output.
// These matches may overlap with the contiguous 
{
   int modified = 0;
   matchlist_t * intervals = *intervp;
   mergesort_mt(matches->match, matches->pos, sizeof(match_t *), 0, 1, compar_readend);
   int gap_start = 0;
   int left_size = 0;
   int nintervals = intervals->pos;
   for (long i = 0, j = 0; i <= nintervals && j < matches->pos; i++) {
      int gap_end, right_size, next_start;
      if (i == intervals->pos) {
         gap_end = seq_len;
         right_size = 1;
         next_start = 0;
      } else {
         match_t * match = intervals->match[i];
         gap_end = match->read_s;
         next_start = match->read_e;
         right_size = next_start - gap_end;
      }

      // Absolute maximum start and end positions.
      int max_end   = gap_end + (int)(right_size * max_overlap);
      match_t * candidate = NULL;
      double e_exp = INFINITY;
      while (j < matches->pos && matches->match[j]->read_e <= max_end) {
         match_t *cmpmatch = matches->match[j++];
         // Conditions:
         // 1. Overlap with gap of at least gap_coverage.
         int gap_overlap = hm_min(cmpmatch->read_e, gap_end) - hm_max(cmpmatch->read_s, gap_start);
         if (gap_overlap*1.0/(gap_end - gap_start) < gap_coverage) continue;

         // 2. No more than max_overlap between contigs.
         int match_size = cmpmatch->read_e - cmpmatch->read_s;
         int ctg_overlap = hm_max(0, gap_start - cmpmatch->read_s);
         if (ctg_overlap*1.0/hm_min(left_size, match_size) > max_overlap) continue;
         ctg_overlap = hm_max(0, cmpmatch->read_e - gap_end);
         if (ctg_overlap*1.0/hm_min(right_size, match_size) > max_overlap) continue;

         // Conditions satisfied, take the one with best identity.
         if (cmpmatch->e_exp < e_exp) {
            candidate = cmpmatch;
            e_exp = cmpmatch->e_exp;
         }
      }

      if (candidate != NULL) {
         candidate->flags |= WARNING_OVERLAP;
         matchlist_add(intervp, candidate);
         intervals = *intervp;
         modified = 1;
      }

      left_size = right_size;
      gap_start = next_start;
   }
   if (modified) {
      mergesort_mt(intervals->match, intervals->pos, sizeof(match_t *), 0, 1, compar_readstart);
   }

   return 0;
}

double
e_value
(
 int L,
 int m,
 long gsize
)
{
   double E = log10(gsize) + log10(2) * L * (m*3.0/L-2);
   for (int i = 0; i < m; i++) E += log10((L-i)/(double)(m-i));
   return E;
}

matchlist_t *
matchlist_new
(
 int elements
)
{
   if (elements < 1) elements = 1;
   matchlist_t * list = malloc(sizeof(matchlist_t) + elements*sizeof(match_t *));
   if (list == NULL) return NULL;
   list->pos  = 0;
   list->size = elements;
   return list;
}

int
matchlist_add
(
 matchlist_t ** listp,
 match_t      * match
)
{
   matchlist_t * list = *listp;

   // Check whether stack is full.
   if (list->pos >= list->size) {
      int newsize = list->size * 2;
      *listp = list = realloc(list, sizeof(matchlist_t) + newsize * sizeof(match_t *));
      if (list == NULL) return -1;
      list->size = newsize;
   }

   // Add new match to the stack.
   list->match[list->pos++] = match;

   return 0;
}

int
compar_long
(
 const void * a,
 const void * b,
 const int param
)
{
   return *(long *)a > *(long *)b ? 1 : -1;
}


int
compar_seqsort
(
 const void * a,
 const void * b,
 const int    val
 )
{
   sub_t * seq_a = (sub_t *) a;
   sub_t * seq_b = (sub_t *) b;
   
   int compar;

   if ((compar = strncmp(seq_a->seq, seq_b->seq, val)) == 0) {
      compar = (seq_a->seqid > seq_b->seqid ? 1 : -1);
   } 
   return (compar < 0 ? -1 : 1);
}

int
compar_matchid
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   if (mb->ident > ma->ident) return 1;
   else return -1;
}

int
compar_readstart
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   if (ma->read_s > mb->read_s) return 1;
   else if (ma->read_s < mb->read_s) return -1;
   else return (ma->read_e > mb->read_e ? 1 : -1);
}

int
compar_readend
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   if (ma->read_e > mb->read_e) return 1;
   else if (ma->read_e < mb->read_e) return -1;
   else return (ma->read_s > mb->read_s ? 1 : -1);
}

int
compar_seedhits
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   if (ma->hits < mb->hits) return 1;
   else return -1;
}

int
compar_matchspan
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);

   int sza = ma->read_e - ma->read_s;
   int szb = mb->read_e - mb->read_s;
   
   if (szb > sza) return 1;
   else if (szb < sza) return -1;
   else return (mb->read_s < ma->read_s ? 1 : -1);
}

int
compar_matcheexp
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);

   if (mb->e_exp < ma->e_exp) return 1;
   else return -1;
}

int
compar_refstart
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   if (ma->ref_s > mb->ref_s) return 1;
   else if (ma->ref_s < mb->ref_s) return -1;
   else return (ma->ref_e > mb->ref_e ? 1 : -1);
}

void
free_match
(
 match_t * match
 )
{
   free(match->repeats);
   free(match);
}
