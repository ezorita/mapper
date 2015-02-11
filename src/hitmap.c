#include "hitmap.h"

int
hitmap
(
 index_t    * index,
 chr_t      * chr,
 seqstack_t * seqs,
 hmargs_t     hmargs
 )
{
   int maxtau = hmargs.maxtau;

   // Poucet score trie.
   trie_t * trie = trie_new(TRIE_SIZE);

   // Initialize pebble and hit stacks.
   pstack_t ** pebbles = malloc((MAX_TRAIL+1)*sizeof(pstack_t *));
   pstack_t ** hits    = malloc((maxtau+1) * sizeof(pstack_t *));
   for (int i = 0 ; i <= MAX_TRAIL ; i++) pebbles[i] = new_pstack(PBSTACK_SIZE);
   for (int i = 0 ; i <= maxtau ; i++)    hits[i]    = new_pstack(HITSTACK_SIZE);

   // Push root node.
   pebble_t root = {.sp = 0, .ep = index->gsize-1, .rowid = 0};
   ppush(pebbles, root);

   // Define number of seqs that will be mixed in each poucet search.
   int seq_block = 10000;

   // Allocate hitmaps.
   vstack_t ** hitmaps = malloc(seq_block * sizeof(vstack_t *));
   for (int i = 0 ; i < seq_block ; i++) hitmaps[i] = new_stack(HITMAP_SIZE);

   // Allocate loci stack.
   matchlist_t * seeds = matchlist_new(hmargs.max_align_per_read);

   // Allocate match buffer.
   matchlist_t ** seqmatches;
   seqmatches = malloc(seq_block * sizeof(matchlist_t *));
   for (int i = 0; i < seq_block; i++) seqmatches[i] = matchlist_new(HITBUF_SIZE);

   // Iterate over sequence blocks.
   for (int s = 0 ; s < seqs->pos; s += seq_block) {

      // Pre-process sequence block: chunk in subsequences and sort.
      int         numseqs = (seq_block < seqs->pos - s ? seq_block : seqs->pos - s);
      sublist_t * subseqs = process_subseq(seqs->seq+s, numseqs, hmargs.kmer_size, hitmaps);

      // Verbose processed blocks.
      if(hmargs.verbose) fprintf(stderr, "processing reads from %d to %d\n", s+1, s+numseqs);

      // Iterate over increasing mismatch tolerance.
      for (int a = 0; a <= maxtau; a++) {
         // Verbose sorting.
         if(hmargs.verbose) fprintf(stderr, "[tau=%d] %d subsequences\nsorting  ...", a, subseqs->size);
         clock_t tstart = clock();
         // Sort subsequences.
         mergesort_mt(subseqs->sub, subseqs->size, sizeof(sub_t), hmargs.kmer_size, 1, compar_seqsort);
         // Verbose sort time.
         if(hmargs.verbose) fprintf(stderr, "  [%.3fms]\n", (clock()-tstart)*1000.0/CLOCKS_PER_SEC);

         // Clear hitmaps.
         for (int i = 0 ; i < numseqs ; i++) hitmaps[i]->pos = 0;

         // Search subsequences.
         poucet_search(subseqs, pebbles, hits, &trie, index, a, hmargs.kmer_size, hmargs.kmer_size, hmargs.seed_max_loci, hmargs.verbose);

         // Reset the subsequence list. It will be filled up with the unmatched regions.
         subseqs->size = 0;

         // Sort all hitmaps
         long recompute = 0, matched = 0, totalnt = 0;
         // Reset verbose variables.
         double progress = 0.0;
         tstart = clock();
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
            hitmap_analysis(hitmaps[i], seeds, hmargs.kmer_size, slen, hmargs);

            // Smith-Waterman alignment.
            int newdata = align_seeds(seqs->seq[i], seeds, seqmatches + i, index, hmargs);

             if (newdata) {
               // Fuse matches.
               fuse_matches(seqmatches+i, slen, hmargs);

               // Find sequence repeats.
               find_repeats(seqmatches[i], hmargs.repeat_min_overlap);
            }
            // Assemble read.
             matchlist_t * intervals = combine_matches(seqmatches[i], hmargs.overlap_tolerance);

            if (a < maxtau) {
               // Feedback gaps between intervals and low identity matches.
               recompute += feedback_gaps(i, seqs->seq[i], intervals, subseqs, hitmaps + i, hmargs);
            } 
            else {
               // Print results if a == maxtau.
               // Fill read gaps allowing greater overlap.
               fill_gaps(&intervals, seqmatches[i], slen, hmargs.match_min_len, 1 - hmargs.overlap_tolerance, hmargs.overlap_max_tolerance);

               // Print intervals.
               if (intervals->pos) {
                  fprintf(stdout, "%s\n", seqs->seq[s+i].tag);
                 print_intervals(intervals, chr, hmargs.repeat_print_num);
               }
            }

            // Compute matched nucleotides.
            for (int i = 0; i < intervals->pos; i++) {
               match_t * current = intervals->match[i];
               long next_s;
               if (i < intervals->pos - 1) next_s = intervals->match[i+1]->read_s;
               else next_s = slen;
               matched += hm_min(next_s, current->read_e) - current->read_s;
            }

            free(intervals);
         }
         if (hmargs.verbose) {
            fprintf(stderr, "aligning ...  [%.3fms]\n", (clock()-tstart)*1000.0/CLOCKS_PER_SEC);
            fprintf(stderr, "matched %ld out of %ld nt (%.1f%%), %ld nt will be queried again.\n", matched, totalnt, matched*100.0/totalnt, recompute);
         }
      }
   }
   

   return 0;
}


int
poucet_search
(
 sublist_t * subseqs,
 pstack_t ** pebbles,
 pstack_t ** hits,
 trie_t   ** trie,
 index_t   * index,
 int         tau,
 int         kmer_size,
 int         max_trail,
 int         max_loci_per_hit,
 int         verbose
 )
{
   // Translate alphabet.
   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                          ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                          ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

   // Poucet variables.
   int start = 0;
   // Poucet algorithm over the k-mers.
   clock_t tstart = clock();
   double  progress = 0.0;
   for (int i = 0; i < subseqs->size; i++) {
      if (verbose && i*100.0/subseqs->size - progress > 0.1) {
         progress = i*100.0/subseqs->size;
         fprintf(stderr, "mapping  ...  %.1f%%\r", progress);
      }
      // Extract subseq.
      sub_t query = subseqs->sub[i];
      int    qlen = kmer_size;
      int   trail = 0;
      if (i < subseqs->size - 1) {
         // Compute trail depth.
         int k = 1;
         while (i+k < subseqs->size && strcmp(query.seq, subseqs->sub[i+k].seq) == 0) k++;
         if (i+k < subseqs->size) {
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
         .triep    = trie,
         .pebbles  = pebbles,
         .hits     = hits
      };

      // Run recursive search from cached pebbles.
      uint row[2*tau+3];
      uint * nwrow = row + tau + 1;
      char path[qlen+tau+1];

      for (int p = 0 ; p < pebbles[start]->pos ; p++) {
         // Next pebble.
         pebble_t pebble = pebbles[start]->pebble[p];
         // Compute current alignment from alignment trie.
         int wingsz;
         if (tau > 0) trie_getrow(*trie, pebble.rowid >> PATH_SCORE_BITS, pebble.rowid & PATH_SCORE_MASK, &wingsz, nwrow);
         else {
            wingsz = 0;
            *nwrow = 0;
         }
         // Recover the current path.
         long gpos = index->pos[pebble.sp];
         for (int j = 0; j < start; j++) {
            path[start-j] = translate[(int)index->genome[gpos+j]];
         }
         poucet(pebble.sp, pebble.ep, wingsz, nwrow, start + 1, path, &arg);
      }
      // Map hits.
      map_hits(hits, query.hitmap, index, tau, query.seqid, max_loci_per_hit);

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
   // Copy to buffer and sort.
   long * sortbuf = malloc(hitmap->pos * sizeof(long));
   memcpy(sortbuf, hitmap->val, hitmap->pos * sizeof(long));
   mergesort_long(hitmap->val, sortbuf, hitmap->pos, 0);
   free(sortbuf);

   long minv = 0;
   int  min  = 0;
   int  rg_ratio = hmargs.read_ref_ratio;
   int  accept_d = hmargs.dist_accept;

   // Reset matchlist.
   matchlist->pos = 0;

   if (hitmap->pos == 0) return 0;

   // Find clusters in hitmap.
   long ref = 0;

   while (ref >= 0) {
      // Reference position.
      long refloc = hitmap->val[ref] >> KMERID_BITS;
      long refkmer = (hitmap->val[ref] & KMERID_MASK) >> 1;
      long refdir = hitmap->val[ref] & 1;
      long refend = refloc + refkmer * rg_ratio;

      // Mark subseq as used.
      hitmap->val[ref] *= -1;

      // Start streak.
      int streak = 1;

      // Explore hitmap.
      long i = ref;
      //long refdbg = ref;
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
         long kmer = (hitmap->val[i] & KMERID_MASK) >> 1;
         long dir = hitmap->val[i] & 1;

         // Distances.
         long r_dist = refkmer - kmer;
         long g_dist = loc - refloc; // Recall that the genome is stored backwards.

         // DEBUG.
         //fprintf(stdout, "[%ld<->%ld]:\tend=%ld(refkmer=%ld, readlen=%d)\tr_dist=%ld\tg_dist=%ld\tfactor=%f\tkmers=[%ld:%c,%ld:%c]\tloc=[%ld,%ld]\n", refdbg, i, refend, refkmer, readlen, r_dist, g_dist, ((float)g_dist)/r_dist, refkmer, (refdir ? '+' : '-'), kmer, (dir ? '+' : '-'), refloc, loc);

         // Check if the compared sequence is too far away.
         if (loc > refend) {
            if (ref == -1) ref = i;
            break;
         }
      
         // Check if the sequences are placed similarly both in the read and the genome.
         if (dir == refdir && ((r_dist > 0 && g_dist < (rg_ratio * r_dist)) || (g_dist < accept_d && r_dist < accept_d && -r_dist < accept_d))) {
            // Save the last index of the streak.
            streakkmer = kmer;
            streakloc = loc;
            streak++;
            // Mark the subsequence as used.
            hitmap->val[i] *= -1;
         }
         else {
            if (ref < 0) ref = i;
         }
      }

      // Save streak.
      if (streak > minv) {
         // Allocate match.
         match_t * match = malloc(sizeof(match_t));
         match->ref_e  = refloc;
         match->read_e = refkmer;
         match->dir    = refdir;
         match->ref_s  = streakloc;
         match->read_s = streakkmer;
         // TODO:
         // Nothing extra to do with reverse strand??
         match->hits   = streak;
         match->flags  = 0;
         match->score  = -1;
         match->repeats = matchlist_new(REPEATS_SIZE);

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
         }
         else {
            match_t * match_min = matchlist->match[min];
            free_match(match_min);
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
 int         max_loci
 )
{
   vstack_t * hmap    = *hitmap;
   long       idstamp = id & KMERID_MASK;

   // Push hits to hitmap.
   for (int a = 0 ; a <= tau ; a++) {
      for (int h = 0; h < hits[a]->pos; h++) {
         pebble_t hit    = hits[a]->pebble[h];
         long     n_hits = hit.ep - hit.sp + 1;
         // The offset between the genome start and end points in the alignment has
         // been stored in rowid during the poucet search.
         long     offset = hit.rowid;
         // Filter out too abundant sequences. (To speed up the sorting).
         if (n_hits <= max_loci) {
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
               hmap->val[hmap->pos++] = ((long)(index->pos[i] + offset) << KMERID_BITS) | idstamp;
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
   int significant = 0;
   int slen = strlen(seq.seq);
   for(long k = 0; k < seeds->pos; k++) {
      // Get next seed.
      match_t * seed = seeds->match[k];
      // Extend alignment.
      align_t align_r, align_l;
      char * read = (seed->dir ? seq.rseq : seq.seq);

      // Will align both from the start point (backward from start, forward from start+1).
      // The previously matched region will be aligned again,
      // but it's the only way to know precisely the total identity of the read.
      // Note that the genome is stored backwards so all signs are changed.
      align_r = nw_align(read + seed->read_s + 1, index->genome + seed->ref_s - 1, slen - seed->read_s - 1, ALIGN_FORWARD, ALIGN_BACKWARD, hmargs.align_likelihood_thr, hmargs.read_match_prob, hmargs.rand_match_prob);
      align_l = nw_align(read + seed->read_s, index->genome + seed->ref_s, seed->read_s + 1, ALIGN_BACKWARD, ALIGN_FORWARD, hmargs.align_likelihood_thr, hmargs.read_match_prob, hmargs.rand_match_prob);

      seed->score  = align_l.score + align_r.score;
      seed->ident  = 1.0 - (seed->score*1.0)/(align_l.pathlen + align_r.pathlen);
      seed->ref_e  = (index->gsize - seed->ref_s) + align_r.end + 2; // Align.end is the genome breakpoint.
      seed->ref_s  = (index->gsize - seed->ref_s) - align_l.end; // Align.end is the genome breakpoint.
      seed->read_e = seed->read_s + align_r.start + 2; // Align.start is the read breakpoint.
      seed->read_s = seed->read_s - align_l.start;     // Align.start is the read breakpoint.

      // Filter out seeds that do not meet the minimum quality.
      if ((seed->ref_e - seed->ref_s < hmargs.match_min_len) || (seed->ident < hmargs.match_min_id)) {
         free_match(seed);
         continue;
      }
            
      // Correct read start and end points if we are aligning the reverse complement.
      if (seed->dir) {
         long end = slen - seed->read_s;
         seed->read_s = slen - seed->read_e;
         seed->read_e = end;
      }

      // Add to significant matchlist.
      significant = 1;
      matchlist_add(seqmatches, seed);
   }

   return significant;
}

int
feedback_gaps
(
 int            seqnum,
 seq_t          seq,
 matchlist_t  * intervals,
 sublist_t    * subseqs,
 vstack_t    ** hitmap,
 hmargs_t       hmargs
)
// Find sequence gaps and feedback them to a sublist_t.
// Both unmapped gaps and intervals with identity below
// FEEDBACK_ID_THR will be mapped again after increasing tau.
{
   int slen = strlen(seq.seq);
   int recompute = 0;
   int gap_start = 0;

   for (long k = 0; k <= intervals->pos; k++) {
      int gap_end;
      int next_start = 0;
      if (k == intervals->pos) {
         gap_end = slen;
      } else {
         match_t * match = intervals->match[k];
         // If maximum identity within repeats is not enough, feedback the whole interval.
         if (match->repeats->pos > 1) {
            mergesort_mt(match->repeats->match, match->repeats->pos, sizeof(match_t *), 0, 1, compar_matchid);
            match = match->repeats->match[0];
         }
         if (match->ident < hmargs.feedback_id_thr) continue;
         // Check the gap length otherwise.
         gap_end = match->read_s;
         next_start = match->read_e;
      }

      if (gap_end - gap_start >= hmargs.match_min_len) {
         recompute += gap_end - gap_start;
         for (int j = gap_start; j <= gap_end - hmargs.kmer_size; j++) {
            sub_t sseq = (sub_t) {
               .seqid = (seqnum << KMERID_BITS | (j*2 & KMERID_MASK)),
               .seq   = seq.seq + j,
               .hitmap = hitmap 
            };
            subseqs->sub[subseqs->size++] = sseq;
            if (seq.rseq != NULL) {
               sseq.seqid = (seqnum << KMERID_BITS | (((slen- j - hmargs.kmer_size)*2+1) & KMERID_MASK));
               sseq.seq = seq.rseq + slen - j - hmargs.kmer_size;
               subseqs->sub[subseqs->size++] = sseq;
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
 matchlist_t * intervals,
 chr_t       * chr,
 int           max_repeats
)
{
   int cnt = 0;
   for (long k = 0; k < intervals->pos; k++) {
      match_t * match = intervals->match[k];
      if (match->repeats->pos > 1) {
         mergesort_mt(match->repeats->match, match->repeats->pos, sizeof(match_t *), 0, 1, compar_matchid);
         fprintf(stdout, "[interval %d]\t*%d* significant loci\n", ++cnt, match->repeats->pos);
         for (int j = 0; j < hm_min(max_repeats, match->repeats->pos); j++) {
            match_t * rmatch = match->repeats->match[j];
            int          dir = rmatch->dir;
            long     g_start = rmatch->ref_s;
            long       g_end = rmatch->ref_e - 1;
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
            fprintf(stdout, "\t\t(%d,%d)\t%s:%ld-%ld:%c\t(%.2f%%)\t%c%c\n",
                    rmatch->read_s, rmatch->read_e - 1,
                    chr->name[chrnum],
                    g_start - chr->start[chrnum]+1,
                    g_end - chr->start[chrnum]+1,
                    dir ? '-' : '+',
                    rmatch->ident*100.0,
                    (match->flags & WARNING_OVERLAP ? 'o' : '-'),
                    (match->flags & FLAG_FUSED ? 'f' : '-'));

         }
         if (match->repeats->pos > max_repeats) fprintf(stdout, "\t\t...\n");
      }
      else {
         int      dir = match->dir;
         long g_start = match->ref_s;
         long   g_end = match->ref_e - 1;
         int   chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "[interval %d]\t(%d,%d)\t%s:%ld-%ld:%c\t(%.2f%%)\t%c%c\n",
                 ++cnt,
                 match->read_s, match->read_e - 1,
                 chr->name[chrnum],
                 g_start - chr->start[chrnum]+1,
                 g_end - chr->start[chrnum]+1,
                 dir ? '-' : '+',
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
   int    lpos  = 0;
   sublist_t * list = malloc(sizeof(sublist_t) + lsize * sizeof(sub_t));
   if (list == NULL) return NULL;

   for (long i = 0; i < numseqs; i++) {
      seq_t s  = seqs[i];
      int slen = strlen(s.seq);
      int subs = slen - k + 1;
      
      // Convert to uppercase.
      for (int k = 0; k < slen; k++)
         if (s.seq[k] > 90) s.seq[k] -= 32;
      if (rev) for (int k = 0; k < slen; k++)
                  if (s.rseq[k] > 90) s.rseq[k] -= 32;

      for (int j = 0; j < subs; j++) {
         // Realloc full list.
         if (lpos + rev >= lsize) {
            lsize *=2;
            list = realloc(list, sizeof(sublist_t) + lsize * sizeof(sub_t));
            if (list == NULL) return NULL;
         }
         // Generate subseq and save to list.
         sub_t subseq = {
            .seqid = (i << KMERID_BITS) | (j*2 & KMERID_MASK),
            .seq = s.seq + j,
            .hitmap = hitmaps + i
         };
         list->sub[lpos++] = subseq;
         // Insert the reverse complement as well.
         if (rev) {
            subseq.seqid++;
            subseq.seq = s.rseq + j;
            list->sub[lpos++] = subseq;
         }
      }
   }
   // Realloc list.
   list->size = lpos;
   list = realloc(list, sizeof(sublist_t) + lpos * sizeof(sub_t));

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
   mergesort_mt(list->match, list->pos, sizeof(match_t *), 0, 1, compar_matchspan);

   // Allocate intervals.
   matchlist_t * interv = matchlist_new(list->pos);

   // Fill read with maximum coverage.
   for (int i = 0; i < list->pos; i++) {
      match_t * match = list->match[i];
      int append = 1;
      for (int j = 0; j < interv->pos; j++) {
         match_t * cmp = interv->match[j];
         if (match->read_s < cmp->read_s) {
            // match must be smaller in size because of the sorting.
            int span    = match->read_e - match->read_s;
            int overlap = hm_max(0, match->read_e - cmp->read_s);
            if (overlap*1.0/span > overlap_tolerance) {
               append = 0;
               break;
            }
         } else if (match->read_e > cmp->read_e) {
            // match must be smaller in size because of the sorting.
            int span    = match->read_e - match->read_s;
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
 int            match_minlen,
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
      // TODO:
      // This overlap must be the minimum 50% between:
      // left: the gap and the previous match.
      // right: the gap and the next match.
      int gap_size = gap_end - gap_start;
      if (gap_size >= match_minlen) {
         // Absolute maximum start and end positions.
         int max_end   = gap_end + (int)(right_size * max_overlap);
         match_t * candidate = NULL;
         double id = 0.0;
         while (j < matches->pos && matches->match[j]->read_e <= max_end) {
            match_t *cmpmatch = matches->match[j++];
            // Conditions:
            // 1. Overlap with gap of at least gap_coverage.
            int gap_overlap = hm_min(cmpmatch->read_e, gap_end) - hm_max(cmpmatch->read_s, gap_start);
            if (gap_overlap*1.0/gap_size < gap_coverage) continue;

            // 2. No more than max_overlap between contigs.
            int match_size = cmpmatch->read_e - cmpmatch->read_s;
            int ctg_overlap = hm_max(0, gap_start - cmpmatch->read_s);
            if (ctg_overlap*1.0/hm_min(left_size, match_size) > max_overlap) continue;
            ctg_overlap = hm_max(0, cmpmatch->read_e - gap_end);
            if (ctg_overlap*1.0/hm_min(right_size, match_size) > max_overlap) continue;

            // Conditions satisfied, take the one with best identity.
            if (cmpmatch->ident > id) {
               candidate = cmpmatch;
               id = cmpmatch->ident;
            }
         }

         if (candidate != NULL) {
            candidate->flags |= WARNING_OVERLAP;
            matchlist_add(intervp, candidate);
            intervals = *intervp;
            modified = 1;
         }
      }
      left_size = right_size;
      gap_start = next_start;
   }
   if (modified) {
      mergesort_mt(intervals->match, intervals->pos, sizeof(match_t *), 0, 1, compar_readstart);
   }

   return 0;
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
