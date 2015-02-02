#include "hitmap.h"

int
hitmap
(
 int          maxtau,
 index_t      index,
 chr_t      * chr,
 seqstack_t * seqs
 )
{
   // Translate alphabet.
   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                          ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                          ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

   // Map vars / structs.
   long   sortbuf_size = SORTBUF_SIZE;
   long * sortbuf      = malloc(sortbuf_size * sizeof(long));
   trie_t * trie = trie_new(TRIE_SIZE);

   // Initialize pebble and hit stacks.
   pstack_t ** pebbles = malloc((MAX_TRAIL+1)*sizeof(pstack_t *));
   pstack_t ** hits    = malloc((maxtau+1) * sizeof(pstack_t *));
   for (int i = 0 ; i <= MAX_TRAIL ; i++) {
      pebbles[i] = new_pstack(PBSTACK_SIZE);
   }
   for (int i = 0 ; i <= maxtau ; i++) {
      hits[i]    = new_pstack(HITSTACK_SIZE);
   }

   // Push root node.
   pebble_t root = {
      .sp = 0,
      .ep = index.gsize-1,
      .rowid = 0
   };
   ppush(pebbles, root);

   // Define number of seqs that will be mixed in each poucet search.
   int seq_block = 10000;
   // Allocate hitmaps.
   vstack_t ** hitmaps = malloc(seq_block * sizeof(vstack_t *));
   for (int i = 0 ; i < seq_block ; i++) {
      hitmaps[i] = new_stack(HITMAP_SIZE);
      if (hitmaps[i] == NULL) return -1;
   }
   // Allocate loci stack.
   matchlist_t * matches = matchlist_new(MATCHLIST_SIZE);
   if (matches == NULL) return -1;
   // Allocate match buffer.
   matchlist_t ** seqmatches;
   seqmatches = malloc(seq_block * sizeof(matchlist_t *));
   for (int i = 0; i < seq_block; i++)
      seqmatches[i] = matchlist_new(HITBUF_SIZE);

   for (int s = 0 ; s < seqs->pos; s += seq_block) {

      /**                                         **
       **         PRE-PROCESS SUBSEQUENCES        **
       **                                         **/

      // Pre-process sequence block: chunk in subsequences and sort.
      int         numseqs = (seq_block < seqs->pos - s ? seq_block : seqs->pos - s);
      sublist_t * subseqs = process_subseq(seqs->seq+s, numseqs, KMER_SIZE, hitmaps);

      fprintf(stderr, "Processing reads from %d to %d\n", s+1, s+numseqs);

      for (int a = 0; a <= maxtau; a++) {

         fprintf(stderr, "[tau=%d] %d subsequences\n", a, subseqs->size);

         // Sort subsequences.
         mergesort_mt(subseqs->sub, subseqs->size, sizeof(sub_t), KMER_SIZE, 1, compar_seqsort);
     
         /**                                         **
          **             MAP SUBSEQUENCES            **
          **                                         **/

         // Clear hitmaps.
         for (int i = 0 ; i < numseqs ; i++) {
            hitmaps[i]->pos = 0;
         }

         // Poucet variables.
         int start = 0;
         clock_t tstart = clock();

         // Poucet algorithm over the k-mers.
         for (int i = 0; i < subseqs->size; i++) {
            if (i%100 == 0) fprintf(stderr, "mapping ... %7d/%7d\r", i, subseqs->size);
            // Extract subseq.
            sub_t query = subseqs->sub[i];
            int    qlen = KMER_SIZE;
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
            trail = min(MAX_TRAIL, trail);

            // Reset hits.
            for (int j = 0; j <= a; j++) {
               hits[j]->pos = 0;
            }

            // Reset the pebbles that will be overwritten.
            for (int j = start+1 ; j <= trail ; j++) {
               pebbles[j]->pos = 0;
            }

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
               .tau      = a,
               .trail    = trail,
               .qlen     = qlen,
               .index    = &index,
               .triep    = &trie,
               .pebbles  = pebbles,
               .hits     = hits
            };

            // Run recursive search from cached pebbles.
            uint row[2*a+3];
            uint * nwrow = row + a + 1;
            char path[qlen+a+1];

            for (int p = 0 ; p < pebbles[start]->pos ; p++) {
               // Next pebble.
               pebble_t pebble = pebbles[start]->pebble[p];
               // Compute current alignment from alignment trie.
               int wingsz;
               if (a > 0) trie_getrow(trie, pebble.rowid >> PATH_SCORE_BITS, pebble.rowid & PATH_SCORE_MASK, &wingsz, nwrow);
               else {
                  wingsz = 0;
                  *nwrow = 0;
               }
               // Recover the current path.
               long gpos = index.pos[pebble.sp];
               for (int j = 0; j < start; j++) {
                  path[start-j] = translate[(int)index.genome[gpos+j]];
               }
               poucet(pebble.sp, pebble.ep, wingsz, nwrow, start + 1, path, &arg);
            }
            // Map hits.
            map_hits(hits, query.hitmap, &index, a, query.seqid);

            start = trail;
         }
         fprintf(stderr, "mapping ... %7d/%7d - hitmap built\t[%.3fms]\n", subseqs->size, subseqs->size, (clock()-tstart)*1000.0/CLOCKS_PER_SEC);

         // Reset the subsequence list. It will be filled up with the unmatched regions.
         subseqs->size = 0;

         tstart = clock();

         // Sort all hitmaps
         for (int i = 0 ; i < numseqs ; i++) {
            fprintf(stderr, "aligning... %7d/%7d\r", i, numseqs);


            /**                                         **
             **             HITMAP ANALYSIS             **
             **                                         **/

            // Sequence length.
            int slen = strlen(seqs->seq[i].seq);

            // Realloc sort buffer if necessary.
            if (sortbuf_size < hitmaps[i]->pos) {
               sortbuf = realloc(sortbuf, hitmaps[i]->pos * sizeof(long));
               if (sortbuf == NULL) return -1;
            }
            // Merge-sort loci.
            memcpy(sortbuf, hitmaps[i]->val, hitmaps[i]->pos * sizeof(long));
            mergesort_long(hitmaps[i]->val, sortbuf, hitmaps[i]->pos, 0);

            // Find matching loci in hitmaps.
            //fprintf(stderr, "starting hitmap analysis: %ld candidate loci.\n", hitmaps[i]->pos);
            hitmap_analysis(hitmaps[i], matches, KMER_SIZE, slen);
            //fprintf(stderr, "hitmap analysis: %d alignments will be performed.\n", matches->pos);
            if (matches->pos == 0) continue;

            /**                                         **
             **           SEQUENCE ALIGNMENT            **
             **                                         **/

            // Smith-Waterman alignment.
            for(long k = 0; k < matches->pos; k++) {
               // Get next match.
               match_t * match = matches->match[k];
               // Extend alignment.
               align_t align_r, align_l;
               char * read = (match->dir ? seqs->seq[i].rseq : seqs->seq[i].seq);

               // Will align both from the start point (backward from start, forward from start+1).
               // The previously matched region will be aligned again,
               // but it's the only way to know precisely the total identity of the read.
               // Note that the genome is stored backwards so all signs are changed.
               align_r = nw_align(read + match->read_s + 1, index.genome + match->ref_s - 1, slen - match->read_s - 1, ALIGN_FORWARD, ALIGN_BACKWARD);
               align_l = nw_align(read + match->read_s, index.genome + match->ref_s, match->read_s, ALIGN_BACKWARD, ALIGN_FORWARD);

               match->score  = align_l.score + align_r.score;
               match->ident  = 1.0 - (match->score*1.0)/(align_l.pathlen + align_r.pathlen);
               match->ref_e  = index.gsize - match->ref_s + 1 + align_r.end; // Align.end is the genome breakpoint.
               match->ref_s  = index.gsize - match->ref_s - align_l.end; // Align.end is the genome breakpoint.
               match->read_e = match->read_s + 1 + align_r.start; // Align.start is the read breakpoint.
               match->read_s = match->read_s - align_l.start;     // Align.start is the read breakpoint.

               // Filter out matches that do not meet the minimum quality.
               if ((match->ref_e - match->ref_s < SEQ_MINLEN) || (match->ident < SEQ_MINID)) {
                  free(match->repeats);
                  free(match);
                  continue;
               }
            
               // Correct read start and end points if we are matching the reverse complement.
               if (match->dir) {
                  long end = slen - match->read_s;
                  match->read_s = slen - match->read_e;
                  match->read_e = end;
               }

               // DEBUG.
               /*
               long g_start = match->ref_s;
               long g_end   = match->ref_e;
               int chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
               fprintf(stdout, "%d\t%d (%.2f%%)\t%d-%d\t%s:%ld-%ld:%c (%ld)\tr=%d\n",
                       match->hits,
                       match->score,
                       match->ident*100.0,
                       match->read_s, match->read_e,
                       chr->name[chrnum],
                       g_start - chr->start[chrnum]+1,
                       g_end - chr->start[chrnum]+1,
                       match->dir ? '-' : '+',
                       g_end - g_start + 1,
                       match->repeats->pos);
               */

               // Add to significant matchlist.
               matchlist_add(seqmatches+i, match);
            }


            /**                                         **
             **              ASSEMBLE READ              **
             **                                         **/
            // DEBUG.
            //            fprintf(stderr, "matches before read re-assembly %d\n", seqmatches[i]->pos);
            
            // Sort candidate matches by genome start position.
            mergesort_mt(seqmatches[i]->match, seqmatches[i]->pos, sizeof(match_t *), 0, 1, compar_refstart);

            // Fuse matches.
            fuse_matches(seqmatches+i, slen);

            // Sort candidate matches by genome start position.
            mergesort_mt(seqmatches[i]->match, seqmatches[i]->pos, sizeof(match_t *), 0, 1, compar_matchstart);

            // Find sequence repeats.
            find_repeats(seqmatches[i]);

            // Sort candidate matches by size.
            mergesort_mt(seqmatches[i]->match, seqmatches[i]->pos, sizeof(match_t *), 0, 1, compar_matchsize);

            // Assemble read.
            matchlist_t * intervals = combine_matches(seqmatches[i]);

            // Compute matched nucleotides.
            long matched = 0;
            for (int i = 0; i < intervals->pos; i++) {
               match_t * current = intervals->match[i];
               long next_s;
               if (i < intervals->pos - 1) next_s = intervals->match[i+1]->read_s;
               else next_s = slen;
               matched += hm_min(next_s, current->read_e) - current->read_s;
            }

            // Intervals:
            if (a == maxtau && seqmatches[i]->pos) fprintf(stdout, "%s\n", seqs->seq[s+i].tag);
            for (long k = 0, cnt = 0 ; k < intervals->pos; k++) {
               match_t * match = intervals->match[k];
               if (a == maxtau) {
                  if (match->repeats->pos > 1) {
                     mergesort_mt(match->repeats->match, match->repeats->pos, sizeof(match_t *), 0, 1, compar_matchid);
                     fprintf(stdout, "[interval %ld]\t*%d* significant loci\n", ++cnt, match->repeats->pos);
                     for (int j = 0; j < hm_min(PRINT_REPEATS_NUM, match->repeats->pos); j++) {
                        match_t * rmatch = match->repeats->match[j];
                        int          dir = rmatch->dir;
                        long     g_start = rmatch->ref_s;
                        long       g_end = rmatch->ref_e;
                        int       chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
                        fprintf(stdout, "\t\t(%d,%d)\t%s:%ld-%ld:%c\t(%.2f%%)\n",
                                rmatch->read_s, rmatch->read_e,
                                chr->name[chrnum],
                                g_start - chr->start[chrnum]+1,
                                g_end - chr->start[chrnum]+1,
                                dir ? '-' : '+',
                                rmatch->ident*100.0);
                     }
                     if (match->repeats->pos > PRINT_REPEATS_NUM) fprintf(stdout, "\t\t...\n");
                     match = match->repeats->match[0];
                  } else {
                     int      dir = match->dir;
                     long g_start = match->ref_s;
                     long   g_end = match->ref_e;
                     int   chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
                     // Print results.
                     fprintf(stdout, "[interval %ld]\t(%d,%d)\t%s:%ld-%ld:%c\t(%.2f%%)\n",
                             ++cnt,
                             match->read_s, match->read_e,
                             chr->name[chrnum],
                             g_start - chr->start[chrnum]+1,
                             g_end - chr->start[chrnum]+1,
                             dir ? '-' : '+',
                             match->ident*100.0);
 
                  }
               }
               // Find sequence gaps and feedback them to a sublist_t.
               // Both unmapped gaps and intervals with identity below
               // INTERVAL_MINID will be mapped again after increasing tau.
               int gapend = (k < intervals->pos - 1 ? intervals->match[k+1]->read_s : slen);
               if (gapend - match->read_e > SEQ_MINLEN || match->ident < INTERVAL_MINID) {
                  for (int j = match->read_e; j <= gapend - KMER_SIZE; j++) {
                     sub_t sseq = (sub_t) {
                        .seqid = (i << KMERID_BITS | (j*2 & KMERID_MASK)),
                        .seq   = seqs->seq[i].seq + j,
                        .hitmap = hitmaps + i 
                     };
                     subseqs->sub[subseqs->size++] = sseq;
                     if (seqs->seq[i].rseq != NULL) {
                        sseq.seqid = (i << KMERID_BITS | (((slen-j)*2+1) & KMERID_MASK));
                        sseq.seq = seqs->seq[i].rseq + slen - j;
                        subseqs->sub[subseqs->size++] = sseq;
                     }
                  }
                  if (a == maxtau && gapend - match->read_e > SEQ_MINLEN) fprintf(stdout, "[interval %ld]\t(%d,%d)\n", ++cnt, match->read_e+1, gapend-1);
               }
            }
            fprintf(stdout, "[%d] matched %ld out of %d (%.1f%%) nucleotides. (%ldnt will be queried again)\n",i, matched, slen, matched*100.0/slen, slen-matched);
         }
         fprintf(stderr, "aligning... %7d/%7d - alignment done\t[%.3fms]\n", numseqs, numseqs, (clock()-tstart)*1000.0/CLOCKS_PER_SEC);
      }
   }
   

   return 0;
}



int
hitmap_analysis
(
 vstack_t    * hitmap,
 matchlist_t * matchlist,
 int           kmer_size,
 int           readlen
 )
{
   if (matchlist->size < 1) return -1;
   long minv = 0;
   int  min  = 0;

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
      long refend = refloc + (refdir ? (readlen - refkmer) : (refkmer))  * READ_TO_GENOME_RATIO;

      // Mark subseq as used.
      hitmap->val[ref] *= -1;

      // Start streak.
      int streak = 1;

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
         long kmer = (hitmap->val[i] & KMERID_MASK) >> 1;
         long dir = hitmap->val[i] & 1;

         // Distances.
         long r_dist = refkmer - kmer;
         long g_dist = loc - refloc; // Recall that the genome is stored backwards.

         // DEBUG.
         //         fprintf(stdout, "[%ld<->%ld]:\tend=%ld(refkmer=%ld, readlen=%d)\tr_dist=%ld\tg_dist=%ld\tfactor=%f\tkmers=[%ld:%c,%ld:%c]\tloc=[%ld,%ld]\n", refdbg, i, refend, refkmer, readlen, r_dist, g_dist, ((float)g_dist)/r_dist, refkmer, (refdir ? '+' : '-'), kmer, (dir ? '+' : '-'), refloc, loc);

         // Check if the compared sequence is too far away.
         if (loc > refend) {
            if (ref == -1) ref = i;
            break;
         }
      
         // Check if the sequences are placed similarly both in the read and the genome.
         if (r_dist > 0 && dir == refdir && (g_dist < (READ_TO_GENOME_RATIO * r_dist) || (g_dist < MIN_DISTANCE_RATIO && r_dist < MIN_DISTANCE_RATIO))) {
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
         match->hits   = streak;
         match->score  = -1;
         match->repeats = matchlist_new(REPEATS_SIZE);

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
            }
         else {
            match_t * match_min = matchlist->match[min];
            free(match_min->repeats);
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

   return 0;
}


int
map_hits
(
 pstack_t ** hits,
 vstack_t ** hitmap,
 index_t   * index,
 int         tau,
 int         id
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
         if (n_hits < HIT_MAX_LOCI) {
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
               hmap->val[hmap->pos++] = ((index->pos[i] + offset) << KMERID_BITS) | idstamp;
            }
         }
      }
   }
   return 0;
}

int
hitmap_push
(
 vstack_t ** hitmap,
 index_t  *  index,
 long     *  fm_ptr,
 int         id
)
{
   long       idstamp = id & KMERID_MASK;
   vstack_t * hmap = *hitmap;
   long n_hits = fm_ptr[1] - fm_ptr[0] + 1;
   // Filter out too abundant sequences. (To speed up the sorting).
   if (n_hits < HIT_MAX_LOCI) {
      // Realloc hitmap if needed.
      if (hmap->pos + n_hits >= hmap->size) {
         long newsize = hmap->size + n_hits + 1;
         hmap = *hitmap = realloc(hmap, sizeof(vstack_t) + newsize * sizeof(long));
         if (hmap == NULL) {
            fprintf(stderr, "error in 'hitmap_push' (realloc): %s\n", strerror(errno));
            return -1;
         }
         hmap->size = newsize;
      }

      // Copy hits.
      for (long i = fm_ptr[0] ; i <= fm_ptr[1] ; i++) {
         hmap->val[hmap->pos++] = (index->pos[i] << KMERID_BITS) | idstamp;
      }
   }
   return 0;
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
      for (int k = 0; k < slen; k++) if (s.seq[k] > 90) s.seq[k] -= 32;

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
 int slen
)
{
   matchlist_t * list = *listp;
   int current_pos = list->pos;
   matchlist_t * to_combine = matchlist_new(current_pos);

   for (int i = 0; i < current_pos - 1; i++) {
      match_t * match = list->match[i];
      if (match == NULL) continue;
      long max_distance = (slen - match->read_s)*READ_TO_GENOME_RATIO;
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
         if (ref_d <= read_d * READ_TO_GENOME_RATIO || (read_d < MIN_DISTANCE_RATIO && ref_d < MIN_DISTANCE_RATIO)) {
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
         totalnt = match->read_e - match->read_s + 1;
         cmid += totalnt * match->ident;

         // Expand combination.
         for (int j = 0; j < to_combine->pos; j++) {
            match_t * comb = to_combine->match[j];
            if (comb->ref_s < newmatch->ref_s) newmatch->ref_s = comb->ref_s;
            if (comb->ref_e > newmatch->ref_e) newmatch->ref_e = comb->ref_e;
            if (comb->read_e > newmatch->read_e) newmatch->read_e = comb->read_e;
            if (comb->read_s < newmatch->read_s) newmatch->read_s = comb->read_s;
            int readspan = comb->read_e - comb->read_s + 1;
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
 matchlist_t * list
)
{
   // Delete previous repeat record.
   for (int i = 0; i < list->pos; i++) list->match[i]->repeats->pos = 0;

   // Iterate over all the matches and compare overlaps.
   for (int i = 0; i < list->pos; i++) {
      match_t * ref = list->match[i];
      // Add a feedback link. This will be helpful when sorting.
      matchlist_add(&(ref->repeats), ref);
      // Compute the position of the last nucleotide.
      int last_nt = ref->read_s + (int)((ref->read_e - ref->read_s)*(1-REPEAT_OVERLAP)) + 1;
      // Iterate over the other matches.
      for (int j = i+1; j < list->pos; j++) {
         match_t * cmp = list->match[j];
         if (cmp->read_s > last_nt) break;
         // Compute combined sequence span and overlap.
         int span    = hm_max(cmp->read_e, ref->read_e) - hm_min(cmp->read_s, ref->read_s);
         int overlap = hm_min(cmp->read_e, ref->read_e) - hm_max(cmp->read_s, ref->read_s);
         // Check whether the two compared matches share at least REPEAT_OVERLAP.
         if (overlap*1.0/span >= REPEAT_OVERLAP) {
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
 matchlist_t * list
)
{
   // Allocate intervals.
   matchlist_t * interv = matchlist_new(list->pos);

   // Fill read with maximum coverage.
   for (int i = 0; i < list->pos; i++) {
      int append = 1;
      for (int j = 0; j < interv->pos; j++) {
         if (list->match[i]->read_s < interv->match[j]->read_s) {
            if (list->match[i]->read_e - interv->match[j]->read_s > OVERLAP_THR) {
               append = 0;
               break;
            }
         } else if (list->match[i]->read_e > interv->match[j]->read_e) {
            if (interv->match[j]->read_e - list->match[i]->read_s > OVERLAP_THR) {
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

   mergesort_mt(interv->match, interv->pos, sizeof(match_t *), 0, 1, compar_matchstart);   

   return interv;
}

void
recursive_free
(
 mnode_t * node
)
{
   // Recursive call.
   for (int i = 0; i < node->nlinks; i++) recursive_free(node->child[i]);
   // Free current node.
   free(node->child);
   free(node);
}


mnode_t *
recursive_build
(
 mnode_t * node,
 match_t * match
)
{
   int add = (node->nlinks == 0);
   mnode_t * best = node;
   for (int i = 0; i < node->nlinks; i++) {
      mnode_t * child = node->child[i];
      if (child->match->read_e - match->read_s < OVERLAP_THR) {
         mnode_t * cmp = recursive_build(child, match);
         if (cmp->matched > best->matched) best = cmp;
      } else {
         add = 1;
      }
   }      

   if (add) {
      mnode_t * newnode = mnode_new(MNODE_SIZE);
      newnode->parent = node;
      newnode->match = match;
      newnode->matched = node->matched + (match->read_e - match->read_s) - hm_max(0, node->match->read_e - match->read_s);
      // Inline function.
      if (node->nlinks >= node->size) {
         int newsize = node->size * 2;
         node->child = realloc(node->child, newsize * sizeof(mnode_t *));
         if (node->child == NULL) return NULL;
         node->size = newsize;
      }
      // Add new child.
      node->child[node->nlinks++] = newnode;
      
      // Check best combination.
      if (newnode->matched > best->matched) best = newnode;
   }

   return best;
}


mnode_t *
mnode_new
(
 int children
)
{
   if (children < 1) children = 1;
   mnode_t * node = malloc(sizeof(mnode_t));
   mnode_t ** child = malloc(children * sizeof(mnode_t *));
   if (node == NULL) return NULL;
   node->match = NULL;
   node->parent = NULL;
   node->size = children;
   node->nlinks = 0;
   node->matched = 0;
   node->child = child;
   
   return node;
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
      compar = (seq_a->seqid >= seq_b->seqid ? 1 : -1);
   } 
   return compar;
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
compar_matchlen
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   long lena = ma->read_e - ma->read_s;
   long lenb = mb->read_e - mb->read_s;
   if (lenb > lena) return 1;
   else if (lenb < lena) return -1;
   else return (mb->ident > ma->ident ? 1 : -1);
}

int
compar_matchstart
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
   else return (mb->ident > ma->ident ? 1 : -1);
}

int
compar_matchsize
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
   else return (ma->ref_s > mb->ref_s ? 1 : -1);
}
