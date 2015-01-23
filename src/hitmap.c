#include "hitmap.h"

int
hitmap
(
 int          tau,
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
   pstack_t ** hits    = malloc((tau+1) * sizeof(pstack_t *));
   for (int i = 0 ; i <= MAX_TRAIL ; i++) {
      pebbles[i] = new_pstack(PBSTACK_SIZE);
   }
   for (int i = 0 ; i <= tau ; i++) {
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
   int max_aligns = 1000;
   matchlist_t * matches = malloc(sizeof(matchlist_t) + max_aligns * sizeof(match_t));
   if (matches == NULL) return -1;
   matches->size = max_aligns;
   matches->pos  = 0;

   for (int s = 0 ; s < seqs->pos; s += seq_block) {
      // Pre-process sequence block: chunk in subsequences and sort.
      int         numseqs = (seq_block < seqs->pos - s ? seq_block : seqs->pos - s);
      sublist_t * subseqs = process_subseq(seqs->seq+s, numseqs, KMER_SIZE, hitmaps);
      
      fprintf(stderr, "subsequences %d to %d\n", s+1, s+numseqs);

      // Clear hitmaps.
      for (int i = 0 ; i < numseqs ; i++) {
         hitmaps[i]->pos = 0;
      }

      // Poucet variables.
      int start = 0;
      clock_t tstart = clock();

      // Poucet algorithm over the k-mers.
      for (int i = 0; i < subseqs->size; i++) {
         if (1%100 == 0) fprintf(stderr, "mapping ... %7d/%7d\r", i, subseqs->size);
         // Extract subseq.
         sub_t query = subseqs->sub[i];
         int    qlen = KMER_SIZE;

         int trail = 0;
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
         for (int j = 0; j <= tau; j++) {
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
            .tau      = tau,
            .trail    = trail,
            .qlen     = qlen,
            .index    = &index,
            .triep    = &trie,
            .pebbles  = pebbles,
            .hits     = hits
         };

         // Run recursive search from cached pebbles.
         uint row[2*MAXTAU+1];
         uint * nwrow = row + MAXTAU;
         char path[qlen+tau+1];

         for (int p = 0 ; p < pebbles[start]->pos ; p++) {
            // Next pebble.
            pebble_t pebble = pebbles[start]->pebble[p];
            // Compute current alignment from alignment trie.
            int wingsz;
            if (tau > 0) trie_getrow(trie, pebble.rowid >> PATH_SCORE_BITS, pebble.rowid & PATH_SCORE_MASK, &wingsz, nwrow);
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
         map_hits(hits, query.hitmap, &index, tau, query.seqid);

         start = trail;
      }
      fprintf(stderr, "mapping ... %7d/%7d - hitmap built\t[%.3fms]\n", subseqs->size, subseqs->size, (clock()-tstart)*1000.0/CLOCKS_PER_SEC);

      tstart = clock();
      // Sort all hitmaps
      for (int i = 0 ; i < numseqs ; i++) {
         if (i % 100 == 0) fprintf(stderr, "aligning... %7d/%7d\r", i, numseqs);
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
         //radix_sort(hitmaps[i]->val, sortbuf, hitmaps[i]->pos, index.gsize);

         // Find matching loci in hitmaps.
         hitmap_analysis(hitmaps[i], matches, KMER_SIZE, tau, slen);
         if (matches == NULL) continue;

         matchlist_t * significant = malloc(sizeof(matchlist_t) + matches->pos * sizeof(match_t));
         significant->size = matches->pos;
         significant->pos = 0;

         // Smith-Waterman alignment.
         for(long k = 0; k < matches->pos; k++) {
            // Get next match.
            match_t * match = matches->match + k;
            // Extend alignment.
            align_t align_r, align_l;
            char * read = (match->dir ? seqs->seq[i].rseq : seqs->seq[i].seq);

            // Will align both from the start point (backward from start, forward from start+1).
            // The previously matched region will be aligned again,
            // but it's the only way to know precisely the total identity of the read.
            // Note that the genome is stored backwards so all signs are changed.
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, index.gsize - match->ref_s +1)-1;
            align_r = nw_align(read + match->read_s + 1, index.genome + match->ref_s - 1, slen - match->read_s - 1, ALIGN_FORWARD, ALIGN_BACKWARD);
            align_l = nw_align(read + match->read_s, index.genome + match->ref_s, match->read_s, ALIGN_BACKWARD, ALIGN_FORWARD);

            match->score  = align_l.score + align_r.score;
            match->ident  = 1.0 - (match->score*1.0)/(align_l.pathlen + align_r.pathlen);
            match->ref_e  = index.gsize - match->ref_s + 1 + align_r.end; // Align.end is the genome breakpoint.
            match->ref_s  = index.gsize - match->ref_s - align_l.end; // Align.end is the genome breakpoint.
            match->read_e = match->read_s + 1 + align_r.start; // Align.start is the read breakpoint.
            match->read_s = match->read_s - align_l.start;     // Align.start is the read breakpoint.

            // Filter out matches that do not meet the minimum quality.
            if ((match->ref_e - match->ref_s < SEQ_MINLEN) || (match->ident < SEQ_MINID)) continue;
            
            // Correct read start and end points if we are matching the reverse complement.
            if (match->dir) {
               long end = slen - match->read_s;
               match->read_s = slen - match->read_e;
               match->read_e = end;
            }

            // Add to significant match list.
            significant->match[significant->pos++] = *match;
         }

         // Sort mapped regions by significance.
         match_t aux[significant->pos];
         memcpy(aux, significant->match, significant->pos * sizeof(match_t));
         mergesort_match(significant->match, aux, significant->pos, 0);

         // Print results.
         if (significant->pos) fprintf(stdout, "%s\n", seqs->seq[s+i].tag);

         for (long k = 0; k < significant->pos; k++) {
            match_t match = significant->match[k];
            int dir = match.dir;
            long g_start = match.ref_s;
            long g_end   = match.ref_e;
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
            fprintf(stdout, "%d\t%d (%.2f%%)\t%d-%d\t%s:%ld-%ld:%c (%ld)\n",
                    match.hits,
                    match.score,
                    match.ident*100.0,
                    match.read_s, match.read_e,
                    chr->name[chrnum],
                    g_start - chr->start[chrnum]+1,
                    g_end - chr->start[chrnum]+1,
                    dir ? '-' : '+',
                    g_end - g_start + 1);
         }
      }
      fprintf(stderr, "aligning... %7d/%7d - alignment done\t[%.3fms]\n", numseqs, numseqs, (clock()-tstart)*1000.0/CLOCKS_PER_SEC);
   }

   return 0;
}



int
hitmap_analysis
(
 vstack_t    * hitmap,
 matchlist_t * matchlist,
 int           kmer_size,
 int           tau,
 int           readlen
 )
{
   if (matchlist->size < 1) return -1;
   long minv = 0;
   int  min  = 0;

   // Initialize matchlist list.
   matchlist->pos = 0;

   // Initialize match.
   match_t match;
   match.read_s = 0;
   match.read_e = 0;
   match.score  = 0;

   // Find clusters in hitmap.
   long ref = 0;

   // DEBUG.
   //   fprintf(stdout, "Hitmap Analysis:\n\n");

   while (ref >= 0) {
      // Reference position.
      long refloc = hitmap->val[ref] >> KMERID_BITS;
      long refkmer = (hitmap->val[ref] & KMERID_MASK) >> 1;
      long refdir = hitmap->val[ref] & 1;
      long refend = refloc + refkmer * READ_TO_GENOME_RATIO;

      // Mark subseq as used.
      hitmap->val[ref] *= -1;

      // Start streak.
      int streak = 1;

      // Store reference start and read offsets.
      match.ref_e  = refloc;
      match.read_e = refkmer;
      match.dir    = refdir;

      // DEBUG.
      long refdbg = ref;
      
      // Explore hitmap.
      long i = ref;
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
         //         fprintf(stdout, "[%ld<->%ld]:\tr_dist=%ld\tg_dist=%ld\tfactor=%f\tkmers=[%ld:%c,%ld:%c]\tloc=[%ld,%ld]\n", refdbg, i,r_dist, g_dist, ((float)g_dist)/r_dist, refkmer, (refdir ? '+' : '-'), kmer, (dir ? '+' : '-'), refloc, loc);

         // Check if the compared sequence is too far away.
         if (loc > refend) {
            if (ref == -1) ref = i;
            break;
         }
      
         // Check if the sequences are placed similarly both in the read and the genome.
         if (r_dist > 0 && dir == refdir && g_dist < READ_TO_GENOME_RATIO * r_dist) {
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
         // Recompute the position of the last subseq in the streak.
         match.ref_s = streakloc;
         match.read_s = streakkmer;
         match.hits = streak;
         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) matchlist->match[matchlist->pos++] = match;
         else matchlist->match[min] = match;
               
         // Find the minimum that will be substituted next.
         if (matchlist->pos == matchlist->size) {
            minv = matchlist->match[0].hits;
            for (int j = 1 ; j < matchlist->pos; j++) {
               if (minv > matchlist->match[j].hits) {
                  min  = j;
                  minv = matchlist->match[j].hits;
               }
            }
         }
      }

   }
   return 0;
}

/*
int
hitmap_analysis
(
 vstack_t    * hitmap,
 matchlist_t * matchlist,
 int           kmer_size,
 int           tau
 )
{
   if (matchlist->size < 1) return -1;
   long minv = 0;
   int  min  = 0;

   // Initialize matchlist list.
   matchlist->pos = 0;

   // Initialize match.
   match_t match;
   match.read_s = 0;
   match.read_e = 0;
   match.score  = 0;

   // Find clusters in hitmap.
   int  streak = 1;
   // First subsequence.
   long loc1, loc2 = hitmap->val[0] >> KMERID_BITS;
   long kmer1, kmer2 = (hitmap->val[0] & KMERID_MASK) >> 1;
   long dir1, dir2 = hitmap->val[0] & 1;

   // DEBUG.
   fprintf(stdout, "Hitmap Analysis:\n\n");

   for (int i = 1; i < hitmap->pos - 1; i++) {
      // Swap subsequences 2 and 1.
      loc1 = loc2;
      kmer1 = kmer2;
      dir1 = dir2;

      // Get new subsequence.
      loc2 = hitmap->val[i+1] >> KMERID_BITS;
      kmer2 = (hitmap->val[i+1] & KMERID_MASK) >> 1;
      dir2 = hitmap->val[i+1] & 1;
      long r_dist = kmer1 - kmer2;
      long g_dist = loc2 - loc1;
      
      // Start streak.
      if (streak == 1) match.ref_s = max(0, loc1 - WINDOW_SIZE);

      // DEBUG.
      fprintf(stdout, "[%d<->%d]:\tr_dist=%ld\tg_dist=%ld\tfactor=%f\tkmers=[%ld:%c,%ld:%c]\tloc=[%ld,%ld]\n", i, i+1,r_dist, g_dist, ((float)g_dist)/r_dist, kmer1, (dir1 ? '+' : '-'), kmer2, (dir2 ? '+' : '-'), loc1, loc2);

      if (r_dist > 0 && dir1 == dir2 && g_dist < READ_TO_GENOME_RATIO * r_dist) streak++;
      else {
         match.ref_e = loc1 + kmer_size + WINDOW_SIZE;
         match.hits = streak;
         // Check significance.
         if (streak > minv) {
            // Append if list not yet full, replace the minimum value otherwise.
            if (matchlist->pos < matchlist->size) matchlist->match[matchlist->pos++] = match;
            else matchlist->match[min] = match;
            
            // Find the minimum that will be substituted next.
            if (matchlist->pos == matchlist->size) {
               minv = matchlist->match[0].hits;
               for (int j = 1 ; j < matchlist->pos; j++) {
                  if (minv > matchlist->match[j].hits) {
                     min  = j;
                     minv = matchlist->match[j].hits;
                  }
               }
            }
         }

         // Reset streak.
         streak = 1;
      }
   }
   return 0;
}
*/

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
   if (n_hits < EXACT_MAX_LOCI) {
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
   // Sort subsequences.
   mergesort_mt(list->sub, list->size, sizeof(sub_t), k, 1, compar_seqsort);

   return list;
}

void
mergesort_match
(
 match_t * data,
 match_t * aux,
 int size,
 int b
 )
// SYNOPSIS:
//   Recursive part of 'seqsort'.
//
// ARGUMENTS:
//   args: a sortargs_t struct.
//
// RETURN:
//   
//
// SIDE EFFECTS:
//   Sorts the array of 'seq_t' specified in 'args'.
{
   
   if (size < 2) return;

   // Next level params.
   int newsize = size/2;
   int newsize2 = newsize + size%2;
   match_t * data2 = data + newsize;
   match_t * aux2 = aux + newsize;
   mergesort_match(data, aux, newsize, (b+1)%2);
   mergesort_match(data2, aux2, newsize2, (b+1)%2);


   // Separate data and buffer (b specifies which is buffer).
   match_t * l = (b ? data : aux);
   match_t * r = (b ? data2 : aux2);
   match_t * buf = (b ? aux : data);

   int i = 0;
   int j = 0;
   int idx = 0;

   // Merge sets
   while (idx < size) {
      // Right buffer is exhausted. Copy left buffer...
      if (j == newsize2) {
         memcpy(buf+idx, l+i, (newsize-i) * sizeof(match_t));
         break;
      }
      // ... or vice versa.
      if (i == newsize) {
         memcpy(buf+idx, r+j, (newsize2-j) * sizeof(match_t));
         break;
      }
      // Do the comparison.
      /*
      if (l[i].hits > r[j].hits || (l[i].hits == r[j].hits && l[i].score > r[j].score)) buf[idx++] = l[i++];
      else                                                                              buf[idx++] = r[j++];
      */
      // Sort by identity.
      if (l[i].ident > r[j].ident) buf[idx++] = l[i++];
      else buf[idx++] = r[j++];
   }
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
