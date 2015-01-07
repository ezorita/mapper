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
   for (int i = 0 ; i < tau+1 ; i++) {
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
   int seq_block = 1000;
   // Allocate hitmaps.
   vstack_t ** hitmaps = malloc(seq_block * sizeof(vstack_t *));
   for (int i = 0 ; i < seq_block ; i++) {
      hitmaps[i] = new_stack(HITMAP_SIZE);
      if (hitmaps[i] == NULL) return -1;
   }
   // Allocate loci stack.
   int max_match = 25;
   matchlist_t * matches = malloc(sizeof(matchlist_t) + max_match * sizeof(match_t));
   if (matches == NULL) return -1;
   matches->size = max_match;
   matches->pos  = 0;

   for (int s = 0 ; s < seqs->pos; s += seq_block) {
      // Pre-process sequence block: chunk in subsequences and sort.
      int         numseqs = (seq_block < seqs->pos - s ? seq_block : seqs->pos - s);
      sublist_t * subseqs = process_subseq(seqs->seq+s, numseqs, KMER_SIZE, hitmaps);

      // Clear hitmaps.
      for (int i = 0 ; i < numseqs ; i++) {
         hitmaps[i]->pos = 0;
      }
         
      // Poucet variables.
      int start = 0;
      char * lastseq = NULL;

      // Poucet algorithm over the k-mers.
      for (int i = 0; i < subseqs->size; i++) {
         // Extract subseq.
         sub_t query = subseqs->sub[i];
         int    qlen = KMER_SIZE;

         // Copy hits for repeated sequences.
         if (lastseq != NULL && strcmp(lastseq, query.seq) == 0) {
            // DEBUG.
            fprintf(stdout, "[repeated seq]\n");
            map_hits(hits, query.hitmap, &index, tau, i);
            continue;
         }

         // Compute trail depth.
         int trail = 0;
         if (i < subseqs->size - 1) {
            // Compute trail wrt the next different sequence.
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

         // DEBUG.
         //fprintf(stdout, "seq: %s", query.seq);
         for (int p = 0 ; p < pebbles[start]->pos ; p++) {
            // Next pebble.
            pebble_t pebble = pebbles[start]->pebble[p];
            // Compute current alignment from alignment trie.
            int wingsz;
            trie_getrow(trie, pebble.rowid >> PATH_SCORE_BITS, pebble.rowid & PATH_SCORE_MASK, &wingsz, nwrow);
            // Recover the current path.
            long gpos = index.pos[pebble.sp];
            for (int j = 0; j < start; j++) {
               path[start-j] = translate[(int)index.genome[gpos+j]];
            }
            poucet(pebble.sp, pebble.ep, wingsz, nwrow, start + 1, path, &arg);
         }

         // Map hits.
         map_hits(hits, query.hitmap, &index, tau, i);

         // Sort hitmap.
         start = trail;
      }

      // Sort all hitmaps
      for (int i = 0 ; i < numseqs ; i++) {
         // Realloc sort buffer if necessary.
         if (sortbuf_size < hitmaps[i]->pos) {
            sortbuf = realloc(sortbuf, hitmaps[i]->pos * sizeof(long));
            if (sortbuf == NULL) return -1;
         }
         // Merge-sort loci.
         memcpy(sortbuf, hitmaps[i]->val, hitmaps[i]->pos * sizeof(long));
         mergesort_long(hitmaps[i]->val, sortbuf, hitmaps[i]->pos, 0);
         //radix_sort(hitmaps[i]->val, sortbuf, hitmaps[i]->pos, index.gsize);
         for (int n = 0; n < hitmaps[i]->pos-1; n++) {
            if (hitmaps[i]->val[n] > hitmaps[i]->val[n+1])
               fprintf(stdout, "sorting error index %d\n", n);
         }
         // Find matching loci in hitmaps.
         hitmap_analysis(hitmaps[i], matches, KMER_SIZE, tau, WINDOW_SIZE);
         if (matches == NULL) continue;

         // Smith-Waterman align of hits.
         for(long k = 0; k < matches->pos; k++) {
            match_t * match = matches->match + k;
            long rlen  = match->ref_e - match->ref_s + 1;               
            // Copy and reverse reference sequence.
            char ref[rlen];
            for (int i = 0 ; i < rlen; i++) ref[i] = index.genome[match->ref_e - i];

            int slen = strlen(seqs->seq[i].seq);
            align_t salign = sw_align(seqs->seq[i].seq, slen, ref, rlen);
            align_t ralign = sw_align(seqs->seq[i].rseq, slen, ref, rlen);
            align_t best = (salign.score > ralign.score ? salign : ralign);
            match->score  = best.score;
            match->read_s = best.start;
            match->read_e = best.max;
         }

         // Sort mapped regions by significance.
         match_t aux[matches->pos];
         memcpy(aux, matches->match, matches->pos * sizeof(match_t));
         mergesort_match(matches->match, aux, matches->pos, 0);

         // Print results.
         for (long k = 0; k < matches->pos; k++) {
            match_t match = matches->match[k];
            long g_start = index.gsize - match.ref_e;
            long g_end   = index.gsize - match.ref_s;
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, g_start+1)-1;
            fprintf(stdout, "%s\t%d\t%d\t%d-%d\t%s:%ld-%ld (%ld)\n",
                    seqs->seq[s+i].tag,
                    match.hits,
                    match.score,
                    match.read_s, match.read_e,
                    chr->name[chrnum],
                    g_start - chr->start[chrnum]+1,
                    g_end - chr->start[chrnum]+1,
                    g_end - g_start + 1);
         }
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
 int           tau,
 int           maxdist
)
{
   if (matchlist->size < 1) return -1;
   long minv = 0;
   int  min  = 0;
   int  mindist = kmer_size - tau;
   
   // Initialize matchlist list.
   matchlist->pos = 0;

   // Initialize match.
   match_t match;
   match.read_s = 0;
   match.read_e = 0;
   match.score  = 0;
   
   // Find clusters in hitmap.
   int   streak = 0;
   long id_mask = (1 << SUBSEQID_BITS) - 1;
   long loc1, loc2 = hitmap->val[0] >> SUBSEQID_BITS;
   long sid1, sid2 = hitmap->val[0] & id_mask;
   for (int i = 1; i < hitmap->pos - 1; i++) {
      loc1 = loc2;
      sid1 = sid2;
      loc2 = hitmap->val[i+1] >> SUBSEQID_BITS;
      sid2 = hitmap->val[i+1] & id_mask;
      long diff = loc2 - loc1;
      if (diff < maxdist) {
         if (diff > mindist && sid1 != sid2) {
            if (!streak) {
               match.ref_s = max(0, loc1 - WINDOW_SIZE);
               streak = 2;
            } else streak++;
         }
      } else {
         if (streak) {
            match.ref_e = loc1 + kmer_size + WINDOW_SIZE;
            match.hits = streak;
            // Check significance.
            if (streak > minv) {
               matchlist->match[min] = match;
               // Extend matchlist list if not yet full.
               if (matchlist->pos == min && matchlist->pos < matchlist->size) min = ++matchlist->pos;
               // Find minimum.
               else {
                  min = 0;
                  minv = matchlist->match[0].hits;
                  for (int j = 1 ; j < matchlist->pos; j++) {
                     if (minv > matchlist->match[j].hits) {
                        min  = j;
                        minv = matchlist->match[j].hits;
                     }
                  }
               }
               // End find minimum.
            }
            streak = 0;
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
   long       idstamp = id & ((1 << SUBSEQID_BITS) - 1);

   // Push hits to hitmap.
   for (int a = 0 ; a < tau ; a++) {
      for (int h = 0; h < hits[a]->pos; h++) {
         pebble_t hit    = hits[a]->pebble[h];
         long     n_hits = hit.ep - hit.sp + 1;
         // Filter out too abundant sequences. (To speed up the sorting).
         if (n_hits < HIT_MAX_LOCI) {
            // Realloc hitmap if needed.
            if (hmap->pos + n_hits >= hmap->size) {
               long newsize = hmap->size + n_hits + 1;
               hmap = *hitmap = realloc(hmap, sizeof(vstack_t) + newsize * sizeof(long));
               if (hmap == NULL) {
                  fprintf(stderr, "error in 'push' (realloc): %s\n", strerror(errno));
                  return -1;
               }
               hmap->size = newsize;
            }

            // Copy hits.
            for (long i = hit.sp; i <= hit.ep ; i++) {
               hmap->val[hmap->pos++] = (index->pos[i] << SUBSEQID_BITS) | idstamp;
            }
         }
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
{
   if (numseqs < 1) return NULL;
   int    rev   = (seqs[0].rseq != NULL);
   char * seq   = seqs[0].seq;
   int    lsize = (strlen(seq)/k+1)*(1+rev)*numseqs;
   int    lpos  = 0;
   sublist_t * list = malloc(sizeof(sublist_t) + lsize * sizeof(sub_t));
   if (list == NULL) return NULL;
   
   for (int i = 0; i < numseqs; i++) {
      seq_t s  = seqs[i];
      int slen = strlen(s.seq);
      int last = (slen%k > LAST_THRESHOLD);
      int subs = slen/k + last;
      int    p = 0;

      for (int j = 0; j < subs; j++) {
         // Copy the last sub if the remaining nucleotides > LAST_THRESHOLD.
         if (j == subs-1 && last) p = slen - k;
         // Realloc full list.
         if (lpos + rev >= lsize) {
            lsize *=2;
            list = realloc(list, sizeof(sublist_t) + lsize * sizeof(sub_t));
            if (list == NULL) return NULL;
         }
         // Generate subseq and save to list.
         sub_t subseq = { .seq = s.seq + p, .hitmap = hitmaps + i};
         list->sub[lpos++] = subseq;
         // Insert the reverse complement as well.
         if (rev) {
            subseq.seq = s.rseq + p;
            list->sub[lpos++] = subseq;
         }
         // Jump k nucleotides.
         p += k;
      }
   }
   // Realloc list.
   list->size = lpos;
   list = realloc(list, sizeof(sublist_t) + lpos * sizeof(sub_t));

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
      if (l[i].hits > r[j].hits || (l[i].hits == r[j].hits && l[i].score > r[j].score)) buf[idx++] = l[i++];
      else                                                                              buf[idx++] = r[j++];
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
   
   return strncmp(seq_a->seq, seq_b->seq, val);
}
