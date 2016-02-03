#include "seed.h"

int
check_kmer
(
 uint8_t * seq,
 uint8_t * query,
 int       slen
)
{
   int strand = 0;
   for (int i = 0; i < slen; i++) {
      // Check 'N'
      if (seq[i] > 3) return -1;
      // Compare forward and reverse strands.
      if (!strand) {
         if (3 - seq[slen-1-i] < seq[i]) strand = -1;
         else if (3 - seq[slen-1-i] > seq[i]) strand = 1;
      }
   }
   if (strand < 0)
      for (int i = 0; i < slen; i++) query[i] = 3-seq[slen-i-1];
   else
      memcpy(query, seq, slen);
   return 0;
}

int
find_uniq_seeds
(
 int            slen,
 int            k,
 uint8_t      * query,
 index_t      * index,
 seedstack_t ** unique,
 seedstack_t ** good
)
{
   uint8_t * q = malloc(k);
   htable_t * htable = (htable_t *)index->seeds;
   for (int i = 0; i <= slen - k; i++) {
      if (check_kmer(query+i, q, 25)) continue;
      uint64_t key = XXH64(q, k, 0);
      uint8_t score = htable_get(htable, key);
      if (score == 1) {
         bwpos_t pos = {0,0,index->size-1};
         for (int j = i+k-1; j >= i && pos.ep >= pos.sp; j--) suffix_extend(query[j],pos,&pos,index);
         if (pos.ep == pos.sp) {
            uint64_t locus = get_sa(pos.sp, index->sa, index->sa_bits);
            if ((index->repeats[locus>>3] >> (locus&7))&1) {
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = i, .ref_pos = pos};
               seedstack_push(seed, unique);
            }
         }
      } else if (score == 2) {
         bwpos_t pos = {0,0,index->size-1};
         for (int j = i+k-1; j >= i && pos.ep >= pos.sp; j--) suffix_extend(query[j],pos,&pos,index);
         if (pos.ep >= pos.sp && pos.ep - pos.sp + 1 <= 20) {
            seed_t seed = (seed_t) {.bulk = 0, .qry_pos = i, .ref_pos = pos};
            seedstack_push(seed, good);
         }
      }
   }
   free(q);
   return 0;
}

seed_t
find_uniq_seed
(
 int            beg,
 int            slen,
 int            k,
 uint8_t      * query,
 index_t      * index,
 int          * zero_cnt
)
{
   uint8_t * q = malloc(k);
   htable_t * htable = (htable_t *)index->seeds;
   //   fprintf(stderr,"unique search: ");
   for (int i = beg; i <= slen - k; i++) {
      if (check_kmer(query+i, q, k)) continue;
      uint64_t key = XXH64(q, k, 0);
      int value = htable_get(htable, key);
      //      fprintf(stderr," %d",value);
      if (value == 1) {
         bwpos_t pos = {0,0,index->size-1};
         for (int j = i+k-1; j >= i && pos.ep >= pos.sp; j--) suffix_extend(query[j],pos,&pos,index);
         //         fprintf(stderr,(pos.ep < pos.sp ? "i" : pos.ep > pos.sp ? "m" : ""));
         if (pos.ep == pos.sp) {
            uint64_t locus = get_sa(pos.sp, index->sa, index->sa_bits);
            //            fprintf(stderr,"locus: %ld (%ld)\n", locus, index->size - 1 - locus - k);
            if (locus >= index->size/2) locus = index->size - 1 - locus - k;
            int unique = (index->repeats[locus>>3] >> (locus&7))&1;
            //            if (!unique) fprintf(stderr,"a");
            // DEBUG.
            /*
            pathstack_t * pstack = pathstack_new(PATHSTACK_DEF_SIZE);
            blocksearch(q, 25, 1, index, &pstack);
            if (pstack->pos == 1 && pstack->path[0].pos.sz == 1 && !unique)
               fprintf(stderr, "annotation error, false negative at seq[%d]\n",i);
            if (unique && pstack->pos == 1 && pstack->path[0].pos.sz == 1) {
            */
            if (unique) {
               //               fprintf(stderr," correct\n");
               free(q);
               return (seed_t) {.bulk = 0, .qry_pos = i, .ref_pos = pos};
            }
            //            free(pstack);
         }
      } else if (value == 0) {
         *zero_cnt += 1;
      }
   }
   //   fprintf(stderr," not found\n");
   free(q);
   return (seed_t) {.bulk = 1, .qry_pos = -1};
}

int
find_thr_seed
(
 int            beg,
 int            slen,
 int            k,
 uint8_t      * query,
 index_t      * index
)
{
   uint8_t * q = malloc(k);
   htable_t * htable = (htable_t *)index->seeds;
   for (int i = beg; i <= slen - k; i++) {
      if (check_kmer(query+i, q, k)) continue;
      uint64_t key = XXH64(q, k, 0);
      // Test.
      int value = htable_get(htable, key);
      if (value == 0 || value == 2) {
      //      if (htable_get(htable, key) == 2) {
         free(q);
         return i;
      }
   }
   free(q);
   return -1;
}

seed_t
seed_locally
(
  char    * seq,
  bwpos_t   pos,
  int       mm,
  index_t * index
)
// SYNOPSIS:
// Find seed locally by extend/shrink.
{

   bwpos_t newpos = {0};

   for (int i = mm-1 ; i > -1 ; i--) {

      suffix_extend(translate[(uint8_t)seq[i]], pos, &newpos, index);

      // Backtrack as nuch as needed.
      while (newpos.ep < newpos.sp) {
         suffix_shrink(pos, &pos, index);
         suffix_extend(translate[(uint8_t)seq[i]], pos, &newpos, index);
      }

      // Tail must not go further than initial mismatch position.
      if (i + newpos.depth - 1 < mm) break;
      if (newpos.depth > 24) {
         // XXX in this crapotype, bulk is just used to say that the seed XXX
         // XXX is valid and it can be used for alignment...              XXX
         return (seed_t) {.bulk = 1, .qry_pos = i+1, .ref_pos = newpos};
      }

      pos = newpos;

   }

   // Nothing found.
   return (seed_t) {0};

}

int
seed_the_right_way
(
 char         * seq,
 seedstack_t ** stack,
 index_t      * index
)
{

   bwpos_t pos = {.depth = 0, .sp = 0, .ep = index->size-1};

   for (int i = strlen(seq)-1 ; i > -1 ; i--) {

      bwpos_t newpos;
      uint8_t c = translate[(uint8_t)seq[i]];

      for (uint8_t a = 0 ; a < 4 ; a++) {
         // Do not add the nucleotide present in the sequence.
         if (a == c) continue;

         // Add a different nucleotide and seed locally.
         suffix_extend(a, pos, &newpos, index);
         if (newpos.ep < newpos.sp) {
            bwpos_t tmp = pos;
            do {
               suffix_shrink(tmp, &tmp, index);
               suffix_extend(a, tmp, &newpos, index);
            }
            while (newpos.ep < newpos.sp);
         }

         seed_t seed = seed_locally(seq, newpos, i, index);
         // XXX See comment above. XXX
         if (seed.bulk) {
            *stack[0]->seed = seed;
            (*stack)->pos = 1;
            return 0;
         }

      }

      // Add the nucleotide present in the sequence.
      suffix_extend(c, pos, &newpos, index);

      // Backtrack as much as needed.
      while (newpos.ep < newpos.sp) {
         suffix_shrink(pos, &pos, index);
         suffix_extend(translate[(uint8_t)seq[i]], pos, &newpos, index);
      }

      pos = newpos;

   }

   // Nothing found.
   return 0;

}

int
count_mismatch
(
 char         * seq,
 size_t         slen,
 index_t      * index
)
{
   int count = 0;
   bwpos_t pos = {.sp = 0, .ep = index->size-1, .depth = 0};
   for (int i = slen-1; i >= 0; i--) {
      int nt = translate[(uint8_t)seq[i]];
      bwpos_t next;
      for (int j = 0; j < 4; j++) {
          if (j == nt) {
            suffix_extend(j, pos, &next, index);
         } else {
            bwpos_t tmp;
            suffix_extend(j, pos, &tmp, index);
            for (int k = i-1; k >= 0 && tmp.ep >= tmp.sp; k--) {
               suffix_extend(translate[(uint8_t)seq[k]], tmp, &tmp, index);
            }
            if (tmp.ep >= tmp.sp) {
               count += tmp.ep - tmp.sp + 1;
            }
         }
      }
      pos = next;
   }
   if (pos.ep >= pos.sp) count += pos.ep - pos.sp + 1;

   return count;
}

int
find_mismatch
(
 char         * seq,
 size_t         slen,
 int            qrypos,
 seedstack_t ** stack,
 index_t      * index
)
{
   bwpos_t pos = {.sp = 0, .ep = index->size-1, .depth = 0};
   for (int i = slen-1; i >= 0; i--) {
      int nt = translate[(uint8_t)seq[i]];
      bwpos_t next;
      for (int j = 0; j < 4; j++) {
          if (j == nt) {
            suffix_extend(j, pos, &next, index);
         } else {
            bwpos_t tmp;
            suffix_extend(j, pos, &tmp, index);
            for (int k = i-1; k >= 0 && tmp.ep >= tmp.sp; k--) {
               suffix_extend(translate[(uint8_t)seq[k]], tmp, &tmp, index);
            }
            if (tmp.ep >= tmp.sp) {
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = tmp};
               seedstack_push(seed, stack);
            }
         }
      }
      pos = next;
   }
   if (pos.ep >= pos.sp) {
      seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = pos};
      seedstack_push(seed, stack);
   }

   return 0;
}

int
seed_mismatch
(
 char         * seq,
 size_t         slen,
 int            qrypos,
 seedstack_t ** stack,
 index_t      * index
)
{
   int c = slen/2;
   fmdpos_t pos = {0,0,index->size};
   // forward extend.
   //  forward extend exact until c-1.
   int i = 0;
   for (; i < c; i++) {
      pos = extend_fw(translate[(int)seq[i]],pos,index);
   }
   //  mismatch all the next positions.
   for (; i < slen; i++) {
      int nt = translate[(int)seq[i]];
      fmdpos_t newpos[5];
      extend_fw_all(pos, newpos, index);
      for (int j = 0; j < 4; j++) {
         if (j == nt) {
            // Extend perfect seed.
            pos = newpos[j];
         } else {
            // Extend mismatched seed.
            for (int k = i+1; k < slen && newpos[j].sz; k++)
               newpos[j] = extend_fw(translate[(int)seq[k]],newpos[j],index);
            if (newpos[j].sz) {
               // Save seed.
               bwpos_t refpos = {.sp = newpos[j].fp, .ep = newpos[j].fp + newpos[j].sz - 1, .depth = slen};
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = refpos};
               seedstack_push(seed, stack);
            }
         }
      }
   }
   // reverse extend.
   pos = (fmdpos_t){0,0,index->size};
   i = slen-1;
   for (; i >= c; i--) {
      pos = extend_bw(translate[(int)seq[i]],pos,index);
   }
   //  mismatch all the next positions.
   for (; i >= 0; i--) {
      int nt = translate[(int)seq[i]];
      fmdpos_t newpos[5];
      extend_bw_all(pos, newpos, index);
      for (int j = 0; j < 4; j++) {
         if (j == nt) {
            // Extend perfect seed.
            pos = newpos[j];
         } else {
            // Extend mismatched seed.
            for (int k = i-1; k >= 0 && newpos[j].sz; k++)
               newpos[j] = extend_bw(translate[(int)seq[k]],newpos[j],index);
            if (newpos[j].sz) {
               // Save seed.
               bwpos_t refpos = {.sp = newpos[j].fp, .ep = newpos[j].fp + newpos[j].sz - 1, .depth = slen};
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = qrypos, .ref_pos = refpos};
               seedstack_push(seed, stack);
            }
         }
      }
   }
   return 0;
}

int
seed_block_all
(
 char *         seq,
 size_t         slen,
 seedstack_t ** stack,
 index_t      * index,
 int            blen,
 int            fw
)
{
   // Wing length.
   // Convert block length to even number.
   int c = (blen + (fw > 0)) / 2;
   int v = blen - c;
   // Allocate seed cache.
   bwpos_t cache[100];
   // Seed:
   // - pre-compute constant part.
   int e = slen-1;
   int i = slen-1;
   bwpos_t pos = (bwpos_t) {.sp = 0, .ep = index->size-1, .depth = 0};
   while (i >= v) {
      bwpos_t newpos;
      suffix_extend(translate[(uint8_t)seq[i]], pos, &newpos, index);
      if (newpos.sp <= newpos.ep) {
         // If got length c, store and shrink.
         if (newpos.depth == c) {
            cache[e] = newpos;
            bwpos_t tmp;
            suffix_shrink(newpos, &tmp, index);
            // This is necessary to avoid shrinking more than 1 nt.
            if (tmp.depth < c-1) {
               newpos.depth--;
            } else {
               newpos = tmp;
            }
            e--;
         }
         pos = newpos;
         i--;
      } else {
         cache[e--] = newpos;
         bwpos_t tmp;
         int depth = pos.depth;
         suffix_shrink(pos, &tmp, index);
         if (tmp.depth < depth-1) {
            pos.depth--;
         } else {
            pos = tmp;
         }
      }
   }

   // - seed for all mismatched positions (i).
   for (int i = slen-1-c; i >= 0; i--) {
      int beg = min(i+blen, slen-1);
      int end = max(i+c, blen-1);
      int nt  = translate[(uint8_t)seq[i]];
      // Interval of this seed (p).
      for (int p = beg; p >= end; p--) {
         bwpos_t current = cache[p];
         // Continue on broken seeds.
         if (current.depth != p - i) continue;
         /*
         else if (current.depth == blen && current.ep - current.sp + 1 > 0 && fw == 1) {
            // Save seed.
            seed_t seed = (seed_t) {.bulk = 0, .qry_pos = p-blen+1, .ref_pos = current};
            seedstack_push(seed, stack);
         }
         */
         // Extend with the 4 possible mutations (m).
         for (int m = 0; m < 4; m++) {
            bwpos_t mpos = current;
            suffix_extend(m, mpos, &mpos, index);
            if (mpos.ep < mpos.sp) continue;
            // Cache extension for the next mismatch.
            if (m == nt) {
               cache[p] = mpos;
               continue;
            } 
            int k = i-1;
            // Extend seed until block length.
            while (mpos.ep >= mpos.sp && mpos.depth < blen)
               suffix_extend(translate[(uint8_t)seq[k--]], mpos, &mpos, index);
            // Save seed if seed exists.
            if (mpos.ep >= mpos.sp) {
               // Save seed.
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = p-blen+1, .ref_pos = mpos};
               seedstack_push(seed, stack);
            }
         }
      }
   }

//      for (int i = 0; i < c; i++) {
//         pos = extend_fw(translate[(uint8_t)seq[b+i]], pos, index);
//         if (pos.sz == 0) break;
//      }    
//seed_t seed = (seed_t) {.bulk = 0, .qry_pos = b, .ref_pos = (bwpos_t) {.sp = mpos.fp, .ep = mpos.fp+mpos.sz-1, .depth = blen}};


   return 0;
}  


char *
rev_comp
(
 char * seq
)
{
   int64_t slen = strlen(seq);
   char * rseq = malloc(slen+1);
   for (int64_t i = 0; i < slen; i++)
      rseq[slen-1-i] = revcomp[(uint8_t)seq[i]];
   rseq[slen] = 0;
   return rseq;
}


int
seed_thr
(
 char         * seq,
 size_t         slen,
 seedstack_t ** stack,
 seedopt_t      opt,
 index_t      * index
)
{
   int32_t i = slen-1;
   uint64_t new_loci;
   // Start position.
   bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
   while (i >= 0) {
      int nt = translate[(uint8_t)seq[i]];
      if (nt > 3) {
         pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
         i--;
         continue;
      }
      bwpos_t newpos;
      // Extend suffix (Backward search).
      suffix_extend(nt, pos, &newpos, index);
      // Count loci.
      new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);
      // Seed found, save and shrink suffix.
      if (new_loci == 0) {
         // Shrink suffix once.
         suffix_shrink(pos, &newpos, index);
         new_loci = newpos.ep - newpos.sp + 1;
      } else if (new_loci <= opt.thr_seed && pos.depth >= opt.min_len) {
         // Build seed.
         int bulk = 0;//newpos.depth < opt.min_len;
         seed_t seed = (seed_t) {.bulk = bulk, .qry_pos = i, .ref_pos = newpos};
         // Save seed.
         seedstack_push(seed, stack);

         // Shrink previous suffix (it had loci_count > thr).
         suffix_shrink(newpos, &newpos, index);
         // Check Complementary matches.
         if (newpos.depth > opt.min_len && newpos.ep - newpos.sp + 1 <= opt.min_loci) {
            bwpos_t tmp;
            do {
               tmp = newpos;
               suffix_shrink(newpos, &newpos, index);
            } while (newpos.depth > opt.min_len && newpos.ep - newpos.sp + 1 <= opt.min_loci);
            // Save seed if it bears new loci.
            if (tmp.ep - tmp.sp + 1 > new_loci) {
               seed_t seed = (seed_t) {.bulk = 0, .qry_pos = i, .ref_pos = tmp};
               seedstack_push(seed, stack);
            }
         }
         // Extend.
         i--;
      } else {
         i--;
      }
      pos = newpos;
   }

   return 0;
}


int
seed_mem
(
 char         * seq,
 int32_t         beg,
 int32_t         end,
 seedstack_t ** stack,
 seedopt_t      opt,
 index_t      * index
)
{

   int32_t i;
   int32_t last_qry_pos = end;
   uint64_t new_loci = index->size;//, last_loci;

   int maxseedlength = 0;
   int seed_pos = (*stack)->pos;
   
   // Start position.
   bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};

   for (i = end-1 ; i >= beg ; i--) {
      bwpos_t newpos;
      int nt = translate[(uint8_t)seq[i]];
      // Extend suffix (Backward search).
      suffix_extend(nt, pos, &newpos, index);
      // Count loci.
      //last_loci = new_loci;
      new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);

      if (new_loci > opt.min_loci) {
         pos = newpos;
         continue;
      }

      // Seed found.
      int seedlength = last_qry_pos - i - 1;

      if (seedlength > 24) {
         if (seedlength == maxseedlength) {
            seed_t seed = (seed_t) {.bulk = 0, .qry_pos = i+1, .ref_pos = pos};
            seedstack_push(seed, stack);
         }
         if (seedlength > maxseedlength) {
            maxseedlength = seedlength;
            (*stack)->seed[seed_pos].bulk = 0;
            (*stack)->seed[seed_pos].qry_pos = i+1;
            (*stack)->seed[seed_pos].ref_pos = pos;
            (*stack)->pos = seed_pos+1;
         }
      }

      if (nt == 4) {
         // Reset 'pos' and jump over N.
         last_qry_pos = i+1;
         pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
      } else if (new_loci < 1){
         do {
            int origdepth = pos.depth;
            suffix_shrink(pos, &pos, index);
            new_loci = pos.ep - pos.sp + 1;
            last_qry_pos -= (origdepth - pos.depth);
         } while (new_loci < opt.min_loci);
         // Redo the loop with updated 'pos'.
         // Quite ugly by the way.
         i++;
      }
   }

   if (last_qry_pos > 24) {
      if (last_qry_pos == maxseedlength) {
         seed_t seed = (seed_t) {.bulk = 0, .qry_pos = 0, .ref_pos = pos};
         seedstack_push(seed, stack);
      }
      if (last_qry_pos > maxseedlength) {
         (*stack)->seed[seed_pos].bulk = 0;
         (*stack)->seed[seed_pos].qry_pos = 0;
         (*stack)->seed[seed_pos].ref_pos = pos;
         (*stack)->pos = seed_pos+1;
      }
   }

   return 0;
}


hit_t *
compute_hits
(
 seedstack_t * seeds,
 index_t     * index,
 uint64_t    * hit_cnt
)
{
   // Count loci.
   uint64_t cnt = 0;
   for (size_t i = 0; i < seeds->pos; i++) cnt += seeds->seed[i].ref_pos.ep - seeds->seed[i].ref_pos.sp + 1;
   *hit_cnt = cnt;

   // Find loci in Suffix Array.
   hit_t * hits = malloc(cnt * sizeof(hit_t));
   size_t l = 0;
   for (size_t i = 0; i < seeds->pos; i++) {
      seed_t seed = seeds->seed[i];
      for (size_t j = seed.ref_pos.sp; j <= seed.ref_pos.ep; j++) {
         int64_t loc = get_sa(j,index->sa,index->sa_bits);
         hits[l++] = (hit_t) {.locus = loc, .qrypos = seed.qry_pos, .depth = seed.ref_pos.depth, .bulk = seed.ref_pos.ep - seed.ref_pos.sp + 1};
      }
   }
   
   return hits;

}

int
chain_seeds
(
 seedstack_t * seeds,
 seedstack_t * rseeds,
 int           seqlen,
 matchlist_t * matchlist,
 index_t     * index,
 int           dist_accept,
 double        read_ref_ratio
)
{
   if (matchlist->size < 1) return -1;

   // Convert seeds to hits.
   uint64_t hit_count, rev_count;
   // Forward seeds.
   hit_t * hits = compute_hits(seeds, index, &hit_count);

   // DEBUG.
   if (VERBOSE_DEBUG) {
      for (int i = 0; i < hit_count; i++) {
         uint64_t g_start;
         int dir = 0;
         if (hits[i].locus >= index->size/2) {
            g_start = index->size - (hits[i].locus+hits[i].depth) - 2;
            dir = 1;
         } else {
            g_start = hits[i].locus + 1;
            dir = 0;
         }
         // Search chromosome name.
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "[%d,%d] (%d) %s,%ld%c\n",
                 hits[i].qrypos, hits[i].qrypos+hits[i].depth-1,
                 hits[i].depth,
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 dir ? '-' : '+'
                 );
      }
   }

   // Reverse seeds.
   if (rseeds->pos > 0) {
      hit_t * rhits = compute_hits(rseeds, index, &rev_count);
      // Revert hits.
      for (int i = 0; i < rev_count; i++) {
         rhits[i].locus = (index->size - 1) - (rhits[i].locus + rhits[i].depth);
         rhits[i].qrypos= seqlen - (rhits[i].qrypos + rhits[i].depth);
      }
      // Merge F/R hits.
      hits = realloc(hits, (hit_count + rev_count)*sizeof(hit_t));
      memcpy(hits+hit_count, rhits, rev_count*sizeof(hit_t));
      hit_count += rev_count;
      free(rhits);
   }
   
   // Merge-sort hits.
   //   mergesort_mt(hits, hit_count, sizeof(hit_t), 0, 1, compar_hit_locus);
   qsort(hits, hit_count, sizeof(hit_t), compar_hit_locus_qsort);

   if (VERBOSE_DEBUG) fprintf(stdout, "** CHAIN SEEDS (total %ld)\n", hit_count);

   // Reset matchlist.
   matchlist->pos = 0;

   // Aux variables.
   long minv = 0;
   int  mini = 0;
   size_t i = 0;

   while (i < hit_count) {
      // Match start.
      int64_t span     = hits[i].depth;
      int64_t lbeg_tmp = hits[i].locus;
      int64_t lend_tmp = lbeg_tmp + span;
      int64_t rbeg_tmp = hits[i].qrypos;
      int64_t rend_tmp = rbeg_tmp + span;
      int     bulk     = hits[i].bulk;
      int64_t lbeg, lend, rbeg, rend;
      lbeg = lend = rbeg = rend = 0;
      i++;

      // Debug seeds...
      if (VERBOSE_DEBUG) {
         uint64_t g_start;
         int dir = 0;
         if (lbeg_tmp >= index->size/2) {
            g_start = index->size - lend_tmp - 2;
            dir = 1;
         } else {
            g_start = lbeg_tmp + 1;
            dir = 0;
         }
         // Search chromosome name.
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "%s,%ld%c/%ld,%ld/%ld,%ld,%d",
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 rbeg_tmp, rend_tmp-1,
                 span, span, bulk
                 );
      }

      // Extend match.
      while (i < hit_count) {
         // Compare genome and read distances.
         int64_t g_dist = (int64_t)(hits[i].locus) - (int64_t)(hits[i-1].locus);
         int64_t r_dist = (int64_t)(hits[i].qrypos) - (int64_t)(hits[i-1].qrypos);
         if (r_dist >= 0 && (g_dist <= (1+read_ref_ratio) * r_dist) && (g_dist >= (1-read_ref_ratio) * r_dist)) {
            uint64_t new_lend = hits[i].locus + hits[i].depth;
            uint64_t new_rend = hits[i].qrypos + hits[i].depth;
            if (new_rend > rend_tmp) {
               // Update hits.
               if (hits[i].qrypos <= rend_tmp) {
                  span += new_rend - rend_tmp;
               }
               else {
                  span += hits[i].depth;
                  // Check if this is the biggest stretch.
                  if (rend_tmp - rbeg_tmp > rend - rbeg) {
                     rbeg = rbeg_tmp;
                     rend = rend_tmp;
                     lbeg = lbeg_tmp;
                     lend = lend_tmp;
                  }
                  // Follow new stretch.
                  rbeg_tmp = hits[i].qrypos;
                  lbeg_tmp = hits[i].locus;
               }
               // Update ends.
               rend_tmp = new_rend;
               lend_tmp = new_lend;

            }
            bulk *= hits[i].bulk;

            // Debug seed chain...
            if (VERBOSE_DEBUG) {
               uint64_t g_start;
               int dir = 0;
               if (hits[i].locus >= index->size/2) {
                  g_start = index->size - hits[i].locus - hits[i].depth - 2;
                  dir = 1;
               } else {
                  g_start = hits[i].locus + 1;
                  dir = 0;
               }
               // Search chromosome name.
               int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
               // Print results.
               fprintf(stdout, "\t%s,%ld%c/%d,%d/%d,%ld,%d",
                       index->chr->name[chrnum],
                       g_start - index->chr->start[chrnum]+1,
                       dir ? '-' : '+',
                       hits[i].qrypos, hits[i].qrypos + hits[i].depth-1,
                       hits[i].depth, span, bulk
                       );
            }

            i++;
         } else break;
      }
      // Check if last was the biggest stretch.
      if (rend_tmp - rbeg_tmp > rend - rbeg) {
         rbeg = rbeg_tmp;
         rend = rend_tmp;
         lbeg = lbeg_tmp;
         lend = lend_tmp;
      }


      if (VERBOSE_DEBUG) fprintf(stdout, "\n");

      // Store match.
      // Non-significant streaks (bulk) will not be saved.
      if (span > minv && bulk == 0) {
         // Allocate match.
         match_t match;
         match.ref_e  = lend - 1;
         match.read_e = rend - 1;
         match.ref_s  = lbeg;
         match.read_s = rbeg;
         match.hits   = span;
         match.flags  = 0;
         match.score  = 0;

         // Debug match
         if (VERBOSE_DEBUG) {
            uint64_t g_start;
            int dir = 0;
            if (match.ref_s >= index->size/2) {
               g_start = index->size - match.ref_e - 2;
               dir = 1;
            } else {
               g_start = match.ref_s + 1;
               dir = 0;
            }
            // Search chromosome name.
            int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
            // Print results.
            fprintf(stdout, "[%d]\t%s,%ld%c\t%ld\t%d,%d\n",
                    match.hits,
                    index->chr->name[chrnum],
                    g_start - index->chr->start[chrnum]+1,
                    dir ? '-' : '+',
                    match.ref_s,
                    match.read_s, match.read_e
                    );
         }

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
         }
         else {
            if(VERBOSE_DEBUG) {
               match_t m = matchlist->match[mini];
               fprintf(stdout, "replaces [%d]\t%ld\t%d,%d\n",m.hits,m.ref_s,m.read_s,m.read_e);
            }
            matchlist->match[mini] = match;
         }
               
         // Find the minimum that will be substituted next.
         if (matchlist->pos == matchlist->size) {
            mini = 0;
            minv = matchlist->match[mini].hits;
            for (int j = 1 ; j < matchlist->size; j++) {
               if (matchlist->match[j].hits < minv) {
                  mini  = j;
                  minv = matchlist->match[j].hits;
               }
            }
         }
      }
   }

   // Free hits.
   free(hits);

   return 0;
}

int
compar_hit_locus_qsort
(
 const void * a,
 const void * b
)
{
   hit_t * ha = (hit_t *)a;
   hit_t * hb = (hit_t *)b;
   if (ha->locus > hb->locus) return 1;
   else if (ha->locus < hb->locus) return -1;
   else return (ha->depth > hb->depth ? 1 : -1);
}

int
get_smem
(
 seedstack_t ** smem,
 seedstack_t  * mem
)
{
   if (mem->pos == 0) return 0;
   // Seeds are sorted by start position, descending.
   seed_t cur = mem->seed[mem->pos-1];
   int    len = cur.ref_pos.depth;
   int    end = cur.qry_pos + cur.ref_pos.depth;
   // Find longest non-overlapping seeds.
   for (int i = mem->pos-2; i >= 0; i--) {
      seed_t s = mem->seed[i];
      if (s.qry_pos >= end) {
         seedstack_push(cur, smem);
         cur = s;
         len = cur.ref_pos.depth;
         end = cur.qry_pos + cur.ref_pos.depth;
      } else if (s.ref_pos.depth > len) {
         cur = s;
         len = cur.ref_pos.depth;
         end = cur.qry_pos + cur.ref_pos.depth;
      }
   }
   seedstack_push(cur, smem);

   return 0;
}


seedstack_t *
naive_smem
(
 char * seq,
 seedopt_t opt,
 index_t * index
)
{
   uint32_t slen = strlen(seq);
   int32_t qry_end = slen - 1;
   int32_t last_qry_pos = slen;
   // Create stack.
   seedstack_t * stack = seedstack_new(SEEDSTACK_SIZE);

   for (int i = slen-1; i > 0; i--) {
      // Start position.
      bwpos_t pos = (bwpos_t){.depth = 0, .sp = 0, .ep = index->size-1};
      qry_end = i;
      int j;
      uint64_t new_loci;
      for (j = 0; j <= i; j++) {
         bwpos_t newpos;
         int nt = translate[(uint8_t)seq[i-j]];
         // Extend suffix (Backward search).
         suffix_extend(nt, pos, &newpos, index);
         // Count loci.
         new_loci = (newpos.ep < newpos.sp ? 0 : newpos.ep - newpos.sp + 1);

         if (new_loci < opt.min_loci || nt == 4) {
            // Check previous suffix.
            uint64_t loci = (pos.ep < pos.sp ? 0 : pos.ep - pos.sp + 1);
            int32_t qry_pos = qry_end + 1 - pos.depth;
            if (loci <= opt.max_loci && pos.depth >= opt.min_len && qry_pos < last_qry_pos) {
               // Seed found.
               seed_t seed = (seed_t) {.qry_pos = qry_pos, .ref_pos = pos};
               seedstack_push(seed, &stack);
               last_qry_pos = qry_pos;
            }
            if (nt == 4) {
               i = i - j;
            }
            break;
         } else if (newpos.depth == opt.max_len || j == i) {
            int32_t qry_pos = qry_end + 1 - newpos.depth;
            if (new_loci <= opt.max_loci && qry_pos < last_qry_pos) {
               // Seed found.
               seed_t seed = (seed_t) {.qry_pos = qry_pos, .ref_pos = newpos};
               seedstack_push(seed, &stack);
               last_qry_pos = qry_pos;
            }
            break;
         }
         pos = newpos;
      }
   }
   return stack;
}


seedstack_t *
seedstack_new
(
 int size
)
{
   seedstack_t * stack = malloc(sizeof(seedstack_t) + size*sizeof(seed_t));
   if (stack == NULL) return NULL;
   stack->pos = 0;
   stack->size = size;
   return stack;
}

int
seedstack_push
(
 seed_t         seed,
 seedstack_t ** stackp
)
{
   seedstack_t * stack = *stackp;
   // Realloc if full.
   if (stack->pos >= stack->size) {
      size_t newsize = stack->size * 2;
      stack = realloc(stack, sizeof(seedstack_t) + newsize * sizeof(seed_t));
      if (stack == NULL) return -1;
      *stackp = stack;
      stack->size = newsize;
   }
   // Push seed.
   stack->seed[stack->pos++] = seed;
   return 0;
}

int
seedstack_append
(
 seedstack_t ** dest,
 seedstack_t  * src
)
{
   if (src->pos < 1) return 0;

   seedstack_t * stack = *dest;
   // Realloc if full.
   size_t newpos = stack->pos + src->pos;
   if (newpos >= stack->size) {
      size_t newsize = newpos + 1;
      stack = realloc(stack, sizeof(seedstack_t) + newsize * sizeof(seed_t));
      if (stack == NULL) return -1;
      *dest = stack;
      stack->size = newsize;
   }
   // Copy seeds.
   memcpy(stack->seed + stack->pos, src->seed, src->pos * sizeof(seed_t));
   stack->pos = newpos;
   return 0;
}

int
compar_hit_locus
(
 const void * a,
 const void * b,
 const int param
)
{
   hit_t * ha = (hit_t *)a;
   hit_t * hb = (hit_t *)b;
   if (ha->locus > hb->locus) return 1;
   else if (ha->locus < hb->locus) return -1;
   else return (ha->depth > hb->depth ? 1 : -1);
}

int
compar_seed_start
(
 const void * a,
 const void * b,
 const int param
)
{
   seed_t * sa = (seed_t *)a;
   seed_t * sb = (seed_t *)b;
   if (sa->qry_pos > sb->qry_pos) return 1;
   else if (sa->qry_pos < sb->qry_pos) return -1;
   else return (sa->ref_pos.depth > sb->ref_pos.depth ? -1 : 1);
}
