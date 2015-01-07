#include "fastseed.h"
#include <time.h>

int
fastseed_filter
(
 seqstack_t * seqs,
 index_t      index,
 chr_t      * chr,
 int          k
)
{
   for (int i = 0; i < seqs->pos; i++) {
      seq_t seq = seqs->seq[i];
      int  slen = strlen(seq.seq);
      vstack_t * loci = new_stack(1024);

      // 1. Find anchors.
      for (int j = 0; j < slen-k; j++) {
         clock_t tstart = clock();
         long l[2];
         char tmp = seq.seq[j+k];
         seq.seq[j+k] = 0;
         if (query_index(seq.seq + j, index.gsize, index.c, l, index.occ))
            pushvec(&loci, index.pos + l[0], l[1] - l[0] + 1);
         seq.seq[j+k] = tmp;
      }
      
      // 2. Sort anchor loci.
      long aux[loci->pos];
      memcpy(aux, loci->val, loci->pos * sizeof(long));
      mergesort_long(loci->val, aux, loci->pos, 0);

      // 3. Filter regions.
      // Grab loci that are close to others and make a single aligment region. Maybe make
      // a stack containing all these regions and sort them by the number of hits that they contain.
      // This can be a stack of align_t structs. (start of the alignment, end and score for the num of hits).
      align_t aligns[loci->pos/2];      
      int j = 0, streak = 0, cnt = 0;
      long DIST = 1000;
      while (j < loci->pos) {
         while (j + streak + 1 < loci->pos && loci->val[j + streak + 1] - loci->val[j + streak] < DIST) streak++;
         if (streak) {
            // Save loci.
            align_t a = {.start = loci->val[j], .max = loci->val[j+streak] + k, .score = streak};
            aligns[cnt++] = a;
            j += streak;
            streak = 0;
         }
         j++;
      }

      // May want to align only the most significant ones.

      // 4. Do local alignments.
      fprintf(stderr, "Alignments that will be computed: %d\n", cnt);
      for (long l = 0; l < cnt; l++) {
         long align_start = aligns[l].start - ALIGNMENT_EXTRA_L;
         long align_end   = aligns[l].max + ALIGNMENT_EXTRA_R;
         align_start = (align_start < 0 ? 0 : align_start);
         align_end   = (align_end >= index.gsize ? index.gsize - 1 : align_end);

         char genomeseq[align_end - align_start + 1];
         for (int m = align_end, n = 0; m >= align_start; m--, n++) genomeseq[n] = index.genome[m];
         align_t align = sw_align(seq.seq + j, slen, genomeseq, align_end - align_start + 1);
         // Identify chromosome.
         long locus = index.gsize - (aligns[l].start + align.start);
         int chrnum = bisect_search(0, chr->nchr-1, chr->start, locus + 1) - 1;
         fprintf(stdout, "\t%s\t%ld\t%ld-%ld\t%s:%ld\n", seq.tag, align.score, align.start, align.max, chr->name[chrnum], locus - chr->start[chrnum]);
      }


   }
   return 0;  
}

int
fastseed
(
 seqstack_t * seqs,
 index_t      index,
 chr_t      * chr,
 int          k
)
{
   for (int i = 0; i < seqs->pos; i++) {
      seq_t seq = seqs->seq[i];
      int  slen = strlen(seq.seq);
      for (int j = 0; j < slen-k; j++) {
         clock_t tstart = clock();
         long l[2];
         char tmp = seq.seq[j+k];
         seq.seq[j+k] = 0;
         if (query_index(seq.seq + j, index.gsize, index.c, l, index.occ)) {
            long align_start = l[0] - j - ALIGNMENT_EXTRA_L;
            long align_end   = l[1] + k - j + ALIGNMENT_EXTRA_R;
            align_start = (align_start < 0 ? 0 : align_start);
            align_end   = (align_end >= index.gsize ? index.gsize - 1 : align_end);

            // 30/12/2014 CHECK CALL TO SW_ALIGN, seems that one call has an excessive SEQ SIZE.
            align_t align = sw_align(seq.seq + j, slen, index.genome + align_start, align_end - align_start);

            // Find genome locus.
            long locus = index.gsize - (l[0] + align.start);

            // Identify chromosome.
            int chrnum = bisect_search(0, chr->nchr-1, chr->start, locus + 1) - 1;
            fprintf(stdout, "\t%s\t%ld\t%ld-%ld\t%s:%ld\n", seq.tag, align.score, align.start, align.max, chr->name[chrnum], locus - chr->start[chrnum]);
         }
         seq.seq[j+k] = tmp;
         fprintf(stderr, "query time: %ld us\n", ((clock()-tstart)*1000000)/CLOCKS_PER_SEC);
      }
   }
   return 0;  
}

int
query_index
(
 char  * query,
 long    gsize,
 long  * c,
 long  * ptr,
 list_t * occs
)
{
   long sp, ep;
   int slen = strlen(query);
   int qval[slen];
   translate_query(query, qval, slen);
   int i = 0;
   // Initialize range.
   sp = c[qval[0]];
   ep = c[qval[0] + 1] - 1;

   // A pebble contains sp, ep for the start node and the current alignment.

   // Iterate over FM-index.
   for (i = 1; i < strlen(query) && sp <= ep; i++) {
      int nt = qval[i];
      long occsp = bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1);
      long occep = bisect_search(0, occs[nt].max-1, occs[nt].val, ep);
      sp = c[nt] + (occs[nt].max ? occsp : 0);
      ep = c[nt] + (occs[nt].max ? occep : 0 ) - 1;
   }
   if (ep < sp) {
      return 0;
   }

   // Return hits.
   ptr[0] = sp;
   ptr[1] = ep;
   
   return ep-sp+1;
}

void
translate_query
(
 char * query,
 int  * qval,
 int    qlen
)
{
   char translate[256] = {[0 ... 255] = 4, ['@'] = 0,
                           ['a'] = 1, ['c'] = 2, ['g'] = 3, ['n'] = 4, ['t'] = 5,
                           ['A'] = 1, ['C'] = 2, ['G'] = 3, ['N'] = 4, ['T'] = 5 };

   // Converts to int and reverses query.
   for (int i = 0; i < qlen; i++) qval[i] = translate[(int)query[qlen - i - 1]];
}
