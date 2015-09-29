#include "format.h"

int
print_and_free
(
 seq_t          seq,
 matchlist_t ** intervals,
 int            n_ints,
 index_t      * index,
 formatopt_t    opt
)
{
   int mapped = 0;
   for (int32_t k = 0; k < n_ints; k++) {
      int itv = (k+1)%n_ints;
      matchlist_t * matches = intervals[itv];
      for (size_t j = 0; j < matches->pos; j++) {
         // Retrieve match.
         match_t match = matches->match[j];
         // Check mapping quality.
         if (match.mapq < opt.mapq_thr) {
            if (itv == 0 || !opt.print_first || j > 0) continue;
         }
         else mapped = 1;
         // Format output.
         int dir;
         uint64_t g_start, g_end;
         if (match.ref_s >= index->size/2) {
            g_start = index->size - match.ref_e - 2;
            g_end   = index->size - match.ref_s - 2;
            dir = 1;
         } else {
            g_start = match.ref_s + 1;
            g_end   = match.ref_e + 1;
            dir = 0;
         }
         // Search chromosome name.
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "%s \t%d\t%d\t%d\t%s:%ld-%ld:%c\t%d\t%.2f\t%.1f%%\t%c%c\n",
                 seq.tag, itv,
                 match.read_s+1, match.read_e+1,
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 g_end - index->chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 (int)match.mapq,
                 match.e_exp,
                 match.ident*100.0,
                 (match.flags & WARNING_OVERLAP ? 'o' : '-'),
                 (match.flags & FLAG_FUSED ? 'f' : '-'));
      }
      free(matches);
   }
   return mapped;
}
