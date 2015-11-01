#include "format.h"

int
print_and_free
(
 seq_t          seq,
 matchlist_t  * matches,
 index_t      * index,
 formatopt_t    opt
)
{

   // Retrieve match.
   if (matches->pos > 0) {
      match_t match = matches->match[0];
      // Check mapping quality.
      // Format output.
      int dir;
      uint64_t g_start, g_end = 0;
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
      int chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;

      // Print results.
      fprintf(stdout, "%s\t%d\t%d\t%s:%ld-%ld:%c\t\n",
              seq.tag,
              match.read_s+1, match.read_e+1,
              index->chr->name[chrnum],
              g_start - index->chr->start[chrnum]+1,
              g_end - index->chr->start[chrnum]+1,
              dir ? '-' : '+');

   }
   else {
      fprintf(stdout, "%s\t0\t0\t*\n", seq.tag);
   }

   // Free match.
   free(seq.seq);
   free(seq.tag);
   free(seq.q);

   return 0;
}
