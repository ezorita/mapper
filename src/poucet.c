#include "poucet.h"

int
poucet
(
 const long   sp,
 const long   ep,
 const char * prow,
 const int    depth,
 char *       path,
 arg_t      * arg
)
{
   int tau = arg->tau;
   // Index.
   uint64_t  * c      = arg->index->c;
   uint64_t ** occs   = arg->index->occ;

   // Penalty for match/mismatch and insertion/deletion resepectively.
   uint mmatch;
   uint shift;
   uint score;

   // Part of the cache that is shared between all the children.
   char r[2*MAXTAU+3];
   for (int i = 0; i < MAXTAU + 2; i++)
      r[i] = r[2*MAXTAU+2-i] = MAXTAU + 1 - i;
   char * row = r + MAXTAU + 1;

   // Upper arm of the L (need the path).
   for (int a = tau ; a > 0 ; a--) {
      mmatch = prow[a] + (path[depth-a] != arg->query[depth]);
      shift = min(prow[a-1], row[a+1]) + 1;
      row[a] = min(mmatch, shift);
   }


   // start at base=1 ('@' is never going to be queried).
   for (int nt = 0 ; nt < NUM_BASES ; nt++) {
      // Check whether child 'i' exists.
      long newsp, newep;
      /*
      long occsp = bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1);
      long occep = bisect_search(0, occs[nt].max-1, occs[nt].val, ep);
      newsp = c[nt] + (occs[nt].max ? occsp : 0);
      newep = c[nt] + (occs[nt].max ? occep : 0) - 1;
      */
      newsp = c[nt] + compute_occ(sp-1, occs[nt]);
      newep = c[nt] + compute_occ(ep  , occs[nt]) - 1;
      // Skip if current node has no child at this position.
      if (newep < newsp) continue;

      // Horizontal arm of the L (need previous characters).
      for (int i = tau ; i > 0 ; i--) {
         mmatch = prow[-i] + (nt != arg->query[depth-i]);
         shift = min(prow[1-i], row[-i-1]) + 1;
         row[-i] = min(mmatch, shift);
      }


      // Center cell (need both arms to be computed).
      mmatch = prow[0] + (nt != arg->query[depth]);
      shift = min(row[-1], row[1]) + 1;
      row[0] = min(mmatch, shift);

      int g_offset = 0;
      score = tau+1;
      for (int i = -tau ; i <= tau; i++) {
         if (row[i] <= score) {
            score = row[i];
            g_offset = max(0, i);
         }
      }

      // Stop searching if 'tau' is exceeded.
      if (score > arg->tau)
         continue;

      // Reached height of the trie: it's a hit!
      if (depth == arg->qlen) {
         // Find number of cummulative in/dels.
         pebble_t hit = {
            .sp = newsp,
            .ep = newep,
            .nwrow[0] = arg->qlen - g_offset - 1
         };
         ppush(arg->hits + score, hit);
         continue;
      }

      // Cache nodes in pebbles when trailing.
      if (depth <= arg->trail) {
         pebble_t pebble = {
            .sp    = newsp,
            .ep    = newep
         };
         memcpy(&(pebble.nwrow), r, 2*MAXTAU+3);
         ppush(arg->pebbles + depth, pebble);
      }

      // Update path.
      path[depth] = nt;

      // Dash path if mismatches exhausted.
      if (depth > arg->trail && score == arg->tau) {
         for (int i = -tau; i <= tau; i++) {
            if (row[i] == score) {
               dash(newsp, newep, depth+1, i, path, arg);
            }
         }
         continue;
      }

      // Recursive call.
      poucet(newsp, newep, row, depth+1, path, arg);
   }

   return 0;
}


void
dash
(
 long    sp,
 long    ep,
 const int     depth,
 const int     align,
 const char  * path,
 const arg_t * arg
)
// TODO:
//  - update.
// SYNOPSIS:                                                              
//   Checks whether the index has the given suffix and reports a hit if this   
//   is the case.                                                         
//                                                                        
// PARAMETERS:                                                            
//   node: the node to test                                               
//   suffix: the suffix to test as a translated sequence                  
//                                                                        
// RETURN:                                                                
//   'void'.                                                              
//                                                                        
// SIDE EFFECTS:                                                          
//   Updates 'arg.hits' if the suffix is found.                           
{
   uint64_t  * c    = arg->index->c;
   uint64_t ** occs = arg->index->occ;
   int i = depth;
   int nt;

   if (align > 0) {
      for (int k = align ; k > 0 ; k--){
         if (arg->query[i] == EOS) break;
         if (path[depth - k] != arg->query[i++]) return;
      }
   } else {
      i += align;
   }
   
   while ((nt = arg->query[i++]) != EOS) {
      sp = c[nt] + compute_occ(sp-1, occs[nt]);
      ep = c[nt] + compute_occ(ep  , occs[nt]) - 1;
      if (ep < sp) return;
   }

   // i has been increased passed EOS.
   i -= 2;

   pebble_t hit = {
      .sp = sp,
      .ep = ep,
      .nwrow[0] = i - 1 - max(align, 0)
   };
   ppush(arg->hits + arg->tau, hit);
}
