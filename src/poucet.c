#include "poucet.h"

int
poucet
(
 const long   sp,
 const long   ep,
 const int    wingsz,
 const uint * prow,
 const int    depth,
 char *       path,
 arg_t      * arg
)
{
   // Index.
   long   * c      = arg->index->c;
   list_t * occs   = arg->index->occ;

   // Penalty for match/mismatch and insertion/deletion resepectively.
   uint mmatch;
   uint shift;
   uint score;

   // Part of the cache that is shared between all the children.
   uint r[MAXTAU*2+1];
   uint * row = r + MAXTAU;
   // Initialize wing ends.
   row[wingsz+1]  = prow[wingsz]  + 1;
   row[-wingsz-1] = prow[-wingsz] + 1;

   // Upper arm of the L (need the path).
   if (wingsz > 0) {
      for (int a = wingsz ; a > 0 ; a--) {
         mmatch = prow[a] + (path[depth-a] != arg->query[depth]);
         shift = min(prow[a-1], row[a+1]) + 1;
         row[a] = min(mmatch, shift);
      }
   }

   // start at base=1 ('@' is never going to be queried).
   for (int nt = 1 ; nt < NUM_BASES ; nt++) {
      // Check whether child 'i' exists.
      long newsp, newep;
      long occsp = bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1);
      long occep = bisect_search(0, occs[nt].max-1, occs[nt].val, ep);
      newsp = c[nt] + (occs[nt].max ? occsp : 0);
      newep = c[nt] + (occs[nt].max ? occep : 0) - 1;
      // Skip if current node has no child at this position.
      if (newep < newsp) continue;

      // Horizontal arm of the L (need previous characters).
      if (wingsz > 0) {
         for (int i = wingsz ; i > 0 ; i--) {
            mmatch = prow[-i] + (nt != arg->query[depth-i]);
            shift = min(prow[1-i], row[-i-1]) + 1;
            row[-i] = min(mmatch, shift);
         }
      }

      // Center cell (need both arms to be computed).
      mmatch = prow[0] + (nt != arg->query[depth]);
      shift = min(row[-1], row[1]) + 1;
      row[0] = min(mmatch, shift);

      int g_offset = 0;
      score = MAXTAU;
      for (int i = -wingsz ; i <= wingsz; i++) {
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
            .rowid = arg->qlen - g_offset - 1
         };
         ppush(arg->hits + score, hit);
         continue;
      }

      // Extend wings.
      int wsz = wingsz;
      if (score > wingsz) wsz++;
         
      // Cache nodes in pebbles when trailing.
      if (depth <= arg->trail) {
         // Compute differential path.
         char nwpath[2*wsz];
         long rowid = row[wsz];
         for (int i = 2*wsz; i > 0 ; i--)
            nwpath[i-1] = 1 + (row[i] > row[i-1]) - (row[i] < row[i-1]);

         // Insert path into path trie.
         unsigned long nodeid = trie_insert(arg->triep, nwpath, 2*wsz);
         pebble_t pebble = {
            .sp    = newsp,
            .ep    = newep,
            .rowid = (rowid & PATH_SCORE_MASK) | (nodeid << PATH_SCORE_BITS)
         };
         ppush(arg->pebbles + depth, pebble);
      }

      // Update path.
      path[depth] = nt;

      // Dash path if mismatches exhausted.
      if (depth > arg->trail && score == arg->tau) {
         for (int i = -wingsz; i <= wingsz; i++) {
            if (row[i] == score) {
               dash(newsp, newep, depth+1, i, path, arg);
            }
         }
         continue;
      }

      // Recursive call.
      poucet(newsp, newep, wsz, row, depth+1, path, arg);
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
   long   * c      = arg->index->c;
   list_t * occs   = arg->index->occ;
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
      long occsp = bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1);
      long occep = bisect_search(0, occs[nt].max-1, occs[nt].val, ep);
      sp = c[nt] + (occs[nt].max ? occsp : 0);
      ep = c[nt] + (occs[nt].max ? occep : 0 ) - 1;
      if (ep < sp) return;
   }

   // i has been increased passed EOS.
   i -= 2;

   pebble_t hit = {
      .sp = sp,
      .ep = ep,
      .rowid = i - 1 - max(align, 0)
   };
   ppush(arg->hits + arg->tau, hit);
}
