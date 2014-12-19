#include "bwmapper.h"
#include <time.h>

int visited = 0;

void SIGSEGV_handler(int sig) {
   void *array[10];
   size_t size;

   // get void*'s for all entries on the stack
   size = backtrace(array, 10);

   // print out all the frames to stderr
   fprintf(stderr, "Error: signal %d:\n", sig);
   backtrace_symbols_fd(array, size, STDERR_FILENO);
   exit(1);
}


int main(int argc, char * argv[]) {

   // TODO: parametrize
   int opt_reverse  = 1;
   int opt_verbose  = 1;
   int tau = 2;

   // Backtrace handler
   signal(SIGSEGV, SIGSEGV_handler); 

   if (argc < 3) {
      fprintf(stderr, "usage: bwmapper {query <file.fastq> | index} <genome file>\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1],"query") == 0) {
      char * filename = malloc(strlen(argv[3])+5);
      strcpy(filename, argv[3]);

      // mmap index file.
      strcpy(filename+strlen(argv[3]), ".fmi");
      int fd = open(filename, O_RDONLY);
      if (fd == -1) {
         fprintf(stderr, "error opening '%s': %s\n", filename, strerror(errno));
         exit(EXIT_FAILURE);
      }

      // Parse query filename.
      FILE * queryfile = fopen(argv[2], "r");
      if (queryfile == NULL) {
         fprintf(stderr, "error: could not open file.\n");
         return EXIT_FAILURE;
      }

      // Read file.
      fprintf(stderr, "reading query file...\n");
      seqstack_t * seqs = read_file(queryfile, opt_reverse, opt_verbose); // TODO: set reverse and verbose options.
      if (seqs == NULL) {
         return EXIT_FAILURE;
      }

      // Sort sequences.
      /*
      fprintf(stderr, "sorting input sequences...\n");
      if(seqsort(&(seqs->seq[0]), seqs->pos, opt_nthreads)) {
         return EXIT_FAILURE;
      }
      */

      int mflags = MAP_PRIVATE;
      /*
      if (seqs->pos < 1000) {
         fprintf(stderr, "reading index (fly mode)...\n");
      } else {
         fprintf(stderr, "pre-loading index (populate)...\n");
         mflags |= MAP_POPULATE;
      }
      */
      long idxsize = lseek(fd, 0, SEEK_END);
      lseek(fd, 0, SEEK_SET);
      long * indexp = mmap(NULL, idxsize, PROT_READ, mflags, fd, 0);
      if (indexp == NULL) {
         fprintf(stderr, "error opening index file: %s.\n", strerror(errno));
         return EXIT_FAILURE;
      }
      
      // Read FM index format.
      index_t index;
      if(format_FMindex(indexp, &index))
         return EXIT_FAILURE;

      // Read chromosome index.
      strcpy(filename+strlen(argv[3]), ".index");
      chr_t * chr = read_CHRindex(filename);
      if (chr == NULL)
         return EXIT_FAILURE;

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
      int seq_block = 100;
      // Allocate hitmaps.
      vstack_t ** hitmaps = malloc(seq_block * sizeof(vstack_t *));
      for (int i = 0 ; i < seq_block ; i++) {
         hitmaps[i] = new_stack(HITMAP_SIZE);
         if (hitmaps[i] == NULL) return -1;
      }
      // Allocate loci stack.
      int max_loci = 25;
      loci_t * loci = malloc(sizeof(loci_t) + max_loci * sizeof(pebble_t));
      if (loci == NULL) return -1;
      loci->size = max_loci;
      loci->pos  = 0;

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

         // DEBUG.
         clock_t istart = clock();

         // Poucet algorithm over the k-mers.
         for (int i = 0; i < subseqs->size; i++) {
            // DEBUG.
            clock_t tstart = clock();
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
            visited = 0;
            for (int p = 0 ; p < pebbles[start]->pos ; p++) {
               // Next pebble.
               pebble_t pebble = pebbles[start]->pebble[p];
               // Compute current alignment from alignment trie.
               int wingsz;
               trie_getrow(trie, pebble.rowid >> SCORE_BITS, pebble.rowid & SCORE_MASK, &wingsz, nwrow);
               // Recover the current path.
               long gpos = index.pos[pebble.sp];
               for (int j = 0; j < start; j++) {
                  path[start-j] = translate[(int)index.genome[gpos+j]];
               }
               poucet(pebble.sp, pebble.ep, wingsz, nwrow, start + 1, path, &arg);
            }
            //            fprintf(stdout, " (visited nodes: %d)\n", visited);

            // Map hits.
            map_hits(hits, query.hitmap, &index, tau, i);

            // Sort hitmap.
            start = trail;
            // DEBUG.
            fprintf(stdout, "query time: %luus\n",((clock()-tstart)*1000000)/CLOCKS_PER_SEC);
         }
         // DEBUG.
         fprintf(stdout, "TOTAL query time: %luus\n",((clock()-istart)*1000000)/CLOCKS_PER_SEC);

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
            hitmap_analysis(hitmaps[i], loci, KMER_SIZE - tau, WINDOW_SIZE);
            // Print results.
            for (int k = 0; k < loci->pos; k++) {
               pebble_t locus = loci->loci[k];
               long start = index.gsize - locus.ep;
               long end   = index.gsize - locus.sp;
               int chrnum = bisect_search(0, chr->nchr-1, chr->start, start+1)-1;
               fprintf(stdout, "%ld\t%s\t%s:%ld-%ld (%ld)\n",
                       locus.rowid,
                       seqs->seq[s+i].tag,
                       chr->name[chrnum],
                       start - chr->start[chrnum]+1,
                       end - chr->start[chrnum]+1,
                       locus.ep - locus.sp + 1);
            }
         }
      }
   }
   else if (strcmp(argv[1],"index") == 0) {
      write_index(argv[2]);
   }
   else 
      fprintf(stderr, "usage: bwmapper {query <seq> | index} <genome file>\n");
   
   return 0;
}


/*********************/
/** query functions **/
/*********************/

int
hitmap_analysis
(
 vstack_t * hitmap,
 loci_t   * loci,
 int        mindist,
 int        maxdist
)
{
   if (loci->size < 1) return -1;
   long minv = 0;
   int  min  = 0;
   
   // Initialize loci list.
   loci->pos = 0;
   
   // Find clusters in hitmap.
   int   streak = 0;
   pebble_t hit;
   for (int i = 0; i < hitmap->pos - 1; i++) {
      long diff = hitmap->val[i+1] - hitmap->val[i];
      if (diff < maxdist) {
         if (diff > mindist) {
            if (!streak) {
               hit.sp = hitmap->val[i];
               streak = 2;
            } else streak++;
         }
      } else {
         if (streak) {
            hit.ep = hitmap->val[i];
            hit.rowid = streak;
            // Check significance.
            if (streak > minv) {
               loci->loci[min] = hit;
               // Extend loci list if not yet full.
               if (loci->pos == min && loci->pos < loci->size) min = ++loci->pos;
               // Find minimum.
               else {
                  min = 0;
                  minv = loci->loci[0].rowid;
                  for (int j = 1 ; j < loci->pos; j++) {
                     if (minv > loci->loci[j].rowid) {
                        min  = j;
                        minv = loci->loci[j].rowid;
                     }
                  }
               }
               // End find minimum.
            }
            streak = 0;
         }
      }
   }
   // Sort mapped regions by significance.
   // mergesort_score(...);
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
   // Push hits to hitmap.
   for (int a = 0 ; a < tau ; a++) {
      for (int h = 0; h < hits[a]->pos; h++) {
         pebble_t hit = hits[a]->pebble[h];
         // Filter out too abundant sequences. (To speed up the sorting).
         if (hit.ep - hit.sp + 1 < HIT_MAX_LOCI) {
            if(pushvec(hitmap, index->pos + hit.sp, hit.ep - hit.sp + 1)) return -1;
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

   seqsort(list->sub, list->size, k, 1);

   return list;
}


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

      visited++;

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

      score = MAXTAU;
      for (int i = -wingsz ; i <= wingsz; i++) {
         if (row[i] < score) score = row[i];
      }

      // DEBUG.
      /*
      for (int k = 1; k <= depth; k++) fprintf(stdout, "%c", revert[(int)arg->query[k]]);
      fprintf(stdout, "\n");
      for (int k = 1 ; k < depth; k++) fprintf(stdout, "%c", revert[(int)path[k]]);
      fprintf(stdout, "%c\t%d\n", revert[nt], score);
      */
      // Stop searching if 'tau' is exceeded.
      if (score > arg->tau)
         continue;

      // Reached height of the trie: it's a hit!
      if (depth == arg->qlen) {
         pebble_t hit = {
            .sp = newsp,
            .ep = newep,
            .rowid = 0
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
            .rowid = (rowid & SCORE_MASK) | (nodeid << SCORE_BITS)
         };
         ppush(arg->pebbles + depth, pebble);
      }

      // Update path.
      path[depth] = nt;

      // Dash path if mismatches exhausted.
      if (depth > arg->trail && score == arg->tau) {
         for (int i = -wingsz; i <= wingsz; i++) {
            if (row[i] == score)
               dash(newsp, newep, depth+1, i, path, arg);
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


   pebble_t hit = {
      .sp = sp,
      .ep = ep,
      .rowid = 0
   };
   ppush(arg->hits + arg->tau, hit);
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
   long sp, ep, vsp=0, vep=0;
   int * qval = translate_query(query);
   int i = 0;
   // Initialize range.
   sp = c[qval[0]];
   ep = c[qval[0] + 1] - 1;

   // A pebble contains sp, ep for the start node and the current alignment.

   // Iterate over FM-index.
   for (i = 1; i < strlen(query) && sp <= ep; i++) {
      vsp = sp;
      vep = ep;
      int nt = qval[i];
      long occsp = bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1);
      long occep = bisect_search(0, occs[nt].max-1, occs[nt].val, ep);
      sp = c[nt] + (occs[nt].max ? occsp : 0);
      ep = c[nt] + (occs[nt].max ? occep : 0 ) - 1;
   }
   if (ep < sp) {
      sp = vsp;
      ep = vep;
      i--;
   }

   // Return hits.
   ptr[0] = sp;
   ptr[1] = ep;
   
   return i;
}


int *
translate_query
(
 char * query
)
{
   // Converts to int and reverses query.
   int   qlen = strlen(query);
   int * qval = malloc(qlen * sizeof(int));
   for (int i = 0; i < qlen; i++) qval[i] = translate[(int)query[qlen - i - 1]];
   return qval;
}


long
bisect_search
(
 long   start,
 long   end,
 long * set,
 long   value
)
{
   if (start >= end - 1) {
      if (value < set[start]) return start;
      return (value >= set[end] ? end : start) + 1;
   }
   long middle = (start + end) / 2;
   if (set[middle] >= value) return bisect_search(start, middle, set, value);
   else return bisect_search(middle, end, set, value);
}



/*********************/
/** index functions **/
/*********************/


chr_t *
read_CHRindex
(
 char * filename
)
{
   // Files
   FILE * input = fopen(filename,"r");
   if (input == NULL) {
      fprintf(stderr, "error in 'read_CHRindex' (fopen): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   int chrcount = 0;
   int structsize = CHRSTR_SIZE;

   long *  start = malloc(structsize*sizeof(long));
   char ** names = malloc(structsize*sizeof(char*));

   // File read vars
   ssize_t rlen;
   size_t  sz = BUFFER_SIZE;
   char  * buffer = malloc(BUFFER_SIZE);
   int lineno = 0;

   // Read chromosome entries.
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      lineno++;
      char *name;
      int i = 0;
      while (buffer[i] != '\t' && buffer[i] != '\n') i++;
      if (buffer[i] == '\n') {
         fprintf(stderr, "illegal format in %s (line %d) - ignoring chromosome: %s\n", filename, lineno, buffer);
         continue;
      }
      
      if (buffer[rlen-1] == '\n') buffer[rlen-1] = 0;

      // Realloc stacks if necessary.
      if (chrcount >= structsize) {
         structsize *= 2;
         start = realloc(start, structsize*sizeof(long));
         names = realloc(names, structsize*sizeof(char*));
         if (start == NULL || names == NULL) {
            fprintf(stderr, "error in 'read_CHRindex' (realloc): %s\n", strerror(errno));
            exit(EXIT_FAILURE);
         }
      }

      // Save chromosome start position.
      buffer[i] = 0;
      start[chrcount] = atol(buffer);
      // Save chromosome name.
      name = buffer + i + 1;
      names[chrcount] = malloc(strlen(name)+1);
      strcpy(names[chrcount], name);
      // Inc.
      chrcount++;
   }

   // Return chromosome index structure.
   chr_t * chrindex = malloc(sizeof(chr_t));
   chrindex->nchr = chrcount;
   chrindex->start= start;
   chrindex->name = names;

   fclose(input);
   return chrindex;
}


int
format_FMindex
(
 long    * index,
 index_t * fmindex
)
{
   // index structure:
   // DATA:
   //   C[NUM_BASES]  (8-byte)
   //   gsize         (8-byte)
   //   genome[gsize] (1-byte)
   //   POS[gsize]    (8-byte)
   //   ooc(1): pos:size:val[0]:val[1]:val[2]:...:val[pos-1]          (8-byte)
   //   ...
   //   occ(NUM_BASES): pos:size:val[0]:val[1]:val[2]:...:val[pos-1]  (8-byte)
   long * start = index;
   long pointer_offset = 0;

   // C values.
   fmindex->c = start;
   pointer_offset += NUM_BASES;

   // Genome size.
   long gsize = fmindex->gsize = start[pointer_offset++];

   // Genome bases. (8-bit data).
   char * genome = (char *)(index + pointer_offset);
   fmindex->genome = genome;
   
   // Go back to (8-byte data).
   pointer_offset = 0;
   start = (long *) (genome + gsize);
   
   // Reverse BWT (pos values).
   fmindex->pos = start + pointer_offset;
   pointer_offset += gsize;

   // Occurrences.
   list_t * occ = fmindex->occ = malloc(NUM_BASES * sizeof(list_t));
   if (occ == NULL)
      return -1;

   for (int i = 0; i < NUM_BASES; i++) {
      occ[i].max = start[pointer_offset++];
      occ[i].val = start + pointer_offset;
      pointer_offset += occ[i].max;
   }

   return 0;
}

ssize_t
write_index
(
 char * filename
)
{
   long gsize;

   // Allocate structures.
   long * C, * pos;
   vstack_t ** occ = malloc(NUM_BASES * sizeof(vstack_t *));

   // Parse genome and convert to integers. (0..NBASES-1)
   fprintf(stdout, "compacting genome...\n");
   char * genome = compact_genome(filename, &gsize);
   C = compute_c(genome, gsize);
   fprintf(stdout, "done\n");

   // Output files.
   char * fmfile = malloc(strlen(filename)+5);
   strcpy(fmfile, filename);
   strcpy(fmfile + strlen(filename), ".fmi");
   int fd = open(fmfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   if (fd == -1) {
      fprintf(stderr, "error in write_index (open %s): %s.\n", fmfile, strerror(errno));
      exit(EXIT_FAILURE);
   }
   free(fmfile);

   // FM index
   bwt_index(genome, gsize, &pos, occ);

   // Write index.
   size_t s = 0;
   while (s < NUM_BASES*sizeof(long)) s += write(fd, C + s/sizeof(long), NUM_BASES*sizeof(long) - s);
   s = 0;
   while (s < sizeof(long)) s += write(fd, &gsize, sizeof(long));
   s = 0;
   while (s < gsize*sizeof(char)) s += write(fd, genome + s/sizeof(char), gsize*sizeof(char) - s);
   s = 0;
   while (s < gsize*sizeof(long)) s += write(fd, pos + s/sizeof(long), gsize*sizeof(long) - s);
   for (int i = 0; i < NUM_BASES; i++) {
      fprintf(stdout, "occ[%d] size = %ld (%ld bytes).\n", i, occ[i]->pos, occ[i]->pos*sizeof(long));
      s = 0;
      while (s < sizeof(long)) s += write(fd, &(occ[i]->pos), sizeof(long));
      fprintf(stdout, "occ[%d]->pos: written %ld bytes.\n", i, s);
      s = 0;
      while (s < occ[i]->pos * sizeof(long)) s += write(fd,((long*) &(occ[i]->val[0])) + s/sizeof(long), occ[i]->pos*sizeof(long) - s);
      fprintf(stdout, "occ[%d]->val: written %ld bytes.\n", i, s);
   }

   return s;
}


char *
compact_genome
(
 char * filename,
 long * genomesize
)
{
   char * fileindex = malloc(strlen(filename)+7);
   strcpy(fileindex, filename);
   strcpy(fileindex + strlen(filename), ".index");

   // Files
   FILE * input = fopen(filename,"r");
   FILE * output = fopen(fileindex, "w");

   if (input == NULL) {
      fprintf(stderr, "error in 'compact_genome' (fopen): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }

   // File read vars
   ssize_t rlen;
   size_t  sz = BUFFER_SIZE;
   char * buffer = malloc(BUFFER_SIZE);

   // Genome storage vars
   char * genome = malloc(GENOME_SIZE);
   long gbufsize = GENOME_SIZE;
   genome[0] = '@';
   *genomesize = 1;
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      if (buffer[0] == '>') {
         buffer[rlen-1] = 0;
         int k = 0;
         while (buffer[k] != ' ' && buffer[k] != 0) k++;
         buffer[k] = 0;
         fprintf(output,"%ld\t%s\n",*genomesize,buffer+1);
      }
      else {
         if (gbufsize < *genomesize + rlen) {
            while (gbufsize < *genomesize + rlen) gbufsize *= 2;
            genome = realloc(genome, gbufsize);
            if (genome == NULL) {
               fprintf(stderr, "error in 'compact_genome' (realloc): %s\n", strerror(errno));
               exit(EXIT_FAILURE);
            }
         }
         int newline = buffer[rlen-1] == '\n';
         strncpy(genome + *genomesize, buffer, rlen-newline);
         *genomesize += rlen-newline;
      }
   }

   // Realloc buffer.
   genome = realloc(genome, *genomesize+1);

   // Reverse genome.
   long end = *genomesize - 1;
   char tmp;
   for (int i = 0 ; i < *genomesize/2; i++) {
      tmp = genome[i];
      genome[i] = genome[end-i];
      genome[end-i] = tmp;
   }
   genome[*genomesize] = 0;

   // Free memory.
   free(buffer);
   free(fileindex);
   
   return genome;
}


void
bwt_index
(
 char      * genome,
 long        gsize,
 long     ** pos,
 vstack_t ** occ
)
{
   // Sort prefixes (will reuse stacks[i] to compute C)
   fprintf(stdout,"suffix sort...\n");
   long * values = malloc((gsize+3)*sizeof(long));
   for (long i = 0; i < gsize; i++) values[i] = (long)translate[(int)genome[i]];
   values[gsize] = values[gsize+1] = values[gsize+2] = 0;
   long * sa = malloc(gsize*sizeof(long));
   suffixArray(values, sa, gsize, NUM_BASES-1);
   free(values);
   //   long * sa = dc3(genome);
   fprintf(stdout,"done.\n");

   // Fill compacted occ.
   vstack_t * stacks[NUM_BASES];
   for (int i = 0; i < NUM_BASES; i++) stacks[i] = new_stack(gsize/NUM_BASES);

   for (long i = 0; i < gsize; i++) push(stacks+translate[(int)genome[(sa[i] == 0 ? gsize-1 : sa[i] - 1)]], i);
   
   // Save stacks to output.
   *pos = sa;
   
   fprintf(stdout, "saving occ structures...\n");
   for (int i = 0; i < NUM_BASES; i++) {
      stacks[i] = realloc(stacks[i], sizeof(vstack_t) + stacks[i]->pos*sizeof(long));
      if (stacks[i] == NULL) {
         fprintf(stderr, "error in bwt_index (realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      stacks[i]->size = stacks[i]->pos;
      occ[i] = stacks[i];
   }
   fprintf(stdout, "done\n");
}


long *
compute_c
(
 char * genome,
 long gsize
)
{
   long * cnt = calloc(NUM_BASES, sizeof(long));
   for (long i = 0; i < gsize; i++) cnt[(int)translate[(int)genome[i]]]++;
   long * c = malloc(NUM_BASES * sizeof(long));
   c[0] = 0;
   for (int i = 1; i < NUM_BASES; i++) c[i] = c[i-1] + cnt[i-1];
   free(cnt);
   return c;
}

/*********************/
/** seq_t functions **/
/*********************/

int
seq_push
(
 seqstack_t ** stackp,
 const char  * tag,
 const char  * seq,
 const int     reverse
)
{
   seqstack_t * stack = *stackp;

   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(seqstack_t) + newsize*sizeof(seq_t));
      if (stack == NULL) {
         return -1;
      }
      stack->size = newsize;
   }

   seq_t * seqt = &(stack->seq[stack->pos++]);
   
   // Copy tag
   seqt->tag = strdup(tag);
   seqt->seq = strdup(seq);

   // Copy sequence (or reverse-complement)
   if (reverse) {
      int len = strlen(seq);
      char * rseq = malloc(len+1);
      for (int i = 0; i <= len; i++)
         rseq[len-1-i] = rcode[(int)seq[i]];
      rseq[len] = 0;
      seqt->rseq = rseq;
   } else {
      seqt->rseq = NULL;
   }
   return 0;
}

seqstack_t *
new_seqstack
(
 int size
 )
{
   if (size < 1) size = 1;
   seqstack_t * stack = malloc(sizeof(seqstack_t) + size * sizeof(seq_t));
   if (stack == NULL)
      return NULL;
   
   stack->pos = 0;
   stack->size = size;
   
   return stack;
}



/*********************/
/** stack functions **/
/*********************/


vstack_t *
new_stack
(
 long size
)
{
   if (size < 2) size = 2;
   vstack_t * stack = malloc(sizeof(vstack_t) + size * sizeof(long));
   if (stack == NULL) return NULL;
   stack->pos  = 0;
   stack->size = size;
   
   return stack;
}


int
push
(
 vstack_t ** stackp,
 long      value
)
{
   vstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(vstack_t) + newsize * sizeof(long));
      if (stack == NULL) {
         fprintf(stderr, "error in 'push' (realloc): %s\n", strerror(errno));
         return -1;
      }
      stack->size = newsize;
   }

   stack->val[stack->pos++] = value;
   return 0;
}

int
pushvec
(
 vstack_t ** stackp,
 long      * vector,
 int         vecsize
)
{
   vstack_t * stack = *stackp;
   if (stack->pos + vecsize >= stack->size) {
      long newsize = stack->size * 2;
      while(newsize <= stack->pos + vecsize) newsize *= 2;
      *stackp = stack = realloc(stack, sizeof(vstack_t) + newsize * sizeof(long));
      if (stack == NULL) {
         fprintf(stderr, "error in 'push' (realloc): %s\n", strerror(errno));
         return -1;
      }
      stack->size = newsize;
   }
   memcpy(stack->val + stack->pos, vector, vecsize * sizeof(long));
   stack->pos += vecsize;
   return 0;
}



pstack_t *
new_pstack
(
 long size
 )
{
   if (size < 2) size = 2;
   pstack_t * stack = malloc(sizeof(pstack_t) + size * sizeof(pebble_t));
   if (stack == NULL) return NULL;
   stack->pos  = 0;
   stack->size = size;

   return stack;
}

void
ppush
(
 pstack_t ** stackp,
 pebble_t pebble
 )
{
   pstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(pstack_t) + newsize * sizeof(pebble_t));
      if (stack == NULL) {
         fprintf(stderr, "error in 'ppush' (realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      stack->size = newsize;
   }

   stack->pebble[stack->pos++] = pebble;
}

/*********************/
/** trie  functions **/
/*********************/

trie_t *
trie_new
(
 int initial_size
)
// TODO: UPDATE.
// SYNOPSIS:                                                              
//   Creates and initializes a new NW-row trie preallocated with the specified
//   number of nodes.
//                                                                        
// PARAMETERS:                                                            
//   initial_size : the number of preallocated nodes.
//   height       : height of the trie. It must be equal to the number of keys
//                  as returned by parse.
//                                                                        
// RETURN:                                                                
//   On success, the function returns a pointer to the new trie_t structure.
//   A NULL pointer is returned in case of error.
//
// SIDE EFFECTS:
//   The returned trie_t struct is allocated using malloc and must be manually freed.
{

   // Allocate at least one node.
   if (initial_size < 1) initial_size = 1;

   trie_t * trie = malloc(sizeof(trie_t) + initial_size*sizeof(node_t));
   if (trie == NULL) {
      fprintf(stderr, "error in 'trie_new' (malloc) trie_t: %s\n", strerror(errno));
      return NULL;
   }

   // Initialize root node.
   memset(&(trie->nodes[0]), 0, initial_size*sizeof(node_t));

   // Initialize trie struct.
   trie->pos = 1;
   trie->size = initial_size;

   return trie;
}


int
trie_getrow
(
 trie_t * trie,
 uint     nodeid,
 int      refval,
 int    * wingsz,
 uint   * nwrow
)
// TODO: UPDATE.
// SYNOPSIS:                                                              
//   Recomputes the NW row that terminates at nodeid.
//                                                                        
// PARAMETERS:                                                            
//   trie   : Pointer to the trie.
//   nodeid : Id of the leaf at which the NW row terminates.
//   refval : The initial condition of the alignment. The score of
//            the rightmost value (top of the inverted L).
//                                                                        
// RETURN:                                                                
//   trie_getrow returns a pointer to the start of the NW row. If an
//   error occurred during the row computation or nodeid did not point
//   to a leaf node, a NULL pointer is returned.
//
// SIDE EFFECTS:
//   An array containing the NW row is allocated using malloc and must be
//   manually freed.
{
   // Return if path starts at root node.
   if (nodeid == 0) {
      *wingsz = 0;
      *nwrow  = refval;
      return 0;
   }

   node_t * nodes  = &(trie->nodes[0]);
   uint     parent = nodeid;
   int      height = 1;   
   while ((parent = nodes[parent].parent) != 0) height++;
   *wingsz = (height-1)/2;

   // Match value.
   int i = *wingsz;
   nwrow[i] = refval;
   uint id = nodeid;
   while (id != 0 && i > -(*wingsz)) {
      uint next_id = nodes[id].parent;
      nwrow[i-1] = nwrow[i] + (nodes[next_id].child[0] == id) - (nodes[next_id].child[2] == id);
      id = next_id;
      i--;
   } 

   // Control.
   if (i != -(*wingsz)) {
      fprintf(stderr, "error in 'trie_getrow': final node != root.\n");
      return -1;
   }

   return 0;
}


uint
trie_insert
(
 trie_t ** triep,
 char   *  path,
 int       pathlen
)
// TODO: UPDATE.
// SYNOPSIS:                                                              
//   Inserts the specified path in the trie and stores the end value and
//   the dfa state at the leaf (last node of the path). If the path already
//   exists, its leaf values will be overwritten.
//                                                                        
// PARAMETERS:                                                            
//   trie     : pointer to a memory space containing the address of the trie.
//   path     : The path as an array of chars containing values {0,1,2}
//   pathlen  : length of the path.
//                                                                        
// RETURN:                                                                
//   On success, dfa_insert returns the id of the leaf where the values were
//   stored, -1 is returned if an error occurred.
//
// SIDE EFFECTS:
//   If the trie has reached its limit of allocated nodes, it will be reallocated
//   doubling its size. The address of the trie may have changed after calling dfa_insert.
{
   trie_t * trie  = *triep;
   node_t * nodes = &(trie->nodes[0]);
   uint id = 0;
   uint initial_pos = trie->pos;

   int i;
   for (i = 0; i < pathlen; i++) {
      if (path[i] < 0 || path[i] >= TRIE_CHILDREN) {
         // Bad path, revert trie and return.
         trie->pos = initial_pos;
         return -1;
      }
      // Walk the tree.
      if (nodes[id].child[(int)path[i]] != 0) {
         id = nodes[id].child[(int)path[i]];
         continue;
      }

      // Create new node.
      if (trie->pos >= trie->size) {
         size_t newsize = trie->size * 2;
         *triep = trie = realloc(trie, sizeof(trie_t) + newsize * sizeof(node_t));
         if (trie == NULL) {
            fprintf(stderr, "error in 'trie_insert' (realloc) trie_t: %s\n", strerror(errno));
            return -1;
         }
         // Update pointers.
         nodes = &(trie->nodes[0]);
         // Initialize new nodes.
         trie->size = newsize;
      }
      
      // Consume one node of the trie.
      uint newid = trie->pos;
      nodes[newid].parent = id;
      nodes[newid].child[0] = nodes[newid].child[1] = nodes[newid].child[2] = 0;
      nodes[id].child[(int)path[i]] = newid;
      trie->pos++;

      // Go one level deeper.
      id = newid;
   }

   return id;
}


void
trie_reset
(
 trie_t * trie
)
// SYNOPSIS:                                                              
//   Resets the trie by pruning the root node. The size of the trie, in terms
//   of preallocated nodes is maintained.
//                                                                        
// PARAMETERS:                                                            
//   trie   : Pointer to the trie.
//                                                                        
// RETURN:                                                                
//   void.
//
// SIDE EFFECTS:
//   None.
{
   trie->pos = 0;
   memset(&(trie->nodes[0]), 0, sizeof(node_t));
}


/*********************/
/** misc  functions **/
/*********************/


seqstack_t *
read_file
(
   FILE      * inputf,
   const int   reverse,
   const int   verbose
)
{
   // TODO: Set a maximum read size in MB, when this limit is reached return.
   // Once the read sequences are processed, read again and process, continuing
   // at the point of the file where the last 'read_file' returned.

   // Read first line of the file to guess format.
   // Store in global variable FORMAT.
   char c = fgetc(inputf);

   if (c == EOF) return NULL;
   if (ungetc(c, inputf) == EOF) {
      return NULL;
   }

   format_t format;

   switch(c) {
   case '>':
      if (verbose) fprintf(stderr, "FASTA format detected\n");
      format = FASTA;
      break;
   case '@':
      if (verbose) fprintf(stderr, "FASTQ format detected\n");
      format = FASTQ;
      break;
   default:
      if (verbose) fprintf(stderr, "raw format detected\n");
      format = RAW;
   }

   ssize_t nread;
   size_t linesz = BUFFER_SIZE;
   size_t sublinesz = BUFFER_SIZE;
   size_t tempsz = BUFFER_SIZE;
   char *line = malloc(BUFFER_SIZE * sizeof(char));
   char *subline = malloc(BUFFER_SIZE * sizeof(char));
   char *temp = malloc(BUFFER_SIZE * sizeof(char));
   if (line == NULL) {
      return NULL;
   }
   seqstack_t *seqstack = new_seqstack(SEQSTACK_SIZE);
   if (seqstack == NULL) {
      return NULL;
   }

   char *seq  = NULL;
   char *tag  = NULL;
   int lineno = 0;

   while ((nread = getline(&line, &linesz, inputf)) != -1) {
      lineno++;
      if (line[nread-1] == '\n') line[nread-1] = '\0';
      switch(format) {
      case RAW:
         seq = line;
         tag = line;
         break;
      case FASTA:
         if (line[0] == '>') {
            if((nread = getline(&subline, &sublinesz, inputf)) == -1) {
               fprintf(stderr, "incorrect FASTA format: line %d\n", lineno);
               continue;
            }
            lineno++;
            if (subline[nread-1] == '\n') subline[nread-1] = '\0';
            tag = line;
            seq = subline;
         } else continue;
         break;
      case FASTQ:
         if (line[0] == '@') {
            if((nread = getline(&subline, &sublinesz, inputf)) == -1) {
               fprintf(stderr, "incorrect FASTQ format: line %d\n", lineno);
               continue;
            }
            lineno++;
            if (subline[nread-1] == '\n') subline[nread-1] = '\0';
            tag = line;
            seq = subline;
            for (int i = 0; i < 2; i++) {
               if((nread = getline(&temp, &tempsz, inputf)) == -1) {
                  fprintf(stderr, "incorrect FASTQ format: line %d\n", lineno);
                  continue;
               }
               lineno++;
            }
         } else continue;
         break;
      }

      size_t seqlen = strlen(seq);
      if (seqlen > MAXSEQLEN) {
         fprintf(stderr, "max sequence length exceeded (%d)\n", MAXSEQLEN);
         fprintf(stderr, "offending sequence:\n%s\n", seq);
         continue;
      }
      
      
      if (seq_push(&seqstack, tag, seq, reverse)) return NULL;
   }

   free(line);
   free(subline);
   free(temp);
   return seqstack;
}


int
seqsort
(
 sub_t * data,
 int     numels,
 int     slen,
 int     thrmax
)
// SYNOPSIS:                                                              
//   Recursive multithreaded merge sort for 'seq_t' arrays.
//   The sequences are sorted alphabetically.
//
// PARAMETERS:                                                            
//   data:   an array of seq_t.                    
//   numels: number of elements in data.                 
//   thrmax: number of threads.                                       
//                                                                        
// RETURN:                                                                
//   
//                                                                        
// SIDE EFFECTS:                                                          
//   Pointers to repeated elements are set to NULL.
{
   // Copy to buffer.
   sub_t *buffer = malloc(numels * sizeof(sub_t));
   memcpy(buffer, data, numels * sizeof(sub_t));

   // Prepare args struct.
   sortargs_t args;
   args.buf0   = data;
   args.buf1   = buffer;
   args.size   = numels;
   // There are two alternating buffers for the merge step.
   // 'args.b' alternates on every call to 'nukesort()' to
   // keep track of which is the source and which is the
   // destination. It has to be initialized to 0 so that
   // sorted elements end in 'data' and not in 'buffer'.
   args.b      = 0;
   args.thread = 0;
   args.slen   = slen;

   // Allocate a number of threads that is a power of 2.
   while ((thrmax >> (args.thread + 1)) > 0) args.thread++;

   nukesort(&args);

   free(buffer);

   return 0;
}

void *
nukesort
(
 void * args
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
   sortargs_t * sortargs = (sortargs_t *) args;
   if (sortargs->size < 2) return 0;

   // Next level params.
   sortargs_t arg1 = *sortargs, arg2 = *sortargs;
   arg1.size /= 2;
   arg2.size = arg1.size + arg2.size % 2;
   arg2.buf0 += arg1.size;
   arg2.buf1 += arg1.size;
   arg1.b = arg2.b = (arg1.b + 1) % 2;

   // Either run threads or DIY.
   if (arg1.thread) {
      // Decrease one level.
      arg1.thread = arg2.thread = arg1.thread - 1;
      // Create threads.
      pthread_t thread1, thread2;
      if ( pthread_create(&thread1, NULL, nukesort, &arg1) ||
           pthread_create(&thread2, NULL, nukesort, &arg2) ) {
         return NULL;
      }
      // Wait for threads.
      pthread_join(thread1, NULL); // TODO: Catch retval and return if retval == -1.
      pthread_join(thread2, NULL);
   }
   else {
      if (nukesort(&arg1)) return NULL;
      if (nukesort(&arg2)) return NULL;
   }

   // Separate data and buffer (b specifies which is buffer).
   sub_t * l = (sortargs->b ? arg1.buf0 : arg1.buf1);
   sub_t * r = (sortargs->b ? arg2.buf0 : arg2.buf1);
   sub_t * buf = (sortargs->b ? arg1.buf1 : arg1.buf0);

   int i = 0;
   int j = 0;
   int idx = 0;
   int cmp = 0;

   // Merge sets
   while (idx < sortargs->size) {
      // Right buffer is exhausted. Copy left buffer...
      if (j == arg2.size) {
         memcpy(buf+idx, l+i, (arg1.size-i) * sizeof(sub_t));
         break;
      }
      // ... or vice versa.
      if (i == arg1.size) {
         memcpy(buf+idx, r+j, (arg2.size-j) * sizeof(sub_t));
         break;
      }
      // Do the comparison.
      sub_t ul = l[i];
      sub_t ur = r[j];
      cmp = strncmp(ul.seq, ur.seq, sortargs->slen);
      if (cmp < 0) buf[idx++] = l[i++];
      else         buf[idx++] = r[j++];
   }

   return NULL;
}

void
mergesort_long
(
 long * data,
 long * aux,
 int    size,
 int    b
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
   long * data2 = data + newsize;
   long * aux2 = aux + newsize;
   mergesort_long(data, aux, newsize, (b+1)%2);
   mergesort_long(data2, aux2, newsize2, (b+1)%2);


   // Separate data and buffer (b specifies which is buffer).
   long * l = (b ? data : aux);
   long * r = (b ? data2 : aux2);
   long * buf = (b ? aux : data);

   int i = 0;
   int j = 0;
   int idx = 0;

   // Merge sets
   while (idx < size) {
      // Right buffer is exhausted. Copy left buffer...
      if (j == newsize2) {
         memcpy(buf+idx, l+i, (newsize-i) * sizeof(long));
         break;
      }
      // ... or vice versa.
      if (i == newsize) {
         memcpy(buf+idx, r+j, (newsize2-j) * sizeof(long));
         break;
      }
      // Do the comparison.
      if (l[i] < r[j]) buf[idx++] = l[i++];
      else             buf[idx++] = r[j++];
   }
}

void
radix_sort
(
 long * a,      // Indices to sort. (may be modified)
 long * b,      // Aux buffer.
 long   n,      // Length of a.
 long   maxval  // Maximum value in a.
)
{
   const int rs_bits = 16;
   const int rs_size = 1 << rs_bits;
   const int rs_mask = rs_size - 1;
 
   int ref = 0, new = 1, it = 0;
   long cnt[rs_size];
   long prf[rs_size];
   long * s[2];

   // Count iterations per value.
   while ((maxval >> (rs_bits*it)) & rs_mask) it++;

   s[0] = a;
   s[1] = b;
   for (long j = 0; j < it; j++) {
      int shift = rs_bits * j;
      // Reset count and prefix.
      memset(cnt, 0, rs_size*sizeof(long));
      prf[0] = 0;
      // Count radix RS_BITS.
      for (long i = 0; i < n; i++) cnt[(s[ref][i] >> shift) & rs_mask]++;
      // Prefix.
      for (int  i = 1; i < rs_size; i++) prf[i] = prf[i-1] + cnt[i-1];
      // Sorted.
      for (long i = 0; i < n; i++) s[new][prf[(s[ref][i] >> shift) & rs_mask]++] = s[ref][i];
      // Swap buffers.
      ref = (ref+1)%2;
      new = (new+1)%2;
   }
   
   // Move data to a.
   if (ref == 1) memcpy(s[0], s[1], n*sizeof(long));
}
