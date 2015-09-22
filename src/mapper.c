#include "mapper.h"

#include <time.h>

int main(int argc, char *argv[])
{
   int opt_verbose = 1;

   if (argc != 3) {
      fprintf(stderr, "usage: bwmapper <query file> <genome index>\n");
      exit(EXIT_FAILURE);
   }

   // Parse query filename.
   FILE * queryfile = fopen(argv[1], "r");
   if (queryfile == NULL) {
      fprintf(stderr, "error: could not open file: %s (%s).\n", argv[1],strerror(errno));
      return EXIT_FAILURE;
   }

   // Read FM index format.
   if (opt_verbose) fprintf(stderr, "loading index...\n");
   index_t * index = index_format(index_open(argv[2]));
   if(index == NULL) return EXIT_FAILURE;

   // Read file.
   if (opt_verbose) fprintf(stderr, "reading query file...\n");
   seqstack_t * seqs = read_file(queryfile); // TODO: set reverse and verbose options.
   if (seqs == NULL) {
      return EXIT_FAILURE;
   }

   // Options.
   seedopt_t seedopt = {.min_len = 1, .max_len = 1000, .min_loci = 1, .max_loci = 1, .aux_loci = 200}; // SMEM
   alignopt_t alignopt = {
      .bp_thr       = 10,
      .bp_max_thr   = 50,
      .bp_resolution = 1,
      .bp_period    = 5,
      .bp_repeats   = 4,
      .read_error   = 0.05,
      .rand_error   = 0.50,
      .width_ratio  = 0.05
   };
   mapopt_t mapopt = {
      .dist_accept = 10,
      .read_ref_ratio = 0.05,
      .align_accept_eexp = -30.0,
      .overlap_max_tolerance = 0.5,
      .align_seed_filter_thr = 0.5,
      .align_filter_ident = 0.9,
      .align_filter_eexp = -10.0,
      .align = alignopt,
      .seed = seedopt
   };

   int max_align_per_read = 100;
   double read_ref_ratio = 0.05;
   int dist_accept = 10;
   // Data structures.
   matchlist_t * seed_matches = matchlist_new(max_align_per_read);
   matchlist_t * map_matches = matchlist_new(64);
   // Process reads.
   for (size_t i = 0; i < seqs->pos; i++) {
      //      clock_t t = clock();
      //      fprintf(stderr, "%ld/%ld\r",i, seqs->pos);
      // Seed.
      seedstack_t * seeds = seed(seqs->seq[i].seq, mapopt.seed, index);
      // Match seeds.
      int slen = strlen(seqs->seq[i].seq);
      seed_matches->pos = 0;
      match_seeds(seeds, slen, seed_matches, index, dist_accept, read_ref_ratio);
      // Smith-Waterman alignment.
      map_matches->pos = 0;
      align_seeds(seqs->seq[i].seq, seed_matches, &map_matches, index, mapopt);
      // Sort matches.
      //      mergesort_mt(map_matches->match, map_matches->pos, sizeof(match_t *), 0, 1, compar_matcheexp);
      // Print matches.

      for (size_t j = 0; j < map_matches->pos; j++) {
         match_t * match = map_matches->match[j];
         int dir;
         uint64_t g_start, g_end;
         if (match->ref_s >= index->size/2) {
            g_start = index->size - match->ref_e - 2;
            g_end   = index->size - match->ref_s - 2;
            dir = 1;
         } else {
            g_start = match->ref_s + 1;
            g_end   = match->ref_e + 1;
            dir = 0;
         }
         int   chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;
         // Print results.
         fprintf(stdout, "%s \t%d\t%d\t%s:%ld-%ld:%c\t%.2f\t%.1f%%\t%c%c\n",
                 seqs->seq[i].tag,
                 match->read_s+1, match->read_e+1,
                 index->chr->name[chrnum],
                 g_start - index->chr->start[chrnum]+1,
                 g_end - index->chr->start[chrnum]+1,
                 dir ? '-' : '+',
                 match->e_exp,
                 match->ident*100.0,
                 (match->flags & WARNING_OVERLAP ? 'o' : '-'),
                 (match->flags & FLAG_FUSED ? 'f' : '-'));
      }
      
      //fprintf(stdout, "%.4f\t%ld\t%ld\t%ld\n", (clock()-t)*1000000.0/CLOCKS_PER_SEC,seeds->pos, seed_matches->pos, map_matches->pos);
      free(seeds);
   }

   return 0;
}

double
e_value
(
 int L,
 int m,
 long gsize
)
{
   double E = log10(gsize) + log10(2) * L * (m*3.0/L-2);
   for (int i = 0; i < m; i++) E += log10((L-i)/(double)(m-i));
   return E;
}

int
align_seeds
(
 char         * read,
 matchlist_t  * seeds,
 matchlist_t ** seqmatches,
 index_t      * index,
 mapopt_t       hmargs
 )
{
   if (seeds->pos == 0) return 0;
   int significant = 0;
   int slen = strlen(read);
   matchlist_t * matches = *seqmatches;

   // Sort by seeded nucleotides.
   mergesort_mt(seeds->match, seeds->pos, sizeof(match_t *), 0, 1, compar_seedhits);

   for(long k = 0; k < seeds->pos; k++) {
      // Get next seed.
      match_t * seed = seeds->match[k];

      // Compute alignment limits.
      long r_min  = seed->read_e - seed->read_s + 1;
      long l_min  = 0;
      long r_qstart = seed->read_s + 1;
      long l_qstart = align_max(r_qstart - 1,0);
      long r_rstart = seed->ref_s + 1;
      long l_rstart = seed->ref_s;
     
      // Extend existing mapping.
      double last_eexp = INFINITY;
      match_t * extend = NULL;
      int extend_score = 0;
      int cancel_align = 0;
      for (int i = 0; i < matches->pos; i++) {
         match_t * match = matches->match[i];
         if (match->e_exp > last_eexp) continue;
         // Check whether the seed is contiguous.
         long r_distr = seed->read_e - match->read_e;
         long g_distr = match->ref_e - seed->ref_e;
         long r_distl = match->read_s - seed->read_s;
         long g_distl = seed->ref_s - match->ref_s;
         long d_accept = hmargs.dist_accept;
         double rg_ratio = hmargs.read_ref_ratio;
         int ext_r = g_distr > -d_accept && r_distr > -d_accept && ((g_distr < (r_distr * rg_ratio)) || (g_distr < d_accept && r_distr < d_accept));
         int ext_l = g_distl > -d_accept && r_distl > -d_accept && ((g_distl < (r_distl * rg_ratio)) || (g_distl < d_accept && r_distl < d_accept));
         if (ext_r || ext_l) {
            r_min  = (ext_r ? r_distr : 0);
            l_min  = (ext_l ? r_distl : 0);
            r_qstart = match->read_e + 1;
            r_rstart = match->ref_e + 1;
            l_qstart = match->read_s - 1;
            l_rstart = match->ref_s - 1;
            last_eexp = match->e_exp;
            extend = match;
            extend_score = match->score;
         } else if (match->e_exp < hmargs.align_accept_eexp) {
            // Check overlap
            int span = align_min(seed->read_e - seed->read_s, match->read_e - match->read_s);
            int overlap = align_max(0, align_min(seed->read_e, match->read_e) - align_max(seed->read_s, match->read_s));
            overlap = overlap > span*hmargs.overlap_max_tolerance;
            int seed_ratio = seed->hits < match->hits*hmargs.align_seed_filter_thr;
            // If overlap is higher than the maximum overlap tolerance, cancel the alignment.
            if (overlap && seed_ratio) {
               cancel_align = 1;
               break;
            }
         } else {
            if (r_distr < 0) {r_distr = -r_distr; g_distr = -g_distr;}
            if (r_distl < 0) {r_distl = -r_distl; g_distl = -g_distl;}
            int same_align = (g_distr > -d_accept && (g_distr < (r_distr*rg_ratio) || g_distr < d_accept)) || (g_distl > -d_accept && (g_distl < (r_distl*rg_ratio) || g_distl < d_accept));
            if (same_align) {
               cancel_align = 1;
               break;
            }
         }
      }
      long r_qlen = slen - r_qstart;
      long l_qlen = l_qstart + 1;
      if (cancel_align || (!r_qlen && !l_qlen)) continue;
      long r_rlen = align_min((long)(r_qlen * (1 + hmargs.align.width_ratio)), index->size - r_rstart);
      long l_rlen = align_min((long)(l_qlen * (1 + hmargs.align.width_ratio)), r_qstart);
      char * r_qry = read + r_qstart;
      char * l_qry = read + l_qstart;
      char * r_ref = index->genome + r_rstart;
      char * l_ref = index->genome + l_rstart;
      
      // Align forward and backward starting from (read_s, ref_s).
      long read_s, read_e, ref_s, ref_e;
      path_t align_r = (path_t){0,0,0}, align_l = (path_t){0,0,0};
      // Forward alignment (right).
      if(r_qlen) {
         align_r = dbf_align(r_qlen, r_qry, r_rlen, r_ref, r_min, ALIGN_FORWARD, ALIGN_FORWARD, hmargs.align);
         read_e = r_qstart + align_r.row;
         ref_e = r_rstart + align_r.col;
      }
      else {
         read_e = r_qstart - 1;
         ref_e  = r_rstart - 1;
      }
      // Backward alignment (left).
      if(l_qlen) {
         align_l = dbf_align(l_qlen, l_qry, l_rlen, l_ref, l_min, ALIGN_BACKWARD, ALIGN_BACKWARD, hmargs.align);
         read_s =  l_qstart - align_l.row;
         ref_s = l_rstart - align_l.col;
      }
      else {
         read_s = l_qstart + 1;
         ref_s  = l_rstart + 1;
      }
      
      // Compute significance.
      long score = extend_score + align_l.score + align_r.score;
      double ident = 1.0 - (score*1.0)/(align_max(read_e-read_s+1,ref_e-ref_s+1));
      double e_exp = INFINITY;
      if (ident > hmargs.align_filter_ident) 
         e_exp = e_value(ref_e - ref_s + 1, extend_score + align_l.score + align_r.score, index->size);

      // If significant, store.
      if (e_exp < hmargs.align_filter_eexp) {
         match_t * hit;
         if (extend != NULL) {
            hit = extend;
            hit->hits += seed->hits;
         }
         else {
            hit = malloc(sizeof(match_t));
            hit->repeats = matchlist_new(REPEATS_SIZE);
            hit->flags  = 0;
            hit->hits   = seed->hits;
         }

         // Fill/Update hit.
         hit->score  = score;
         hit->ident  = ident;
         hit->ref_e  = ref_e;
         hit->ref_s  = ref_s;
         hit->read_e = read_e;
         hit->read_s = read_s;
         hit->e_exp  = e_exp;

         // Add to significant matchlist.
         significant = 1;
         if (extend == NULL) {
            matchlist_add(seqmatches, hit);
            matches = *seqmatches;
         }
      }
            
      free(seed);
   }

   return significant;
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
         hits[l++] = (hit_t) {.locus = get_sa(j,index->sa,index->sa_bits), .qrypos = seed.qry_pos, .depth = seed.ref_pos.depth, .bulk = seed.bulk};
      }
   }
   
   return hits;
}

int
match_seeds
(
 seedstack_t * seeds,
 int           seqlen,
 matchlist_t * matchlist,
 index_t     * index,
 int           dist_accept,
 double        read_ref_ratio
)
{
   if (matchlist->size < 1) return -1;

   // Convert seeds to hits.
   uint64_t hit_count;
   hit_t * hits = compute_hits(seeds, index, &hit_count);

   // Merge-sort hits.
   mergesort_mt(hits, hit_count, sizeof(hit_t), 0, 1, compar_hit_locus);

   // Reset matchlist.
   matchlist->pos = 0;

   // Aux variables.
   long minv = 0;
   int  min  = 0;
   size_t i = 0;

   while (i < hit_count) {
      // Match start.
      int64_t span   = hits[i].depth;
      int64_t lstart = hits[i].locus;
      int64_t lend   = lstart + span;
      int64_t rstart = hits[i].qrypos;
      int64_t rend   = rstart + span;
      int     bulk   = hits[i].bulk;
      i++;
      // Extend match.
      while (i < hit_count) {
         // Compare genome and read distances.
         int64_t g_dist = (int64_t)hits[i].locus - (int64_t)hits[i-1].locus;
         int64_t r_dist = (int64_t)hits[i].qrypos - (int64_t)hits[i-1].qrypos;
         if (r_dist > 0 && (g_dist < (read_ref_ratio * r_dist) || (g_dist < dist_accept && r_dist < dist_accept))) {
            span += hits[i].depth - (lend - (int64_t)hits[i].locus > 0 ? lend - (int64_t)hits[i].locus : 0);
            lend = hits[i].locus + hits[i].depth;
            rend = hits[i].qrypos + hits[i].depth;
            bulk *= hits[i].bulk;
            i++;
         } else break;
      }
      // Store match.
      // Non-significant streaks (bulk) will not be saved.
      if (span > minv && bulk == 0) {
         // Allocate match.
         match_t * match = malloc(sizeof(match_t));
         match->ref_e  = lend;
         match->read_e = rend;
         match->ref_s  = lstart;
         match->read_s = rstart;
         match->hits   = span;
         match->flags  = 0;
         match->score  = -1;

         // Append if list not yet full, replace the minimum value otherwise.
         if (matchlist->pos < matchlist->size) {
            matchlist->match[matchlist->pos++] = match;
         }
         else {
            match_t * match_min = matchlist->match[min];
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
   else return (ha->qrypos > hb->qrypos ? 1 : -1);
}

idxfiles_t *
index_open
(
 char * file
)
{
   // Alloc index struct.
   idxfiles_t * files = malloc(sizeof(idxfiles_t));
   
   // Read chromosome index.
   char * chr_file = malloc(strlen(file)+5);
   strcpy(chr_file, file);
   strcpy(chr_file+strlen(file), ".chr");
   files->chr = read_CHRindex(chr_file);
   if (files->chr == NULL) {
      fprintf(stderr, "error opening '%s': %s\n", chr_file, strerror(errno));      
      return NULL;
   }
   free(chr_file);

   // Open OCC file.
   char * occ_file = malloc(strlen(file)+5);
   strcpy(occ_file, file);
   strcpy(occ_file+strlen(file), ".occ");
   int fd_occ = open(occ_file, O_RDONLY);
   if (fd_occ == -1) {
      fprintf(stderr, "error opening '%s': %s\n", occ_file, strerror(errno));
      return NULL;
   }
   free(occ_file);

   // Open SA file.
   char * sa_file = malloc(strlen(file)+5);
   strcpy(sa_file, file);
   strcpy(sa_file+strlen(file), ".sar");
   int fd_sa = open(sa_file, O_RDONLY);
   if (fd_sa == -1) {
      fprintf(stderr, "error opening '%s': %s\n", sa_file, strerror(errno));
      return NULL;
   }
   free(sa_file);

   // Open GEN file.
   char * gen_file = malloc(strlen(file)+5);
   strcpy(gen_file, file);
   strcpy(gen_file+strlen(file), ".gen");
   int fd_gen = open(gen_file, O_RDONLY);
   if (fd_gen == -1) {
      fprintf(stderr, "error opening '%s': %s\n", gen_file, strerror(errno));
      return NULL;
   }
   free(gen_file);

   // Open LCP file.
   char * lcp_file = malloc(strlen(file)+5);
   strcpy(lcp_file, file);
   strcpy(lcp_file+strlen(file), ".lcp");
   int fd_lcp = open(lcp_file, O_RDONLY);
   if (fd_lcp == -1) {
      fprintf(stderr, "error opening '%s': %s\n", lcp_file, strerror(errno));
      return NULL;
   }
   free(lcp_file);

   // Open LUT file.
   /*
   char * lut_file = malloc(strlen(file)+5);
   strcpy(lut_file, file);
   strcpy(lut_file+strlen(file), ".lut");
   int fd_lut = open(lut_file, O_RDONLY);
   if (fd_lut == -1) {
      fprintf(stderr, "error opening '%s': %s\n", lut_file, strerror(errno));
      return NULL;
   }
   free(lut_file);
   */
   // Load OCC index.
   long idxsize = lseek(fd_occ, 0, SEEK_END);
   lseek(fd_occ, 0, SEEK_SET);
   files->occ_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_occ, 0);
   close(fd_occ);
   if (files->occ_file == NULL) {
      fprintf(stderr, "error mmaping .occ index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load SA index.
   idxsize = lseek(fd_sa, 0, SEEK_END);
   lseek(fd_sa, 0, SEEK_SET);
   files->sa_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_sa, 0);
   close(fd_sa);
   if (files->sa_file == NULL) {
      fprintf(stderr, "error mmaping .sar index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load GEN index.
   idxsize = lseek(fd_gen, 0, SEEK_END);
   lseek(fd_gen, 0, SEEK_SET);
   files->gen_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_gen, 0);
   close(fd_gen);
   if (files->gen_file == NULL) {
      fprintf(stderr, "error mmaping .gen index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load LCP index.
   idxsize = lseek(fd_lcp, 0, SEEK_END);
   lseek(fd_lcp, 0, SEEK_SET);
   files->lcp_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_lcp, 0);
   close(fd_lcp);
   if (files->lcp_file == NULL) {
      fprintf(stderr, "error mmaping .lcp index file: %s.\n", strerror(errno));
      return NULL;
   }
   // Load LUT index.
   /*
   idxsize = lseek(fd_lut, 0, SEEK_END);
   lseek(fd_lut, 0, SEEK_SET);
   files->lut_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_lut, 0);
   close(fd_lut);
   if (files->lut_file == NULL) {
      fprintf(stderr, "error mmaping .lut index file: %s.\n", strerror(errno));
      return NULL;
   }
   */
   return files;

}

index_t *
index_format
(
 idxfiles_t * files
)
{
   if (files == NULL) return NULL;
   // Alloc index struct.
   index_t * index = malloc(sizeof(index_t));
   // GEN.
   index->genome = (char *) files->gen_file;
   // OCC.
   index->occ_mark_int = *(uint64_t *) files->occ_file;
   index->c = ((uint64_t *) files->occ_file + 1);
   index->occ = index->c + NUM_BASES + 1;
   // Genome size.
   index->size = index->c[NUM_BASES]; // Count forward strand only.
   // SAR.
   index->sa_bits = *((uint64_t *) files->sa_file);
   index->sa      = ((uint64_t *) files->sa_file) + 1;
   // LUT.
   //   index->lut_kmer = *((uint64_t *) files->lut_file);
   //   index->lut = ((uint64_t *) files->lut_file) + 1;
   // LCP.
   // params
   index->lcp_mark_int = *((uint64_t *)files->lcp_file);
   index->lcp_min_depth =  *((uint64_t *)files->lcp_file + 1);
   uint64_t lcp_idx_size = *((uint64_t *)files->lcp_file + 2);//((2*index->size + LCP_WORD_SIZE - 1)/LCP_WORD_SIZE + index->lcp_mark_int - 1)/index->lcp_mark_int * (index->lcp_mark_int + 1) + 1;
   // index
   index->lcp_sample_idx = ((uint64_t *) files->lcp_file + 3);
   uint64_t ext_idx_size = *(index->lcp_sample_idx + lcp_idx_size);
   index->lcp_extend_idx = index->lcp_sample_idx + lcp_idx_size + 1;
   // alloc structures
   index->lcp_sample = malloc(sizeof(lcpdata_t));
   if (index->lcp_sample == NULL) return NULL;
   index->lcp_extend = malloc(sizeof(list64_t));
   if (index->lcp_extend == NULL) return NULL;
   // data
   index->lcp_sample->size = *(index->lcp_extend_idx + ext_idx_size)/2;
   index->lcp_sample->lcp = (lcpval_t *)(index->lcp_extend_idx + ext_idx_size + 1);
   uint64_t * lcpext_size = (uint64_t *)(index->lcp_sample->lcp + index->lcp_sample->size);
   index->lcp_extend->size = *lcpext_size;
   index->lcp_extend->val = (int64_t *)(lcpext_size + 1);
   //CHR.
   index->chr = files->chr;
   

   return index;
}

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


seqstack_t *
read_file
(
   FILE      * inputf
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
      format = FASTA;
      break;
   case '@':
      format = FASTQ;
      break;
   default:
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
   char *q    = NULL;
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
            q = temp;
         } else continue;
         break;
      }

      /*size_t seqlen = strlen(seq);
      if (seqlen > MAXSEQLEN) {
         fprintf(stderr, "max sequence length exceeded (%d)\n", MAXSEQLEN);
         fprintf(stderr, "offending sequence:\n%s\n", seq);
         continue;
         }*/
      
      
      if (seq_push(&seqstack, tag, seq, q)) return NULL;
   }

   free(line);
   free(subline);
   free(temp);
   return seqstack;
}

int
compar_seedhits
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);
   
   if (ma->hits < mb->hits) return 1;
   else return -1;
}

int
compar_matcheexp
(
 const void * a,
 const void * b,
 const int   param
)
{
   match_t * ma = *((match_t **) a);
   match_t * mb = *((match_t **) b);

   if (mb->e_exp < ma->e_exp) return 1;
   else return -1;
}
