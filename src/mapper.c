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

   // Test seed function.

   clock_t t = clock();
   seedopt_t opt = {.min_len = 1, .max_len = 100, .min_loci = 1, .max_loci = 1}; // SMEM
   int nseeds = 0;
   for (int i = 0; i < seqs->pos; i++) {
      // Backward-Shrink SMEMS.
      seedstack_t * seeds = seed(seqs->seq[i].seq, opt, index);
      // Naive SMEMS.
      seedstack_t * nvseeds = naive_smem(seqs->seq[i].seq, opt, index);

      if (seeds->pos != nvseeds->pos) {
         fprintf(stdout,"different seed number.\n");
         break;
      }
      for (int i = 0; i < seeds->pos; i++) {
         seed_t ns = nvseeds->seed[i];
         seed_t s = seeds->seed[i];
         if (s.qry_pos != ns.qry_pos || s.ref_pos.sp != ns.ref_pos.sp || s.ref_pos.ep != ns.ref_pos.ep || s.ref_pos.depth != ns.ref_pos.depth) {
            fprintf(stdout,"different seed info.\n");
            break;
         }
      } 

      free(seeds);
      free(nvseeds);
      /*
        for (int i = 0; i < seeds->pos; i++) {
        seed_t seed = seeds->seed[i];
        fprintf(stdout, "seed %d: qry_pos = %d, ref_pos = %ld,%ld,%d\n", i, seed.qry_pos, seed.ref_pos.sp, seed.ref_pos.ep, seed.ref_pos.depth);
        }
      */
   }
   fprintf(stdout, "seeds: %d (%.3fms)\n", nseeds, (clock()-t)*1000.0/CLOCKS_PER_SEC);

   /*
   bwpos_t pos = {.depth = 0, .sp = 0, .ep = index->size-1};
   // Extend 'AGCG' then shrink.
   suffix_extend(2, pos, &pos, index); // 16,22
   suffix_extend(1, pos, &pos, index); // 10,11
   suffix_extend(2, pos, &pos, index); // 17,17
   suffix_extend(0, pos, &pos, index); // 3,3
   // Pos should be (4,3,3)
   suffix_shrink(pos, &pos, index); // 3,3,5
   // Pos should be (3,3,5)

   pos = (bwpos_t) {.depth = 0, .sp = 0, .ep = index->size-1};
   // Extend 'CTAGT' then shrink.
   suffix_extend(3, pos, &pos, index); // 23,30
   suffix_extend(2, pos, &pos, index); // 21,22
   suffix_extend(0, pos, &pos, index); // 6,6
   suffix_extend(3, pos, &pos, index); // 27,27
   suffix_extend(1, pos, &pos, index); // 15,15
   // Pos should be (5,15,15)
   suffix_shrink(pos, &pos, index); // 4,12,15
   // Pos should be (4,12,15)

   pos = (bwpos_t) {.depth = 0, .sp = 0, .ep = index->size-1};
   // Extend 'TAGCT' then shrink twice.
   suffix_extend(3, pos, &pos, index); // 22,29
   suffix_extend(1, pos, &pos, index); // 11,14
   suffix_extend(2, pos, &pos, index); // 17,19
   suffix_extend(0, pos, &pos, index); // 3,4
   suffix_extend(3, pos, &pos, index); // 24,25
   // Pos should be (24,25,5).
   suffix_shrink(pos, &pos, index); // 23,25
   suffix_shrink(pos, &pos, index); // 23,26
   // Pos should be (23,26,3).

   pos = (bwpos_t) {.depth = 0, .sp = 0, .ep = index->size-1};
   // Extend 'CTAGCTAGC' then shrink.
   suffix_extend(1, pos, &pos, index); // 8,14
   suffix_extend(2, pos, &pos, index); // 16,19
   suffix_extend(0, pos, &pos, index); // 2,4
   suffix_extend(3, pos, &pos, index); // 23,25
   suffix_extend(1, pos, &pos, index); // 11,13
   suffix_extend(2, pos, &pos, index); // 17,18
   suffix_extend(0, pos, &pos, index); // 3,3
   suffix_extend(3, pos, &pos, index); // 24,24
   suffix_extend(1, pos, &pos, index); // 12,12
   // Pos should be (12,12,9).
   suffix_shrink(pos, &pos, index); // 12,13,8
   suffix_shrink(pos, &pos, index); // 11,13,5
   suffix_shrink(pos, &pos, index); // 11,14,4
   suffix_shrink(pos, &pos, index); // 8.14,1
   // Pos should be (8,14,1).

   pos = (bwpos_t) {.depth = 0, .sp = 0, .ep = index->size-1};
   // Extend 'GCTAGCTA' then shrink.
   suffix_extend(0, pos, &pos, index); // 0,7
   suffix_extend(3, pos, &pos, index); // 23,27
   suffix_extend(1, pos, &pos, index); // 11,14
   suffix_extend(2, pos, &pos, index); // 17,19
   suffix_extend(0, pos, &pos, index); // 3,4
   suffix_extend(3, pos, &pos, index); // 24,25
   suffix_extend(1, pos, &pos, index); // 12,13
   suffix_extend(2, pos, &pos, index); // 18,18
   // Pos should be (12,12,9).
   suffix_shrink(pos, &pos, index); // 17,18,6
   suffix_shrink(pos, &pos, index); // 17,19,5
   suffix_shrink(pos, &pos, index); // 16,19,2
   suffix_shrink(pos, &pos, index); // 15,21,1
   // Pos should be (15,21,1).
   suffix_shrink(pos, &pos, index); // 0,29,0
   */
   return 0;
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
   char * lut_file = malloc(strlen(file)+5);
   strcpy(lut_file, file);
   strcpy(lut_file+strlen(file), ".lut");
   int fd_lut = open(lut_file, O_RDONLY);
   if (fd_lut == -1) {
      fprintf(stderr, "error opening '%s': %s\n", lut_file, strerror(errno));
      return NULL;
   }
   free(lut_file);

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
   idxsize = lseek(fd_lut, 0, SEEK_END);
   lseek(fd_lut, 0, SEEK_SET);
   files->lut_file = mmap(NULL, idxsize, PROT_READ, MMAP_FLAGS, fd_lut, 0);
   close(fd_lut);
   if (files->lut_file == NULL) {
      fprintf(stderr, "error mmaping .lut index file: %s.\n", strerror(errno));
      return NULL;
   }

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
   index->lut_kmer = *((uint64_t *) files->lut_file);
   index->lut = ((uint64_t *) files->lut_file) + 1;
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
   index->lcp_extend = malloc(sizeof(list32_t));
   if (index->lcp_extend == NULL) return NULL;
   // data
   index->lcp_sample->size = *(index->lcp_extend_idx + ext_idx_size)/2;
   index->lcp_sample->lcp = (lcpval_t *)(index->lcp_extend_idx + ext_idx_size + 1);
   uint64_t * lcpext_size = (uint64_t *)(index->lcp_sample->lcp + index->lcp_sample->size);
   index->lcp_extend->size = *lcpext_size;
   index->lcp_extend->val = (int32_t *)(lcpext_size + 1);
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
