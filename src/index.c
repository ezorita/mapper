#include "index.h"

char *
add_suffix
(
 const char * string,
 const char * suffix
)
{
   char * new = malloc(strlen(string) + strlen(suffix) + 1);
   if (new == NULL) return NULL;
   strcpy(new,string);
   strcpy(new+strlen(string),suffix);
   return new;
}

annlist_t *
ann_index_read
(
 char * index_file
)
{
   // Annotation files
   char * pattern = add_suffix(index_file, ".ann.*");
   glob_t gbuf;
   glob(pattern,GLOB_TILDE,NULL,&gbuf);
   free(pattern);

   // Alloc list.
   annlist_t * list = malloc(sizeof(annlist_t) + gbuf.gl_pathc * sizeof(ann_t));
   list->count = 0;
   
   // Iterate matching files.
   for (int i = 0; i < gbuf.gl_pathc; i++) {
      // Open file.
      char * fname = gbuf.gl_pathv[i];
      fprintf(stderr,"[info] loading glob result: '%s'.\n",fname);
      int fd = open(fname, O_RDONLY);
      if (fd == -1) {
         fprintf(stderr, "[error] opening '%s': %s\n", fname, strerror(errno));
         free(list);
         return NULL;
      }
      // mmap.
      struct stat sb;
      fstat(fd, &sb);
      if (sb.st_size <= 24) {
         fprintf(stderr, "[error] Found empty annotation in '%s', file size is %ld bytes. Ignoring file.\n",fname,sb.st_size);
         continue;
      }
      anndata_t * data = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);

      if (list->ann[i].data == NULL) {
         fprintf(stderr, "[error] mmaping '%s' index file: %s.\n", fname, strerror(errno));
         free(list);
         return NULL;
      }

      // Check file content.
      if (data->magic != ANN_MAGICNO) {
         fprintf(stderr, "[error] Wrong magic numer for '%s'. Ignoring file.\n",fname);
         continue;
      }
      list->ann[list->count++] = (ann_t){strdup(fname), data};
      close(fd);
   }
   // Sort annotations by ascending k.
   qsort(list->ann, list->count, sizeof(ann_t), compar_ann_k_asc);
   return list;
}


char *
ann_new_filename
(
 int    k,
 int    d,
 char * index_file
)
{
   // Default suffix.
   char suffix[12];
   sprintf(suffix, ".ann.%d%d",k,d);
   char * path = add_suffix(index_file,suffix);
   
   // Check whether file exists.
   int i = 0;
   while (access(path,F_OK) == 0) {
      free(path);
      sprintf(suffix, ".ann.%d%d.%d",k,d,i++);
      path = add_suffix(index_file,suffix);
   }

   return path;
}

void
ann_print_index
(
 annlist_t * list
)
{
   if (! list->count) {
      fprintf(stderr,"[info] 0 existing annotations found\n");
   } else {
      fprintf(stderr,"[info] existing annotations:\nk\td\tpath\n");
      for (int i = 0; i < list->count; i++) {
         fprintf(stderr, "%d\t%d\t%s\n", list->ann[i].data->kmer, list->ann[i].data->tau, list->ann[i].file);
      }
      fprintf(stderr, "[info] total: %d annotations.\n", list->count);
   }
}


int
index_add_annotation
(
 int       kmer,
 int       tau,
 int       threads,
 index_t * index,
 char    * index_file
)
{
   annlist_t * annlist = NULL;
   // Load annotation files.
   annlist = ann_index_read(index_file);
   if (annlist == NULL)
      return EXIT_FAILURE;
   
   // Print index.
   ann_print_index(annlist);

   // Check existing annotations.
   for (int i = 0; i < annlist->count; i++) {
      if (kmer == annlist->ann[i].data->kmer) {
         if (tau <= annlist->ann[i].data->tau) {
            fprintf(stderr, "[warning] the index has a (%d,%d) annotation in '%s'.\n",kmer,tau,annlist->ann[i].file);
         }
      }
   }
   free(annlist);

   // Generate file name.
   char * fname  = ann_new_filename(kmer,tau,index_file);
   if (fname == NULL)
      return EXIT_FAILURE;

   // Open file.
   fprintf(stderr,"[info] opening annotation file: '%s'.\n", fname);
   int fd = open(fname,O_WRONLY | O_CREAT | O_TRUNC,0644);
   if (fd == -1) {
      fprintf(stderr, "error opening file to write: '%s'\n", fname);
      return EXIT_FAILURE;
   }
   free(fname);

   // Compute neighbors.
   fprintf(stderr, "[info] computing (%d,%d) genomic neighborhood using %d threads.\n", kmer, tau, threads);
   annotation_t ann = annotate(kmer,tau,index,threads);

   // Write header. (magic number, k, tau, size).
   size_t bytes = 0;
   uint64_t magicno = ANN_MAGICNO;
   bytes = write(fd,&magicno,sizeof(uint64_t));
   bytes = write(fd,&(ann.kmer),sizeof(uint32_t));
   bytes = write(fd,&(ann.tau),sizeof(uint32_t));
   bytes = write(fd,&(ann.size),sizeof(size_t));

   // Write annotation.
   fprintf(stderr,"[info] writing annotation (%ld bytes)... ", ann.size);
   bytes = 0;
   while ((bytes += write(fd,ann.info+bytes,ann.size-bytes)) < ann.size);
   fprintf(stderr,"%ld bytes written.\n",bytes);

   // Close file descriptor.
   close(fd);

   fprintf(stderr,"done.\n");

   return EXIT_SUCCESS;
}

index_t *
index_load_base
(
 char * index_file
)
{
   // Alloc index struct.
   index_t * index = malloc(sizeof(index_t));
   if (index == NULL) return NULL;
   // Load genome bases.
   index->genome = index_load_gen(index_file);
   if (index->genome == NULL) return NULL;
   // Load bwt table.
   index->bwt = index_load_bwt(index_file);
   if (index->bwt == NULL) return NULL;
   // Load sar.
   index->sar = index_load_sar(index_file);
   if (index->sar == NULL) return NULL;
   // Load chr list.
   index->chr = index_load_chr(index_file);
   if (index->chr == NULL) return NULL;
   // Set genome size.
   index->size = index->bwt->c[NUM_BASES];

   return index;
}

char *
index_load_gen
(
 char * index_file
)
{
   // Open GEN file.
   char * gen_file = add_suffix(index_file, ".gen");
   int fd_gen = open(gen_file, O_RDONLY);
   if (fd_gen == -1) {
      fprintf(stderr, "error opening '%s': %s\n", gen_file, strerror(errno));
      return NULL;
   }
   free(gen_file);

   // Load GEN index.
   size_t gen_len = lseek(fd_gen, 0, SEEK_END);
   lseek(fd_gen, 0, SEEK_SET);
   void * gen_map = mmap(NULL, gen_len, PROT_READ, MMAP_FLAGS, fd_gen, 0);
   close(fd_gen);
   if (gen_file == NULL) {
      fprintf(stderr, "error mmaping .gen index file: %s.\n", strerror(errno));
      return NULL;
   }

   return (char *) gen_map;
}

bwt_t *
index_load_bwt
(
 char * index_file
)
{
   bwt_t * bwt = malloc(sizeof(bwt_t));
   if (bwt == NULL) return NULL;
   // Open OCC file.
   char * bwt_file = add_suffix(index_file, ".bwt");
   if (bwt_file == NULL) return NULL;
   int fd_bwt = open(bwt_file, O_RDONLY);
   if (fd_bwt == -1) {
      fprintf(stderr, "error opening '%s': %s\n", bwt_file, strerror(errno));
      return NULL;
   }
   free(bwt_file);

   // Load BWT index.
   size_t bwt_len = lseek(fd_bwt, 0, SEEK_END);
   lseek(fd_bwt, 0, SEEK_SET);
   void * bwt_map = mmap(NULL, bwt_len, PROT_READ, MMAP_FLAGS, fd_bwt, 0);
   close(fd_bwt);
   if (bwt_map == NULL) {
      fprintf(stderr, "error mmaping .bwt index file: %s.\n", strerror(errno));
      return NULL;
   }

   // Fill OCC struct.
   bwt->occ_mark_int = *(uint64_t *) bwt_map;
   bwt->c = ((uint64_t *) bwt_map + 1);
   bwt->occ = bwt->c + NUM_BASES + 1;

   // Base positions.
   bwt->fmd_base = (fmdpos_t){.fp = 0, .rp = 0, .sz = bwt->c[NUM_BASES], .dp = 0};
   bwt->bwt_base = (bwpos_t){.sp = 0, .ep = bwt->c[NUM_BASES] - 1, .depth = 0};

   return bwt;
}

sar_t *
index_load_sar
(
 char * index_file
)
{
   sar_t * sar = malloc(sizeof(sar_t));
   if (sar == NULL) return NULL;
   // Open SAR file.
   char * sa_file = add_suffix(index_file,".sar");
   if (sa_file == NULL) return NULL;
   int fd_sa = open(sa_file, O_RDONLY);
   if (fd_sa == -1) {
      fprintf(stderr, "error opening '%s': %s\n", sa_file, strerror(errno));
      return NULL;
   }
   free(sa_file);
 
   // Load SAR index.
   size_t sa_len = lseek(fd_sa, 0, SEEK_END);
   lseek(fd_sa, 0, SEEK_SET);
   void * sa_map = mmap(NULL, sa_len, PROT_READ, MMAP_FLAGS, fd_sa, 0);
   close(fd_sa);
   if (sa_map == NULL) {
      fprintf(stderr, "error mmaping .sar index file: %s.\n", strerror(errno));
      return NULL;
   }

   // Fill SAR struct.
   sar->bits = *((uint64_t *) sa_map);
   sar->sa   =  ((uint64_t *) sa_map) + 1;

   return sar;
}

chr_t *
index_load_chr
(
 char * index_file
)
{

   // Read chromosome index.
   char * chr_file = add_suffix(index_file,".chr");
   if (chr_file == NULL) return NULL;
   // Files
   FILE * input = fopen(chr_file,"r");
   if (input == NULL) {
      fprintf(stderr, "error in 'read_CHRindex' (fopen): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   free(chr_file);

   int chrcount = 0;
   int structsize = CHRSTR_SIZE;

   long  * start = malloc(structsize*sizeof(long));
   char ** names = malloc(structsize*sizeof(char*));
   char  * buffer = malloc(BUFFER_SIZE);
   if (start == NULL || names == NULL || buffer == NULL) return NULL;

   // File read vars
   ssize_t rlen;
   size_t  sz = BUFFER_SIZE;
   int lineno = 0;

   // Read chromosome entries.
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      lineno++;
      char *name;
      int i = 0;
      while (buffer[i] != '\t' && buffer[i] != '\n') i++;
      if (buffer[i] == '\n') {
         fprintf(stderr, "illegal format in %s (line %d) - ignoring chromosome: %s\n", index_file, lineno, buffer);
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
      if (names[chrcount] == NULL) return NULL;
      strcpy(names[chrcount], name);
      // Inc.
      chrcount++;
   }

   // Return chromosome index structure.
   chr_t * chrindex = malloc(sizeof(chr_t));
   if (chrindex == NULL) return NULL;
   chrindex->nchr = chrcount;
   chrindex->start= start;
   chrindex->name = names;

   // Close files.
   fclose(input);

   // Free memory.
   free(buffer);

   return chrindex;
}

int
ann_read
(
 ann_t      ann,
 uint64_t   locus,
 int      * cnt
)
{
   return 0;
/* TO COMPILE
   int bits = 0;
   while ((ann.d >> bits) > 0) bits++;
   bits += 2;
   uint16_t mask = 0xFFFF;
   mask >>= (16-bits);
   uint64_t bit = bits*locus;
   uint64_t word = bit >> 3;
   uint8_t shift = bit & 7;
   uint16_t w = (((uint16_t)ann.data[word]) >> shift) | (((uint16_t)ann.data[word+1]) << (8-shift));
   int d        = w & (mask >> 2);
   int log10cnt = (w >> (bits-2)) & 3;
   if (d == 0 && log10cnt == 3) {
      log10cnt = 0;
      d = ann.d+1;
   }
   if (cnt != NULL) *cnt = log10cnt;
   return d;
*/
}


ann_t *
ann_find
(
 int         k,
 annlist_t * list
)
{
   if (list->count == 0) return NULL;
   int idx = -1;
   for (int i = 0; i < list->count; i++) {
      if (list->ann[i].data->kmer > k)
         break;
      else
         idx = i;
   }
   if (idx < 0) return NULL;
   return list->ann + idx;
}


int
compar_ann_k_asc
(
 const void * aa,
 const void * ab
)
{
   ann_t * a = (ann_t *)aa;
   ann_t * b = (ann_t *)ab;
   return (a->data->kmer > b->data->kmer ? 1 : -1);
}
