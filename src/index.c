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

int
ann_index_write
(
 annlist_t * list,
 char * index_file
)
{
   // File name.
   char * fname = add_suffix(index_file, ".ann");
   if (fname == NULL) return -1;
   // Open file to write.
   int fd = open(fname,O_WRONLY | O_CREAT | O_TRUNC,0644);
   free(fname);
   if (fd == -1) {
      fprintf(stderr, "error opening '%s': %s\n", fname, strerror(errno));
      return -1;
   }
   // Write table.
   if (write(fd, &(list->count), 1) < 1) {
      fprintf(stderr,"[error] could not write annotation index.\n");
      return -1;
   }
   size_t bytes = 0, total = list->count*sizeof(ann_t);
   while ((bytes += write(fd, ((uint8_t *)list->ann) + bytes, total-bytes)) < total);

   return 0;
}

annlist_t *
ann_index_read
(
 char * index_file
)
{
   // File name.
   char * fname = add_suffix(index_file, ".ann");

   if (fname == NULL) return NULL;
   
   // If file does not exist, create a new annotation index.
   if (access(fname, F_OK)) {
      annlist_t * list = malloc(sizeof(annlist_t) + 10 * sizeof(ann_t));
      list->count = 0;
      return list;
   }
   // Open file.
   int fd = open(fname, O_RDONLY);
   free(fname);
   if (fd == -1) {
      fprintf(stderr, "error opening '%s': %s\n", fname, strerror(errno));
      return NULL;
   }
   // Read table size.
   uint8_t count = 0;
   if (read(fd, &count, sizeof(uint8_t)) < 1) {
      fprintf(stderr, "[error] could not read data from annotation index.\n");
      return NULL;
   }
   // Alloc structure.
   annlist_t * list = malloc(sizeof(annlist_t) + count * sizeof(ann_t));
   if (list == NULL) return NULL;
   // Read headers.
   list->count = count;
   size_t bytes = 0, total = count*sizeof(ann_t);
   while ((bytes += read(fd,((uint8_t*)list->ann)+bytes,total-bytes)) < total);

   // Sort annotations by ascending k.
   qsort(list->ann, list->count, sizeof(ann_t), compar_ann_k_asc);
   return list;
}

int
ann_index_load
(
 annlist_t * index,
 char      * index_file
)
{
   for (int i = 0; i < index->count; i++) {
      // Generate file names.
      char * suffix = malloc(9);
      if (suffix == NULL)
         return -1;
      sprintf(suffix, ".ann.%d", index->ann[i].id);
      char * fname = add_suffix(index_file, suffix);
      if (fname == NULL)
         return -1;
      // Open file.
      int fd = open(fname, O_RDONLY);
      if (fd == -1) {
         fprintf(stderr, "[error] opening '%s': %s\n", fname, strerror(errno));
         return -1;
      }
      // mmap.
      size_t f_len = lseek(fd, 0, SEEK_END);
      lseek(fd, 0, SEEK_SET);
      index->ann[i].data = mmap(NULL, f_len, PROT_READ, MMAP_FLAGS, fd, 0);
      if (index->ann[i].data == NULL) {
         fprintf(stderr, "[error] mmaping '%s' index file: %s.\n", fname, strerror(errno));
         return -1;
      }
      close(fd);
      free(fname);
      free(suffix);
   }
      
   return 0;
}

int
ann_find_slot
(
 annlist_t * list
)
{
   uint8_t * idlist = calloc(257,1);
   for (int i = 0; i < list->count; i++)
      idlist[list->ann[i].id] = 1;
   int id = 0;
   while (idlist[id]) id++;
   return id;
}

void
ann_print_index
(
 annlist_t * list
)
{
   fprintf(stderr,"[info] annotation index content:\n");
   for (int i = 0; i < list->count; i++) {
      fprintf(stderr, "[info] {%d} k:%d, d:%d\n", list->ann[i].id, list->ann[i].k, list->ann[i].d);
   }
   fprintf(stderr, "[info] %d annotations.\n", list->count);
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
   int ann_slot = -1, ann_pos = -1;
   annlist_t * annlist = NULL;
   // Check current annotations.
   annlist = ann_index_read(index_file);
   ann_print_index(annlist);
   if (annlist == NULL)
      return EXIT_FAILURE;
   for (int i = 0; i < annlist->count; i++) {
      if (kmer == annlist->ann[i].k) {
         if (tau <= annlist->ann[i].d) {
            fprintf(stderr, "[warning] the index already has a (%d,%d) annotation.\n",kmer,tau);
            return EXIT_SUCCESS;
         } else {
            ann_slot = annlist->ann[i].id;
            ann_pos  = i;
            fprintf(stderr, "[warning] the existing (%d,%d) annotation will be overwritten (id=%d).\n",kmer,annlist->ann[i].d,ann_slot);
            break;
         }
      }
   }
   // Find slot to store annotation.
   if (ann_slot < 0) ann_slot = ann_find_slot(annlist);

   fprintf(stderr, "[info] computing (%d,%d) genomic neighborhood using %d threads.\n", kmer, tau, threads);
   annotation_t ann = annotate(kmer,tau,index,threads);

   // Write output file.
   // Generate file name.
   char * suffix = malloc(9);
   if (suffix == NULL)
      return EXIT_FAILURE;
   sprintf(suffix, ".ann.%d", ann_slot);
   char * fname = add_suffix(index_file, suffix);
   if (fname == NULL)
      return EXIT_FAILURE;
   free(suffix);
   // Open file.
   int fd = open(fname,O_WRONLY | O_CREAT | O_TRUNC,0644);
   free(fname);
   size_t bytes = 0;
   // Write annotation.
   fprintf(stderr,"[info] writing annotation (%ld bytes)... ", ann.size);
   while ((bytes += write(fd,ann.info+bytes,ann.size-bytes)) < ann.size);
   fprintf(stderr,"%ld bytes written.\n",bytes);
   // Close file descriptor.
   close(fd);
   // Update table.
   fprintf(stderr,"[info] updating annotation index... ");
   annlist = realloc(annlist, sizeof(annlist_t) + (annlist->count + 1)*sizeof(ann_t));
   if (ann_pos < 0) ann_pos = annlist->count++;
   annlist->ann[ann_pos] = (ann_t) {
      .id = ann_slot,
      .k = kmer,
      .d = tau,
      .size = ann.size,
      .data = NULL
   };
   ann_index_write(annlist, index_file);
   free(annlist);
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
      if (list->ann[i].k > k)
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
   return (a->k > b->k ? 1 : -1);
}
