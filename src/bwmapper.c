#include "bwmapper.h"

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

  // Backtrace handler
   signal(SIGSEGV, SIGSEGV_handler); 

   if (argc < 3) {
      fprintf(stderr, "usage: bwmapper {query <seq> | index} <genome file>\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1],"query") == 0) {
      // Parse filename.
      char * filename = malloc(strlen(argv[3] + 5));
      strcpy(filename, argv[3]);
      // Read index.
      long     gsize;
      long   * c;
      long   * pos;
      list_t * occ;
      strcpy(filename+strlen(argv[3]), ".fmi");
      read_FMindex(filename, &gsize, &c, &pos, &occ);
      strcpy(filename+strlen(argv[3]), ".index");
      chr_t    chr   = read_CHRindex(filename);

      // Query index.
      vstack_t * hits;
      hits = query_index(argv[2], gsize, c, pos, occ);
      
      fprintf(stdout, "[Query: %s]\n", argv[2]);
      for (long i = 0; i < hits->pos; i++) {
         int chrnum = bisect_search(0, chr.nchr-1, chr.start, hits->val[i]+1)-1;
         fprintf(stdout, "  %s:%ld\n", chr.name[chrnum], hits->val[i]-chr.start[chrnum]+1);
      }
      fprintf(stdout,"\nTotal: %ld matches.\n", hits->pos);
   }
   else if (strcmp(argv[1],"index") == 0) {
      write_index(argv[2]);
   }
   else 
      fprintf(stderr, "usage: bwmapper {query <seq> | index} <genome file>\n");
   
   return 0;
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
   stack->pos  = 0;
   stack->size = size;
   
   return stack;
}


void
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
         exit(EXIT_FAILURE);
      }
      stack->size = newsize;
   }

   stack->val[stack->pos++] = value;
}



/*********************/
/** query functions **/
/*********************/


chr_t
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

   // Read chromosome entries.
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      char *name;
      int i = 0;
      while (buffer[i] != '\t' && buffer[i] != '\n') i++;
      if (buffer[i] == '\n') {
         fprintf(stderr, "detected corrupted value in index file, ignoring chromosome: %s\n", buffer);
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
   chr_t chrindex;
   chrindex.nchr = chrcount;
   chrindex.start= start;
   chrindex.name = names;
   return chrindex;
}

long *
read_FMindex
(
 char   *  filename,
 long   *  gsize,
 long   ** Cp,
 long   ** posp,
 list_t ** occp
)
{
   // mmap index.
   int fd = open(filename, O_RDONLY);
   if (fd == -1) {
      fprintf(stderr, "error opening '%s': %s\n", filename, strerror(errno));
      exit(EXIT_FAILURE);
   }

   long idxsize = lseek(fd, 0, SEEK_END);
   lseek(fd, 0, SEEK_SET);
   long * index = mmap(NULL, idxsize, PROT_READ, MAP_PRIVATE, fd, 0);

   // index structure: (8 byte alignment)
   // DATA: C[NUM_BASES]:gsize:POS[gsize]:ooc(1)[sz_occ(1)]:...:occ(NUM_BASES)[sz_occ(NUM_BASES)]
   // OCC:  pos:size:val[0]:val[1]:val[2]:...:val[pos-1]
   long pointer_offset = 0;
   long * c = *Cp = malloc((NUM_BASES+1)*sizeof(long));
   for (int i = 0; i < NUM_BASES; i++) c[i] = index[i];
   pointer_offset += NUM_BASES;

   *gsize = index[pointer_offset++];
   c[NUM_BASES] = *gsize;
   
   *posp = index + pointer_offset;
   pointer_offset += *gsize;

   list_t * occ = *occp = malloc(NUM_BASES * sizeof(list_t));
   for (int i = 0; i < NUM_BASES; i++) {
      occ[i].max = index[pointer_offset++];
      occ[i].val = index + pointer_offset;
      pointer_offset += occ[i].max;
   }

   return index;
}


vstack_t *
query_index
(
 char  * query,
 long    gsize,
 long  * c,
 long  * pos,
 list_t * occs
)
{
   long sp, ep;
   int * qval = translate_query(query);
   // Initialize range.
   sp = c[qval[0]];
   ep = c[qval[0] + 1] - 1;

   // Iterate over FM-index.
   for (int i = 1; i < strlen(query) && sp <= ep; i++) {
      int nt = qval[i];
      long occsp =  bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1);
      long occep = bisect_search(0, occs[nt].max-1, occs[nt].val, ep);
      sp = c[nt] + (occs[nt].max ? occsp : 0);
      ep = c[nt] + (occs[nt].max ? occep : 0 ) - 1;
   }
   if (ep < sp) return new_stack(0);

   // Save hits.
   vstack_t * hits = new_stack((ep-sp)+1);
   hits->pos = ep-sp+1;
   for (int i = 0; i < hits->pos; i++) {
      hits->val[i] = pos[sp+i];
   }

   return hits;
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
   while(s < NUM_BASES*sizeof(long)) s += write(fd, C + s/sizeof(long), NUM_BASES*sizeof(long) - s);
   s = 0;
   while(s < sizeof(long)) s += write(fd, &gsize, sizeof(long));
   s = 0;
   while(s < gsize*sizeof(long)) s += write(fd, pos + s/sizeof(long), gsize*sizeof(long) - s);
   for (int i = 0; i < NUM_BASES; i++) {
      fprintf(stdout, "occ[%d] size = %ld (%ld bytes).\n", i, occ[i]->pos, occ[i]->pos*sizeof(long));
      s = 0;
      while(s < sizeof(long)) s += write(fd, &(occ[i]->pos), sizeof(long));
      fprintf(stdout, "occ[%d]->pos: written %ld bytes.\n", i, s);
      s = 0;
      while(s < occ[i]->pos * sizeof(long)) s += write(fd,((long*) &(occ[i]->val[0])) + s/sizeof(long), occ[i]->pos*sizeof(long) - s);
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
   
   genome = realloc(genome, *genomesize+1);
   genome[*genomesize] = 0;
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


