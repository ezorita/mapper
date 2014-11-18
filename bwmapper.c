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
      fprintf(stderr, "usage: bwmapper {read | write}  <file>\n");
      exit(EXIT_FAILURE);
   }

   if (strcmp(argv[1],"query") == 0) {
      // Read index.
      long    gsize;
      long  * index;
      long  * c;
      long  * pos;
      list_t * occ;
      index = read_index(argv[3], &gsize, &c, &pos, &occ);

      // Query index.
      vstack_t * hits;
      hits = query_index(argv[2], gsize, c, pos, occ);
      
      fprintf(stdout, "[Query: %s]", argv[2]);
      if (hits == NULL) fprintf(stdout, " not found.\n");
      else {
         for (int i = 0; i < hits->pos; i++) fprintf(stdout, " %ld", hits->val[i]);
         fprintf(stdout,"\n");
      }
   }
   else if (strcmp(argv[1],"index") == 0) {
      write_index(argv[2]);
   }
   else 
      fprintf(stderr, "usage: bwmapper {query <query> <file> | index <genomefile>}\n");
   
   return 0;
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
   sp = c[qval[0]];
   ep = c[qval[0] + 1] - 1;
   //   fprintf(stdout, "[iteration 0 %d]\tsp=%ld\tep=%ld\n",qval[0],sp,ep);
   for (int i = 1; i < strlen(query) && sp <= ep; i++) {
      int nt = qval[i];
      sp = c[nt] + (occs[nt].max ? bisect_search(0, occs[nt].max-1, occs[nt].val, sp-1) : 0);
      ep = c[nt] + (occs[nt].max ? bisect_search(0, occs[nt].max-1, occs[nt].val, ep) : 0 ) - 1;
      // fprintf(stdout, "[iteration %d %d]\tsp=%ld\tep=%ld\n",i,qval[i],sp,ep);
   }
   if (ep < sp) return NULL;
   vstack_t * hits = new_stack((ep-sp)+1);
   hits->pos = ep-sp+1;
   for (int i = 0; i < hits->pos; i++) {
      hits->val[i] = pos[sp+i];
   }
   return hits;
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

int *
translate_query
(
 char * query
)
{
   int   qlen = strlen(query);
   int * qval = malloc(qlen * sizeof(int));
   for (int i = 0; i < qlen; i++) qval[i] = translate[query[qlen - i - 1]];
   return qval;
}


long *
read_index
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


   long * C = *Cp;
   long * pos = *posp;
   /*
   for (int i = 0; i < NUM_BASES; i++) fprintf(stdout, " %ld", C[i]);
   fprintf(stdout, "\nPOS[] =");
   for (int i = 0; i < *gsize; i++) fprintf(stdout, " %ld", pos[i]);
   fprintf(stdout, "\n");
   for (int i = 0; i < NUM_BASES; i++) {
      fprintf(stdout, "occ[%d] =", i);
      for (int j = 0; j < occ[i].max; j++) fprintf(stdout, " %ld", occ[i].val[j]);
      fprintf(stdout, "\n");
   }
   */
   return index;
}

void
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
   char * genome = compact_genome(filename, &gsize);
   C = compute_c(genome, gsize);

   // Output files.
   char * fmfile = malloc(strlen(filename)+5);
   strcpy(fmfile, filename);
   strcpy(fmfile + strlen(filename), ".fmi");
   int fd = open(fmfile, O_WRONLY | O_CREAT | O_TRUNC, 0644);
   free(fmfile);

   // FM index
   bwt_index(genome, gsize, &pos, occ);

   // Write index.
   write(fd, C, NUM_BASES*sizeof(long));
   write(fd, &gsize, sizeof(long));
   write(fd, pos, gsize*sizeof(long));
   for (int i = 0; i < NUM_BASES; i++) {
      write(fd, &(occ[i]->pos), sizeof(long));
      write(fd, &(occ[i]->val[0]), occ[i]->pos * sizeof(long));
   }

   /*
   fprintf(stdout, "BWT(%s) = %s\nC[] =", genome, btext);
   for (int i = 0; i < NUM_BASES; i++) fprintf(stdout, " %ld", C[i]);
   fprintf(stdout, "\nPOS[] =");
   for (int i = 0; i < strlen(genome); i++) fprintf(stdout, " %ld", pos[i]);
   fprintf(stdout, "\n");
   for (int i = 0; i < NUM_BASES; i++) {
      fprintf(stdout, "occ[%d] =", i);
      for (int j = 0; j < occ[i]->pos; j++) fprintf(stdout, " %ld", occ[i]->val[j]);
      fprintf(stdout, "\n");
   }
   */
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
         fprintf(output,"%s %ld\n",buffer+1,*genomesize);
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


long *
compute_c
(
 char * genome,
 long gsize
)
{
   long * cnt = calloc(NUM_BASES, sizeof(long));
   for (int i = 0; i < gsize; i++) cnt[translate[genome[i]]]++;
   long * c = malloc(NUM_BASES * sizeof(long));
   c[0] = 0;
   for (int i = 1; i < NUM_BASES; i++) c[i] = c[i-1] + cnt[i-1];
   free(cnt);
   return c;
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
   fprintf(stdout,"Suffix sort...");
   long * sa = dc3(genome);
   fprintf(stdout,"done.\n");

   // Fill compacted occ.
   vstack_t * stacks[NUM_BASES];
   for (int i = 0; i < NUM_BASES; i++) stacks[i] = new_stack(gsize/NUM_BASES);

   for (long i = 0; i < gsize; i++) push(stacks+translate[genome[(sa[i] - 1 < 0 ? gsize-1 : sa[i] - 1)]], i);
   
   // Save stacks to output.
   *pos = sa;

   for (int i = 0; i < NUM_BASES; i++) {
      stacks[i] = realloc(stacks[i], sizeof(vstack_t) + stacks[i]->pos*sizeof(long));
      if (stacks[i] == NULL) {
         fprintf(stderr, "error in bwt_index (realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      stacks[i]->size = stacks[i]->pos;
      occ[i] = stacks[i];
   }
}





/*
void
prefix_sort
(
 char     *  text,
 long        tlen,
 long     ** cp,
 vstack_t ** sorted
)
{
   jstack_t * jstack = new_jstack(128);
   vstack_t * stacks[NUM_BASES];

   for (int i = 0; i < NUM_BASES; i++) {
      stacks[i] = new_stack(tlen/NUM_BASES);
      stacks[i]->offset = 1;
   }
   
   
   // Classify prefixes
   for (long i = 0; i < job->pos; i++) push(stacks+translate[text[(job->val[i]+job->offset)%tlen]], job->val[i]);
   for (int i = 0; i < NUM_BASES; i++) stacks[i] = realloc(stacks[i], sizeof(vstack_t) + stacks[i]->pos*sizeof(long));

   // Compute C[] from stacks[i]->pos.
   long * c = malloc(NUM_BASES*sizeof(long));
   *cp = c;
   c[0] = 0;
   for (int i = 1; i < NUM_BASES; i++) c[i] = c[i-1] + end[i-1];

   // Commit jobs.
   for (long i = NUM_BASES-1; i >=0; i--) 
      if (end[i] - start[i] > 0)
         push_job(&jstack, (void *)new_job(start[i], end[i], 0));
  
   // Recursive call.
   job_t * job;
   while ((job = (job_t *)pop_job(jstack)) != NULL) {
      if (job->end - job->start == 1) {
         push(sorted, stack[job->buffer]->val[0]);
         free(job);
         continue;
      }
      // Job data.
      int src  = job->buffer;
      int dest = (job->buffer+1)%2;
      // Classify prefixes
      k = job->start;
      for (int n = 0; n < NUM_BASES; n++) {
         start[n] = k;
         for (long i = job->start; i < job->end; i++) {
            if (translate[text[(stacks[src]->val[i] + job->offset)%tlen]] == n) stacks[dest]->val[k++] = stacks[src]->val[i];
         }
         end[n] = k;
      }

      for (long i = 0; i < job->pos; i++) push(stacks+translate[text[(job->val[i]+job->offset)%tlen]], job->val[i]);
      for (int i = 0; i < NUM_BASES; i++) stacks[i] = realloc(stacks[i], sizeof(vstack_t) + stacks[i]->pos*sizeof(long));
      // Free job.
      free(job);
      // Commit new jobs.
      for (long i = NUM_BASES-1; i >=0; i--) {
         if (stacks[i]->pos > 0)
            push_job(&jstack, (void *)stacks[i]);
         else {
            free(stacks[i]);
         }
      }
   }
}
*/

jstack_t *
new_jstack
(
 long size
)
{
   jstack_t * stack = malloc(sizeof(jstack_t) + size * sizeof(void *));
   if (stack == NULL) {
      fprintf(stderr, "error in 'new_jstack' (malloc): %s\n", strerror(errno));
      exit(EXIT_FAILURE);
   }
   stack->size = size;
   stack->pos = 0;
   return stack;
}

void
push_job
(
 jstack_t ** stackp,
 void     *  job
)
{
   jstack_t * stack = *stackp;
   if (stack->pos >= stack->size) {
      long newsize = stack->size * 2;
      *stackp = stack = realloc(stack, sizeof(jstack_t) + newsize * sizeof(void *));
      if (stack == NULL) {
         fprintf(stderr, "error in 'push_job' (realloc): %s\n", strerror(errno));
         exit(EXIT_FAILURE);
      }
      stack->size = newsize;
   }
   stack->job[stack->pos++] = job;
}

void *
pop_job
(
 jstack_t * stack
 )
{
   if (stack->pos - 1 >= 0) return stack->job[--stack->pos];
   else return NULL;
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
