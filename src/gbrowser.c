#include "algs.h"
#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <unistd.h>

#define CHRSTR_SIZE 100
#define BUFFER_SIZE 100

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


int main(int argc, char *argv[])
{
   if (argc != 5) exit(EXIT_FAILURE);
   
   char * chromosome = argv[2];   
   long locus = atol(argv[3]);
   int len = atoi(argv[4]);

   char * filename = malloc(strlen(argv[1])+5);
   strcpy(filename, argv[1]);

   // mmap index file.
   strcpy(filename+strlen(argv[1]), ".fmi");
   int fd = open(filename, O_RDONLY);
   if (fd == -1) {
      fprintf(stderr, "error opening index\n");
      exit(EXIT_FAILURE);
   }

   int mflags = MAP_PRIVATE;
   long idxsize = lseek(fd, 0, SEEK_END);
   lseek(fd, 0, SEEK_SET);
   long * indexp = mmap(NULL, idxsize, PROT_READ, mflags, fd, 0);
   if (indexp == NULL) {
      fprintf(stderr, "error opening index file\n");
      return EXIT_FAILURE;
   }
      
   // Read FM index format.
   index_t index;
   if(format_FMindex(indexp, &index))
      return EXIT_FAILURE;

   // Read chromosome index.
   strcpy(filename+strlen(argv[1]), ".index");
   chr_t * chr = read_CHRindex(filename);
   if (chr == NULL) return EXIT_FAILURE;

   for (int i = 0; i < chr->nchr; i++) {
      if (strcmp(chromosome, chr->name[i]) == 0) {
         long g_start = index.gsize - chr->start[i] - locus;
         fprintf(stdout, "%s:%ld-%ld\n", chr->name[i], locus, locus+len);
         for (int i = 0; i < len; i++) fprintf(stdout, "%c", index.genome[g_start-i]);
         fprintf(stdout, "\n");
      }

   }


   return 0;
}
