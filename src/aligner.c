#define _GNU_SOURCE
#include "align.h"

int main(int argc, char *argv[])
{
   if (argc != 3) {
      fprintf(stderr, "usage: aligner query.txt ref.txt\n");
      exit(1);
   }
   FILE * a = fopen(argv[1],"r");
   FILE * b = fopen(argv[2],"r");
   
   char * sa = NULL, * sb = NULL;
   size_t sza, szb;
   ssize_t alen, blen;

   alen = getline(&sa, &sza, a);
   blen = getline(&sb, &szb, b);
   
   if (sa[alen] == '\n') sa[alen--] = 0;
   if (sb[blen] == '\n') sb[blen--] = 0;

   // Do not count '\0' character.
   alen--;
   blen--;

   // Debug
   fprintf(stdout, "Query     [%ld nt]\n", alen);
   fprintf(stdout, "Reference [%ld nt]\n", blen);

   // Align.
   alignopt_t opt = {0.168, 26, 5, 20, 3, 0.15, 0.50, 0.05, 0.30};
   nw_align(sa, sb, alen, 1, 1, opt);
   
   return 0;
}
