#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdint.h>


int main(int argc, char *argv[])
{
   size_t hist[256];
   for (int i = 0; i < 256; i++) hist[i]=0;
   int fd = open(argv[1],O_RDONLY);
   size_t fsize = lseek(fd, 0, SEEK_END);
   uint8_t * repeats = mmap(NULL,fsize,PROT_READ,MAP_PRIVATE,fd,0);
   size_t uniq = 0;
   for (size_t i = 0; i < fsize/2; i++) {
      uniq += repeats[i] < 3;
      hist[repeats[i]]++;
   }
   fprintf(stderr,"%ld out of %ld (%.2f)\n",uniq,fsize/2, uniq*2.0/fsize);
   for (int i = 0; i < 256; i++) fprintf(stdout, "%ld\n",hist[i]);
   return 0;
}
