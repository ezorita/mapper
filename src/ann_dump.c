#include "mapper.h"

int main(int argc, char *argv[])
{
   if (argc != 2) {
      fprintf(stderr, "usage: %s index-file\n", argv[0]);
   }

   // Load index.
   index_t * index = index_load_base(argv[1]);
   if (index == NULL) {
      fprintf(stderr, "error loading index.\n");
      exit(1);
   }

   // Load annotations.
   annlist_t * anns = ann_index_read(argv[1]);
   if (anns == NULL) {
      fprintf(stderr, "error loading annotations\n");
      exit(1);
   }

   // Dump annotations.
   for (int i = 0; i < anns->count; i++) {
      ann_t ann = anns->ann[i];
      fprintf(stderr, "Annotation info of: %s\n", ann.file);
      fprintf(stderr, "magic = 0x%llx\nkmer = %d\ntau = %d\nsize = %ld\n", ann.data->magic, ann.data->kmer, ann.data->tau, ann.data->size);
      char * buf = calloc(51,sizeof(char));
      for (int j = 0; j < ann.data->size; j += 50) {
         size_t len = min(50, ann.data->size-j-1);
         buf[len] = 0;
         strncpy(buf, index->genome+j, min(50, index->size-j-1));
         fprintf(stderr, "seq  %s\n", buf);

         fprintf(stderr, "cnt  ");
         for (int k = j; k < j+len; k++)
            fprintf(stderr, "%d", ann.data->data[k] & 0x0F);

         fprintf(stderr, "\ntau  ");
         for (int k = j; k < j+len; k++) {
            if (ann.data->data[k] & 0x0F)
               fprintf(stderr, "%d", ((ann.data->data[k] >> 4) & 0x03) + 1); 
            else
               fprintf(stderr, " ");
         }

         fprintf(stderr, "\nflag ");
         for (int k = j; k < j+len; k++) {
            if (ann.data->data[k] & 0x0F) 
               fprintf(stderr, "%d", (ann.data->data[k] >> 6) & 0x01);
            else 
               fprintf(stderr, " ");
         }

         fprintf(stderr, "\naln  ");
         for (int k = j; k < j+len; k++)
            fprintf(stderr, "%d", (ann.data->data[k] >> 7) & 0x01);
         fprintf(stderr, "\n\n");
      }
      free(buf);
   }
   return 0;
}
