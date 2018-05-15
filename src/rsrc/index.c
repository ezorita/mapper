#include "index.h"

// Private function headers.
txt_t       * read_fasta       (char * filename);


index_t *
index_read
(
  char * filename_base
)
/* Reads all index structure assuming a common base name and the following extensions:
** .sym
** .txt
** .sar
** .bwt
** .ann.*
*/
{
   if (filename_base == NULL)
      return NULL;

   char * file_path = NULL;

   // Alloc index structure.
   index_t * index = malloc(sizeof(index_t));
   if (index == NULL)
      return NULL;

   index->fname_base = strdup(filename_base);
   index->ann_cnt = 0;
   index->sym = NULL;
   index->txt = NULL;
   index->sar = NULL;
   index->bwt = NULL;
   index->ann = NULL;
   
   file_path = malloc(strlen(filename_base)+7);
   if (file_path == NULL)
      goto failure_return;

   // sym
   sprintf(file_path, "%s.sym", filename_base);
   index->sym = sym_file_read(file_path);
   if (index->sym == NULL)
      goto failure_return;

   // txt
   sprintf(file_path, "%s.txt", filename_base);
   index->txt = txt_file_read(file_path, index->sym);
   if (index->txt == NULL)
      goto failure_return;

   // sar
   sprintf(file_path, "%s.sar", filename_base);
   index->sar = sar_file_read(file_path);
   if (index->sar == NULL)
      goto failure_return;

   // bwt
   sprintf(file_path, "%s.bwt", filename_base);
   index->bwt = bwt_file_read(file_path, index->txt);
   if (index->bwt == NULL)
      goto failure_return;

   // ann glob
   glob_t gbuf;
   sprintf(file_path, "%s.ann.*", filename_base);   
   glob(file_path, 0, NULL, &gbuf);

   // alloc ann list
   index->ann = malloc(gbuf.gl_pathc * sizeof(void *));
   if (index->ann == NULL)
      goto failure_return;

   // ann
   index->ann_cnt = gbuf.gl_pathc;
   for (int i = 0; i < index->ann_cnt; i++) {
      index->ann[i] = ann_file_read(gbuf.gl_pathv[i]);
      if (index->ann[i] == NULL)
         goto failure_return;
   }


   free(file_path);

   return index;
   

 failure_return:
   sym_free(index->sym);
   txt_free(index->txt);
   sar_free(index->sar);
   bwt_free(index->bwt);
   for (int i = 0; i < index->ann_cnt; i++) {
      ann_free(index->ann[i]);
   }
   free(index->ann);
   free(index);
   free(file_path);
   return NULL;
   
}


index_t *
index_build
(
  char * ref_txt,
  char * filename_base
)
{
   if (ref_txt == NULL || filename_base == NULL)
      return NULL;

   char * file_path = NULL;

   index_t * index = malloc(sizeof(index_t));
   if (index == NULL)
      goto failure_return;

   index->sym = NULL;
   index->txt = NULL;
   index->sar = NULL;
   index->bwt = NULL;
   index->ann = NULL;
   index->ann_cnt = 0;

   index->fname_base = strdup(filename_base);

   // Compute index structures.
   index->txt = read_fasta(ref_txt);
   if (index->txt == NULL)
      goto failure_return;

   index->sym = txt_get_symbols(index->txt);
   if (index->sym == NULL)
      goto failure_return;

   index->sar = sar_build(index->txt);
   if (index->sar == NULL)
      goto failure_return;

   index->bwt = bwt_build(index->txt, index->sar);
   if (index->bwt == NULL)
      goto failure_return;


   // Write to output files.
   file_path = malloc(strlen(filename_base)+5);

   sprintf(file_path, "%s.sym", filename_base);
   if (sym_file_write(file_path, index->sym))
      goto failure_return;

   sprintf(file_path, "%s.txt", filename_base);
   if (txt_file_write(file_path, index->txt))
      goto failure_return;

   sprintf(file_path, "%s.sar", filename_base);
   if (sar_file_write(file_path, index->sar))
      goto failure_return;

   sprintf(file_path, "%s.bwt", filename_base);
   if (bwt_file_write(file_path, index->bwt))
      goto failure_return;

   free(file_path);

   return index;
   
 failure_return:
   index_free(index);
   free(file_path);
   return NULL;
}


int
index_ann_new
(
 int        kmer,
 int        tau,
 int        threads,
 index_t  * index
)
{
   if (index == NULL)
      return -1;
   
   // Check whether annotation exists.
   for (int i = 0; i < index->ann_cnt; i++) {
      if (ann_get_kmer(index->ann[i]) == kmer && ann_get_dist(index->ann[i]) == tau)
         return 1;
   }

   // Realloc annotation list.
   index->ann = realloc(index->ann, (index->ann_cnt+1) * sizeof(void *));
   if (index->ann == NULL)
      return -1;

   // Compute annotation.
   index->ann[index->ann_cnt] = ann_build(kmer, tau, index->bwt, index->sar, threads);
   if (index->ann[index->ann_cnt] == NULL)
      return -1;
   
   // Write to file.
   char * file_path = malloc(strlen(index->fname_base) + 100);
   sprintf(file_path, "%s.ann.%d.%d", index->fname_base, kmer, tau);
   ann_file_write(file_path, index->ann[index->ann_cnt]);

   index->ann_cnt++;

   return 0;
}


void
index_free
(
  index_t  * index
)
{
   if (index != NULL) {
      free(index->fname_base);
      sym_free(index->sym);
      txt_free(index->txt);
      sar_free(index->sar);
      bwt_free(index->bwt);
      if (index->ann != NULL) {
         for (int i = 0; i < index->ann_cnt; i++) {
            ann_free(index->ann[i]);
         }
      }
      free(index->ann);
      free(index);
   }

   return;
}



txt_t *
read_fasta
(
  char * filename
)
{
   // Files
   FILE * input = fopen(filename,"r");
   if (input == NULL)
      return NULL;

   // Check FASTA format
   if (fgetc(input) != '>')
      return NULL;
   
   ungetc('>', input);

   // Alloc vars.
   char * seqname = NULL;
   char * buffer  = NULL;

   // Alloc txt structure.
   txt_t * txt = txt_new(sym_new_dna());
   if (txt == NULL)
      goto failure_return;

   // File read vars
   buffer  = malloc(INDEX_FASTA_BUFFER);
   if (buffer == NULL)
      goto failure_return;

   ssize_t   rlen;
   size_t    sz     = INDEX_FASTA_BUFFER;
   uint64_t  lineno = 0;
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      lineno++;
      if (buffer[0] == '>') {
         // Commit previous sequence
         if (seqname != NULL) {
            if (txt_commit_seq(seqname, txt) == -1)
               goto failure_return;
         }
         // Parse new sequence name (from '>' to first space or \0).
         buffer[rlen-1] = 0;
         // Remove leading spaces.
         int beg = 1;
         while (buffer[beg] == ' ') beg++;
         // Determine end of sequence name
         int k = beg;
         while (buffer[k] != ' ' && buffer[k] != 0) k++;
         buffer[k] = 0;
         if (beg == k)
            goto failure_return;
         // Store sequence name.
         seqname = strdup(buffer+beg);

      } else {
         // Remove newline character.
         if (buffer[rlen-1] == '\n') {
            buffer[rlen-1] = 0;
         }
         // Append sequence.
         if (txt_append(buffer, txt) == -1)
            goto failure_return;
      }
   }

   // Commit last sequence.
   txt_commit_seq(seqname, txt);
   // Generate reverse complement.
   if (txt_commit_rc(txt) == -1)
      goto failure_return;

   free(buffer);
   free(seqname);
   fclose(input);
   
   return txt;
   
 failure_return:
   fclose(input);
   free(buffer);
   free(seqname);
   txt_free(txt);
   return NULL;
}
