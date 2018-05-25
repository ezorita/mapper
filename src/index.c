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
   // Declare variables.
   char * file_path = NULL;
   index_t * index = NULL;
   glob_t gbuf;
   int glob_flag = 0;

   // Check arguments.
   error_test_msg(filename_base == NULL, "argument 'filename_base' is NULL.");

   // Alloc index structure.
   index = malloc(sizeof(index_t));
   error_test_mem(index);

   index->fname_base = strdup(filename_base);
   error_test_mem(index->fname_base);
   index->ann_cnt = 0;
   index->sym = NULL;
   index->txt = NULL;
   index->sar = NULL;
   index->bwt = NULL;
   index->ann = NULL;
   
   file_path = malloc(strlen(filename_base)+7);
   error_test_mem(file_path);

   // sym
   sprintf(file_path, "%s.sym", filename_base);
   index->sym = sym_file_read(file_path);
   error_test(index->sym == NULL);

   // txt
   sprintf(file_path, "%s.txt", filename_base);
   index->txt = txt_file_read(file_path, index->sym);
   error_test(index->txt == NULL);

   // sar
   sprintf(file_path, "%s.sar", filename_base);
   index->sar = sar_file_read(file_path);
   error_test(index->sar == NULL);

   // bwt
   sprintf(file_path, "%s.bwt", filename_base);
   index->bwt = bwt_file_read(file_path, index->txt);
   error_test(index->bwt == NULL);

   // ann glob
   sprintf(file_path, "%s.ann.*", filename_base);   
   glob(file_path, 0, NULL, &gbuf);
   glob_flag = 1;

   // alloc ann list
   if (gbuf.gl_pathc > 0) {
      index->ann = malloc(gbuf.gl_pathc * sizeof(void *));
      error_test(index->ann == NULL);
   }

   // ann
   index->ann_cnt = gbuf.gl_pathc;
   for (int i = 0; i < index->ann_cnt; i++) {
      index->ann[i] = ann_file_read(gbuf.gl_pathv[i]);
      error_test(index->ann[i] == NULL);
   }

   globfree(&gbuf);
   free(file_path);
   return index;

 failure_return:
   if (glob_flag)
      globfree(&gbuf);
   index_free(index);
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
   char * file_path = NULL;
   index_t * index = NULL;

   error_test_msg(ref_txt == NULL, "argument 'ref_txt' is NULL.");
   error_test_msg(filename_base == NULL, "argument 'filename_base' is NULL.");

   index = malloc(sizeof(index_t));
   error_test_mem(index);

   index->sym = NULL;
   index->txt = NULL;
   index->sar = NULL;
   index->bwt = NULL;
   index->ann = NULL;
   index->ann_cnt = 0;

   index->fname_base = strdup(filename_base);
   error_test_mem(index->fname_base);

   // Compute index structures.
   index->txt = read_fasta(ref_txt);
   error_test(index->txt == NULL);

   index->sym = txt_get_symbols(index->txt);
   error_test(index->sym == NULL);

   index->sar = sar_build(index->txt);
   error_test(index->sar == NULL);

   index->bwt = bwt_build(index->txt, index->sar);
   error_test(index->bwt == NULL);

   // Write to output files.
   file_path = malloc(strlen(filename_base)+5);
   error_test_mem(file_path);

   sprintf(file_path, "%s.sym", filename_base);
   error_test(sym_file_write(file_path, index->sym) == -1);

   sprintf(file_path, "%s.txt", filename_base);
   error_test(txt_file_write(file_path, index->txt) == -1);

   sprintf(file_path, "%s.sar", filename_base);
   error_test(sar_file_write(file_path, index->sar) == -1);

   sprintf(file_path, "%s.bwt", filename_base);
   error_test(bwt_file_write(file_path, index->bwt) == -1);

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
   // Declare variables.
   char * file_path = NULL;

   // Check arguments.
   error_test_msg(index == NULL, "argument 'index' is NULL.");
   error_test_msg(kmer < 1, "argument 'kmer' must be greater than 0.");
   error_test_msg(threads < 1, "argument 'threads' must be greater than 0.");
   error_test_msg(tau < 1, "argument 'tau' must be greater than 0");
   
   // Check whether annotation exists.
   for (int i = 0; i < index->ann_cnt; i++) {
      if (ann_get_kmer(index->ann[i]) == kmer && ann_get_dist(index->ann[i]) == tau)
         return 1;
   }

   // Realloc annotation list.
   ann_t ** new_ann = realloc(index->ann, (index->ann_cnt+1) * sizeof(void *));
   error_test_mem(new_ann);
   index->ann = new_ann;

   // Compute annotation.
   index->ann[index->ann_cnt] = ann_build(kmer, tau, index->bwt, index->sar, threads);
   error_test(index->ann[index->ann_cnt] == NULL);
   
   // Write to file.
   file_path = malloc(strlen(index->fname_base) + 100);
   error_test_mem(file_path);
   sprintf(file_path, "%s.ann.%d.%d", index->fname_base, kmer, tau);
   error_test(ann_file_write(file_path, index->ann[index->ann_cnt]) == -1);

   index->ann_cnt++;

   free(file_path);

   return 0;

 failure_return:
   free(file_path);
   return -1;
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
   // Declare variables.
   FILE * input   = NULL;
   char * seqname = NULL;
   char * buffer  = NULL;
   txt_t * txt    = NULL;
   sym_t * sym    = NULL;

   // Check arguments.
   error_test_msg(filename == NULL, "argument 'filename' is NULL.");

   // Files
   input = fopen(filename,"r");
   error_test_def(input == NULL);

   // Check FASTA format
   error_test_msg(fgetc(input) != '>', "incorrect input format (fgetc).");
   
   ungetc('>', input);

   // Alloc txt structure.
   sym = sym_new_dna();
   error_test(sym == NULL);
   txt = txt_new(sym);
   error_test(txt == NULL);

   // File read vars
   buffer  = malloc(INDEX_FASTA_BUFFER);
   error_test_mem(buffer);

   ssize_t   rlen;
   size_t    sz     = INDEX_FASTA_BUFFER;
   uint64_t  lineno = 0;
   while ((rlen = getline(&buffer, &sz, input)) != -1) {
      lineno++;
      if (buffer[0] == '>') {
         // Commit previous sequence
         if (seqname != NULL) {
            error_test(txt_commit_seq(seqname, txt) == -1);
            free(seqname);
            seqname = NULL;
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
         error_test_msg(beg == k, "found empty sequence name.");
         // Store sequence name.
         seqname = strdup(buffer+beg);
         error_test_mem(seqname);
      } else {
         // Remove newline character.
         if (buffer[rlen-1] == '\n') {
            buffer[rlen-1] = 0;
         }
         // Append sequence.
         error_test(txt_append(buffer, txt) == -1);
      }
   }

   // Commit last sequence.
   error_test(txt_commit_seq(seqname, txt) == -1);
   // Generate reverse complement.
   error_test(txt_commit_rc(txt) == -1);

   free(buffer);
   free(seqname);
   fclose(input);
   
   return txt;
   
 failure_return:
   if (input != NULL)
      fclose(input);
   free(buffer);
   free(seqname);
   txt_free(txt);
   sym_free(sym);
   return NULL;
}
