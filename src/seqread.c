#include "seqread.h"

// Define data structures.

struct seqread_t {
   size_t    seqlen;
   char    * tag;
   char    * seq;
   char    * qscore;
};

// Interface functions.

seqread_t *
seqread_new
(
 char  * tag,
 char  * seq,
 char  * qscore
)
{
   // Declare variables.
   seqread_t * read = NULL;

   // Check arguments.
   error_test_msg(seq == NULL, "argument 'seq' is NULL.");
   
   // Alloc structure.
   read = malloc(sizeof(seqread_t));
   error_test_mem(read);

   read->tag    = NULL;
   read->seq    = NULL;
   read->qscore = NULL;
   read->seqlen = strlen(seq);
   
   // Duplicate strings
   read->seq = strdup(seq);
   error_test_def(read->seq == NULL);

   if (tag != NULL) {
      read->tag = strdup(tag);
      error_test_def(read->tag == NULL);
   }

   if (qscore != NULL) {
      if (strlen(qscore) != read->seqlen) {
         error_throw_msg("inconsistent length of sequence quality.");
      }
      read->qscore = strdup(qscore);
      error_test_def(read->qscore == NULL);
   }

   return read;

 failure_return:
   seqread_free(read);
   return NULL;
}

gstack_t *
seqread_stack
(
  size_t max_elm
)
{
   gstack_t * stack = gstack_new(max_elm, seqread_free);
   error_test(stack == NULL);
   
   return stack;

 failure_return:
   return NULL;
}

seqread_t *
seqread_pop
(
  gstack_t * stack
)
{
   seqread_t * read = (seqread_t *) gstack_pop(stack);
   return read;
}

seqread_t *
seqread_get
(
  size_t     index,
  gstack_t * stack
)
{
   seqread_t * read = (seqread_t *) gstack_get(index, stack);
   error_test(read == NULL);

   return read;
   
 failure_return:
   return NULL;
}

void
seqread_free
(
  void * ptr
)
{
   seqread_t * read = (seqread_t *) ptr;
   if (read != NULL) {
      free(read->tag);
      free(read->seq);
      free(read->qscore);
      free(read);
   }

   return;
}

seqread_t *
seqread_parse_fastq
(
  char ** lines
)
{
   // Declare variables.
   seqread_t * read = NULL;

   // Check arguments.
   error_test_msg(lines == NULL, "argument 'lines' is NULL.");

   // Alloc seqread_t.
   read = seqread_new(lines[0], lines[1], lines[3]);
   error_test(read);

   return read;
   
 failure_return:
   seqread_free(read);
   return NULL;
}

seqread_t *
seqread_parse_fasta
(
  char ** lines
)
{
   
}

seqread_t *
seqread_parse_raw
(
  char ** lines
)
{
   
}

// Helper functions.

int64_t
seqread_len
(
  seqread_t * read
)
{
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   
   return (int64_t)read->seqlen;
   
 failure_return:
   return -1;
}

char *
seqread_seq
(
  seqread_t * read
)
{
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   
   return read->seq;
   
 failure_return:
   return NULL;
}

char *
seqread_tag
(
  seqread_t * read
)
{
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   
   return read->tag;
   
 failure_return:
   return NULL;
}

char *
seqread_qscore
(
  seqread_t * read
)
{
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   
   return read->qscore;
   
 failure_return:
   return NULL;
}
