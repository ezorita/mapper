#include "seqread.h"

// Define data structures.

struct seqread_t {
   size_t    seqlen;
   char    * tag;
   char    * seq;
   char    * qscore;
   uint8_t * sym;
};

// Interface functions.

seqread_t *
seqread_new
(
 char  * tag,
 char  * seq,
 char  * qscore,
 sym_t * sym
)
{
   // Declare variables.
   seqread_t * read = NULL;

   // Check arguments.
   error_test_msg(seq == NULL, "argument 'seq' is NULL.");
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");
   
   // Alloc structure.
   read = malloc(sizeof(seqread_t));
   error_test_mem(read);

   read->tag    = NULL;
   read->seq    = NULL;
   read->qscore = NULL;
   read->sym    = NULL;
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

   // Compute symbols.
   read->sym = sym_str_index(seq, sym);
   error_test(read->sym == NULL);

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
      free(read->sym);
      free(read);
   }

   return;
}

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

uint8_t *
seqread_sym
(
  seqread_t * read
)
{
   error_test_msg(read == NULL, "argument 'read' is NULL.");
   
   return read->sym;
   
 failure_return:
   return NULL;
}
