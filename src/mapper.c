#include "mapper.h"

int
mapper
(
 char    * qfile,
 index_t * index
)
{
   // Declare variables.
   iostream_t * qinput = NULL;
   gstack_t   * rstack = NULL;

   // Check arguments.
   error_test_msg(qfile == NULL, "argument 'qfile' is NULL.");
   error_test_msg(index == NULL, "argument 'index' is NULL.");

   // Open input stream.
   qinput = io_stream_open(qfile);
   error_test(qinput == NULL);

   // Read file by chunks and report sequences in output stream.
   while ((rstack = io_stream_read_seq(qinput)) != NULL) {
      // Break when no more sequences are found.
      if (gstack_num_elm(rstack) == 0) {
	 break;
      }
      
      // Stream info.
      fprintf(stderr, "read: %ld sequences.\n", gstack_num_elm(rstack));
      
      // Iterate over sequences.
      for (int i = 0; i < gstack_num_elm(rstack); i++) {
	 seqread_t * read = seqread_get(i, rstack);
	 fprintf(stdout, "%s\n", seqread_seq(read));
      }
      
      // Free stack.
      gstack_free(rstack);
   }

   io_stream_free(qinput);
   return 0;

 failure_return:
   return -1;
}
