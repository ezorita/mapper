#include "unittest.h"
#include "io.h"


void
test_io_stream_open
(void)
{
   iostream_t * stream;

   // Test invalid arguments.
   redirect_stderr();
   stream = io_stream_open("fake-file.txt");
   test_assert(stream == NULL);
   unredirect_stderr();

   // Load files.
   stream = io_stream_open("examples/io_input.raw");
   test_assert(stream != NULL);
   io_stream_free(stream);
 
   stream = io_stream_open("examples/io_input.fasta");
   test_assert(stream != NULL);
   io_stream_free(stream);

   stream = io_stream_open("examples/io_input.fastq");
   test_assert(stream != NULL);
   io_stream_free(stream);
}

void
test_io_stream_open_buf
(void)
{
   iostream_t * stream;

   // Test invalid arguments.
   redirect_stderr();
   stream = io_stream_open_buf("fake-file.txt", 100);
   test_assert(stream == NULL);
   unredirect_stderr();

   // Load files.
   stream = io_stream_open_buf("examples/io_input.raw", 100);
   test_assert(stream != NULL);
   io_stream_free(stream);
 
   stream = io_stream_open_buf("examples/io_input.fasta", 100);
   test_assert(stream != NULL);
   io_stream_free(stream);

   stream = io_stream_open_buf("examples/io_input.fastq", 100);
   test_assert(stream != NULL);
   io_stream_free(stream);   
}

void
test_io_stream_read_seq
(void)
{
   iostream_t * stream;
   gstack_t   * stack;
   seqread_t  * read;

   // Invalid arguments.
   redirect_stderr();
   test_assert(io_stream_read_seq(NULL) == NULL);
   unredirect_stderr();

   // Load files.
   // Read RAW file.
   stream = io_stream_open("examples/io_input.raw");
   test_assert_critical(stream != NULL);
   
   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 6);
   
   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "NNNNNNLLYYYYJDFLS") == 0);
   test_assert(strcmp(seqread_tag(read), "6") == 0);
   seqread_free(read);
   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "NNNNNNNNACGTACGCC") == 0);
   test_assert(strcmp(seqread_tag(read), "5") == 0);
   seqread_free(read);
   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "ACGTACATGTATGACAC") == 0);
   test_assert(strcmp(seqread_tag(read), "4") == 0);
   seqread_free(read);

   gstack_free(stack);
   io_stream_free(stream);

   // Read FASTA file.
   stream = io_stream_open_buf("examples/io_input.fasta", 80);
   test_assert_critical(stream != NULL);

   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 4);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "ACGTACATGTATGACAC") == 0);
   test_assert(strcmp(seqread_tag(read), "seq4") == 0);
   seqread_free(read);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "AGTCGANTATACNTACG") == 0);
   test_assert(strcmp(seqread_tag(read), "seq3") == 0);
   seqread_free(read);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "GTATCGACTACGAGCTA") == 0);
   test_assert(strcmp(seqread_tag(read), "seq2") == 0);
   seqread_free(read);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "ATGCGTACGTCGTATCA") == 0);
   test_assert(strcmp(seqread_tag(read), "seq1") == 0);
   seqread_free(read);
   gstack_free(stack);
   
   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 2);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "NNNNNNLLYYYYJDFLS") == 0);
   test_assert(strcmp(seqread_tag(read), "seq6") == 0);
   seqread_free(read);
   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "NNNNNNNNACGTACGCC") == 0);
   test_assert(strcmp(seqread_tag(read), "seq5") == 0);
   seqread_free(read);


   gstack_free(stack);
   io_stream_free(stream);

   // Read FASTQ file.
   stream = io_stream_open_buf("examples/io_input.fastq", 50);
   test_assert_critical(stream != NULL);
   
   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 2);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "GTATCGACTACGAGCTA") == 0);
   test_assert(strcmp(seqread_tag(read), "seq2") == 0);
   test_assert(strcmp(seqread_qscore(read), "BACBABCA0ACBB00AC") == 0);   
   seqread_free(read);

   read = seqread_pop(stack);
   test_assert(read != NULL);
   test_assert(strcmp(seqread_seq(read), "ATGCGTACGTCGTATCA") == 0);
   test_assert(strcmp(seqread_tag(read), "seq1") == 0);
   test_assert(strcmp(seqread_qscore(read), "12391284194819241") == 0);   
   seqread_free(read);
   gstack_free(stack);

   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 2);
   gstack_free(stack);

   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 2);
   gstack_free(stack);

   stack = io_stream_read_seq(stream);
   test_assert_critical(stack != NULL);
   test_assert(gstack_num_elm(stack) == 0);
   gstack_free(stack);

   stack = io_stream_read_seq(stream);
   test_assert(stack == NULL);
   
   io_stream_free(stream);   
}



// Define test cases to be run (for export).
const test_case_t test_cases_io[] = {
   {"io/io_stream_open",         test_io_stream_open},
   {"io/io_stream_open_buf",     test_io_stream_open_buf},
   {"io/io_stream_read_seq",     test_io_stream_read_seq},
   {NULL, NULL}, // Sentinel. //
};
