#include "unittest.h"
#include "io.h"

void
test_mem_io_stream_open
(void)
{
   iostream_t * stream = NULL;
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      stream = io_stream_open("examples/io_input.raw");
      io_stream_free(stream);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      stream = io_stream_open("examples/io_input.raw");
      io_stream_free(stream);
   }
   reset_alloc();
   unredirect_stderr();
}

void
test_mem_io_stream_read_seq
(void)
{
   iostream_t * stream;
   gstack_t   * stack;

   // Invalid arguments.
   redirect_stderr();
   test_assert(io_stream_read_seq(NULL) == NULL);
   unredirect_stderr();

   // Load files.
   // Read RAW file.
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   for (int i = 0; i < 1000; i++) {
      stream = io_stream_open("examples/io_input.fastq");
      test_assert_critical(stream != NULL);

      set_alloc_failure_rate_to(0.1);
      stack = io_stream_read_seq(stream);
      reset_alloc();
      
      gstack_free(stack);
      io_stream_free(stream);
   }
   reset_alloc();
   
   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      stream = io_stream_open("examples/io_input.fasta");
      test_assert_critical(stream != NULL);

      set_alloc_failure_countdown_to(i);
      stack = io_stream_read_seq(stream);
      reset_alloc();
      
      gstack_free(stack);
      io_stream_free(stream);
   }
   reset_alloc();
   
}

// Define test cases to be run (for export).
const test_case_t test_mem_io[] = {
   {"mem/io/io_stream_open",         test_mem_io_stream_open},
   {"mem/io/io_stream_read_seq",     test_mem_io_stream_read_seq},
   {NULL, NULL}, // Sentinel. //
};
