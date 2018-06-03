#include "unittest.h"
#include "mapper.h"


void
test_mapper
(void)
{
   index_t * index = index_read("ui_test00");
   test_assert_critical(index != NULL);

   redirect_stderr();
   redirect_stdout();
   test_assert(mapper("examples/io_input.fastq", index) == 0);
   unredirect_stdout();
   unredirect_stderr();
   test_assert(strcmp("ATGCGTACGTCGTATCA\nGTATCGACTACGAGCTA\nAGTCGANTATACNTACG\nACGTACATGTATGACAC\nNNNNNNNNACGTACGCC\nNNNNNNLLYYYYJDFLS\n", caught_in_stdout()));

   index_free(index);
   
}

// Define test cases to be run (for export).
const test_case_t test_cases_mapper[] = {
   {"mapper/mapper",         test_mapper},
   {NULL, NULL}, // Sentinel. //
};


