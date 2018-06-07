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
   test_assert(mapper("examples/seed_query.fa", index) == 0);
   unredirect_stdout();
   unredirect_stderr();
   test_assert(strcmp("tag:q1, mems:2\ntag:q2, mems:2\ntag:q3, mems:1\ntag:q4, mems:5\n", caught_in_stdout()) == 0);

   index_free(index);
   
}

// Define test cases to be run (for export).
const test_case_t test_cases_mapper[] = {
   {"mapper/mapper",         test_mapper},
   {NULL, NULL}, // Sentinel. //
};


