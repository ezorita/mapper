#include "unittest.h"

int
main(
   int argc,
   char **argv
)
{

   // Import the test cases from linked files. //
   extern test_case_t test_cases_index_sym[];
   extern test_case_t test_cases_index_txt[];

   const test_case_t *test_case_list[] = {
      test_cases_index_sym,
      test_cases_index_txt,
      NULL, // Sentinel. //
   };

   return run_unittest(argc, argv, test_case_list);

}
