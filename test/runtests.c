#include "unittest.h"

int
main(
   int argc,
   char **argv
)
{

   // Import the test cases from linked files. //
   extern test_case_t test_cases_from_file_1[];

   const test_case_t *test_case_list[] = {
      test_cases_from_file_1,
      NULL, // Sentinel. //
   };

   return run_unittest(argc, argv, test_case_list);

}
