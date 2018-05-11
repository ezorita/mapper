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
   extern test_case_t test_cases_index_sar[];
   extern test_case_t test_cases_index_bwt[];
   extern test_case_t test_cases_index_ann[];

   const test_case_t *test_case_list[] = {
      test_cases_index_sym,
      test_cases_index_txt,
      test_cases_index_sar,
      test_cases_index_bwt,
      test_cases_index_ann,
      NULL, // Sentinel. //
   };

   return run_unittest(argc, argv, test_case_list);

}
