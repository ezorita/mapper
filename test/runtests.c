#include "unittest.h"

int
main(
   int argc,
   char **argv
)
{

   // Import the test cases from linked files. //

   // Unit tests.
   extern test_case_t test_cases_index_sym[];
   extern test_case_t test_cases_index_txt[];
   extern test_case_t test_cases_index_sar[];
   extern test_case_t test_cases_index_bwt[];
   extern test_case_t test_cases_index_ann[];
   extern test_case_t test_cases_index[];
   extern test_case_t test_cases_ui[];
   extern test_case_t test_cases_blocksearch[];
   extern test_case_t test_cases_gstack[];
   extern test_case_t test_cases_seqread[];
   extern test_case_t test_cases_io[];
   extern test_case_t test_cases_mapper[];

   // Memory tests.
   extern test_case_t test_mem_index_sym[];
   extern test_case_t test_mem_index_txt[];
   extern test_case_t test_mem_index_sar[];
   extern test_case_t test_mem_index_bwt[];
   extern test_case_t test_mem_index_ann[];
   extern test_case_t test_mem_index[];
   extern test_case_t test_mem_ui[];
   extern test_case_t test_mem_blocksearch[];
   extern test_case_t test_mem_gstack[];
   extern test_case_t test_mem_seqread[];
   extern test_case_t test_mem_io[];

   const test_case_t *test_case_list[] = {
      test_cases_index_sym,
      test_cases_index_txt,
      test_cases_index_sar,
      test_cases_index_bwt,
      test_cases_index_ann,
      test_cases_index,
      test_cases_ui,
      test_cases_blocksearch,
      test_cases_gstack,
      test_cases_seqread,
      test_cases_io,
      test_cases_mapper,
      test_mem_index_sym,
      test_mem_index_txt,
      test_mem_index_sar,
      test_mem_index_bwt,
      test_mem_index_ann,
      test_mem_index,
      test_mem_ui,
      test_mem_blocksearch,
      test_mem_gstack,
      test_mem_seqread,
      test_mem_io,
      NULL, // Sentinel. //
   };

   return run_unittest(argc, argv, test_case_list);

}
