#include "unittest.h"
// Include target source file to test in order to
// define all the declared symbols.
#include "indexbuild.h"

// Write the test cases as function with prototype 'void (*) (void)'.
void
test_compute_sa
(void)
{
   char   * seq_1   = "TAGCNTCGACA$";
   int64_t  sar_1[] = {11,10,8,1,9,6,3,7,2,0,5,4};
   char   * seq_2   = "TGNGA$TCGAT$";
   int64_t  sar_2[] = {11,5,4,9,7,3,8,1,10,6,0,2};
   char   * seq_3   = "AAaaAAaaA$";
   int64_t  sar_3[] = {9,8,7,6,5,4,3,2,1,0};
   char   * seq_4   = "acgTPGaKcnN$";
   int64_t  sar_4[] = {11,0,6,1,8,5,2,3,10,7,4,9};
   char   * seq_5   = "$TgAgcatGGC$acTCGA$";
   int64_t  sar_5[] = {18,11,0,17,12,3,6,10,5,15,13,16,2,9,4,8,14,1,7};

   int64_t * sa;
   sa = compute_sa(seq_1, strlen(seq_1));
   test_assert_critical(sa != NULL);
   test_assert(memcmp(sa, sar_1, strlen(seq_1)*sizeof(int64_t)) == 0);
   free(sa);
   sa = compute_sa(seq_2, strlen(seq_2));
   test_assert_critical(sa != NULL);
   test_assert(memcmp(sa, sar_2, strlen(seq_2)*sizeof(int64_t)) == 0);
   free(sa);
   sa = compute_sa(seq_3, strlen(seq_3));
   test_assert_critical(sa != NULL);
   test_assert(memcmp(sa, sar_3, strlen(seq_3)*sizeof(int64_t)) == 0);
   free(sa);
   sa = compute_sa(seq_4, strlen(seq_4));
   test_assert_critical(sa != NULL);
   test_assert(memcmp(sa, sar_4, strlen(seq_4)*sizeof(int64_t)) == 0);
   free(sa);
   sa = compute_sa(seq_5, strlen(seq_5));
   test_assert_critical(sa != NULL);
   test_assert(memcmp(sa, sar_5, strlen(seq_5)*sizeof(int64_t)) == 0);
   free(sa);
}

void
test_compute_occ
(void)
{
   char   * seq_1   = "TAGCNTCGACA$";
   int64_t  sar_1[] = {11,10,8,1,9,6,3,7,2,0,5,4};
   char   * seq_2   = "TGNGA$TCGAT$";
   int64_t  sar_2[] = {11,5,4,9,7,3,8,1,10,6,0,2};
   char   * seq_3   = "$TgAgcatGGC$acTCGA$";
   int64_t  sar_3[] = {18,11,0,17,12,3,6,10,5,15,13,16,2,9,4,8,14,1,7};
   
   
}

/*
void
test_case_1
(void)
{
   // Use 'test_assert()' to assert results are as expected.
   test_assert(return_0() == 0);
   test_assert(return_1() == 1);
}


void
test_case_2
(void)
{
   test_assert(should_return_0() == 0);
   test_assert(return_0() == 0);
   test_assert(return_1() == 1);
}


void
test_case_3
(void)
{
   redirect_stderr();
   errmsg();
   unredirect_stderr();
   test_assert_stderr("error message");
}


void
test_case_4
(void)
{
   // Use 'redirect_stderr()' to record stderr.
   redirect_stderr();
   errmsg();
   // Use 'unredirect_stderr()' to stop recording.
   unredirect_stderr();
   // Use 'test_assert_stderr()' to compare the recorded
   // error message with the target.
   test_assert_stderr("wrong error message");
   // The recorded text is in 'caught_in_stderr()' .
   fprintf(stderr, "%s", caught_in_stderr());
}


void
test_case_5
(void)
{
   // In this test a segmentation fault will happen.
   // This will cause the test to fail, but the execution
   // will nevertheless continue to the next test.

   if (*(int*) 0x0 == 0) {
      printf("won't be executed (segfault before)\n");
   }

   return;

}


void
test_case_6
(void)
{
   char * c = call_malloc(); 
   // Use 'test_assert_critical()' to interrupt the test
   // in case assertion fails. In the example below 'c'
   // must not be NULL otherwise 'free(c)' will fail.
   test_assert_critical(c != NULL);
   free(c);

   // Use 'set_alloc_failure_rate_to(p)' to cause 'malloc()'
   // (including 'calloc()' and 'realloc()') to fail randomly
   // with probability 'p'.
   set_alloc_failure_rate_to(1.0);
   c = call_malloc();
   test_assert(c == NULL);
   // Set 'malloc()' to normal mode.
   reset_alloc();

   // Use 'set_alloc_failure_countdown_to(n)' to cause 'malloc()'
   // (including 'calloc()' and 'realloc()') to systematically
   // fail after 'n' calls.
   set_alloc_failure_countdown_to(1);
   c = call_malloc();
   test_assert_critical(c != NULL);
   free(c);
   // The second call to 'malloc()' will fail.
   c = call_malloc();
   test_assert(c == NULL);
   reset_alloc();

}
*/

// Define test cases to be run (for export).
const test_case_t test_cases_from_file_1[] = {
   {"indexbuild/compute_sa",  test_compute_sa},
   {NULL, NULL}, // Sentinel. //
};
