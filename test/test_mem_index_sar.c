#include "unittest.h"
#include "index_sar.h"

void
test_mem_sar_build
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   sar_t * sar;

   test_assert(txt_append("TAGCNTCGACA", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sar = sar_build(txt);
      sar_free(sar);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sar = sar_build(txt);
      sar_free(sar);
   }
   reset_alloc();

   sym_free(sym);
   txt_free(txt);
   unredirect_stderr();
}

void
test_mem_sar_get
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append("TAGCNTCGACA", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sar_get(0,sar);
      sar_get(10,sar);
      sar_get(19,sar);
      sar_get(-1,sar);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sar_get(0,sar);
      sar_get(10,sar);
      sar_get(19,sar);
      sar_get(-1,sar);
   }
   reset_alloc();

   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_sar_get_range
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append("TAGCNTCGACA", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   int64_t * sa_buf = malloc(30*sizeof(int64_t));
   test_assert_critical(sa_buf != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sar_get_range(0,10,sa_buf,sar);
      sar_get_range(10,5,sa_buf,sar);
      sar_get_range(19,30,sa_buf,sar);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sar_get_range(0,10,sa_buf,sar);
      sar_get_range(10,5,sa_buf,sar);
      sar_get_range(19,30,sa_buf,sar);
   }
   reset_alloc();

   free(sa_buf);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_sar_file
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append("TAGCNTCGACA", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   sar_t * sar_i;

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sar_file_write("test30.sar", sar);
      sar_i = sar_file_read("test30.sar");
      sar_free(sar_i);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sar_file_write("test30.sar", sar);
      sar_i = sar_file_read("test30.sar");
      sar_free(sar_i);
   }
   reset_alloc();

   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_mem_index_sar[] = {
   {"mem/index_sar/sar_build",             test_mem_sar_build},
   {"mem/index_sar/sar_get",               test_mem_sar_get},
   {"mem/index_sar/sar_get_range",         test_mem_sar_get_range},
   {"mem/index_sar/sar_file",              test_mem_sar_file},
   {NULL, NULL}, // Sentinel. //
};

