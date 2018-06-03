#include "unittest.h"
#include "index_bwt.h"

void
test_mem_bwt_build
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt;

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      bwt = bwt_build(txt, sar);
      bwt_free(bwt);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      bwt = bwt_build(txt, sar);
      bwt_free(bwt);
   }
   reset_alloc();

   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_new_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      bwtquery_t * q = bwt_new_query(bwt);
      free(q);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      bwtquery_t * q = bwt_new_query(bwt);
      free(q);
   }
   reset_alloc();

   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}


void
test_mem_bwt_dup_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      bwtquery_t * qdup = bwt_dup_query(q);
      free(qdup);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      bwtquery_t * qdup = bwt_dup_query(q);
      free(qdup);
   }
   reset_alloc();

   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_new_vec
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      bwtquery_t ** qv = bwt_new_vec(bwt);
      bwt_free_vec(qv);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      bwtquery_t ** qv = bwt_new_vec(bwt);
      bwt_free_vec(qv);
   }
   reset_alloc();

   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_dup_vec
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t ** qv = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      bwtquery_t ** qvdup = bwt_dup_vec(qv);
      bwt_free_vec(qvdup);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      bwtquery_t ** qvdup = bwt_dup_vec(qv);
      bwt_free_vec(qvdup);
   }
   reset_alloc();

   bwt_free_vec(qv);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   bwtquery_t * qo = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      if (bwt_query(0, BWT_QUERY_SUFFIX, q, qo) == 0)
         bwt_query(3, BWT_QUERY_PREFIX, qo, qo);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      if (bwt_query(0, BWT_QUERY_SUFFIX, q, qo) == 0)
         bwt_query(3, BWT_QUERY_PREFIX, qo, qo);
   }
   reset_alloc();

   free(q);
   free(qo);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}


void
test_mem_bwt_query_all
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   bwtquery_t ** qv = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      if (bwt_query_all(BWT_QUERY_SUFFIX, q, qv) == 0) {
         bwt_query_all(BWT_QUERY_PREFIX, qv[0], qv);
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      if (bwt_query_all(BWT_QUERY_SUFFIX, q, qv) == 0)
         bwt_query_all(BWT_QUERY_PREFIX, qv[0], qv);
   }
   reset_alloc();

   free(q);
   bwt_free_vec(qv);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_prefix
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   bwtquery_t * qo = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      if (bwt_prefix(0, q, qo) == 0)
         bwt_prefix(3, qo, qo);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      if (bwt_prefix(0, q, qo) == 0)
         bwt_prefix(3, qo, qo);
   }
   reset_alloc();

   free(q);
   free(qo);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_prefix_all
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   bwtquery_t ** qv = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      if (bwt_prefix_all(q, qv) == 0) {
         bwt_prefix_all(qv[0], qv);
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      if (bwt_prefix_all(q, qv) == 0) {
         bwt_prefix_all(qv[0], qv);
      }
   }
   reset_alloc();

   free(q);
   bwt_free_vec(qv);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_bwt_file
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwt_t * bwt_i;

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      bwt_file_write("test30.bwt", bwt);
      bwt_i = bwt_file_read("test30.bwt", txt);
      bwt_free(bwt_i);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      bwt_file_write("test30.bwt", bwt);
      bwt_i = bwt_file_read("test30.bwt", txt);
      bwt_free(bwt_i);
   }
   reset_alloc();

   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_mem_index_bwt[] = {
   {"mem/index_bwt/bwt_build",            test_mem_bwt_build},
   {"mem/index_bwt/bwt_new_query",        test_mem_bwt_new_query},
   {"mem/index_bwt/bwt_dup_query",        test_mem_bwt_dup_query},
   {"mem/index_bwt/bwt_new_vec",          test_mem_bwt_new_vec},
   {"mem/index_bwt/bwt_dup_vec",          test_mem_bwt_dup_vec},
   {"mem/index_bwt/bwt_query",            test_mem_bwt_query},
   {"mem/index_bwt/bwt_query_all",        test_mem_bwt_query_all},
   {"mem/index_bwt/bwt_prefix",           test_mem_bwt_prefix},
   {"mem/index_bwt/bwt_prefix_all",       test_mem_bwt_prefix_all},
   {"mem/index_bwt/bwt_file",             test_mem_bwt_file},
   {NULL, NULL}, // Sentinel. //
};
