#include "unittest.h"
#include "index_sar.h"

void
test_sar_build
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

   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(sar_build(NULL) == NULL);

   int64_t * sa_buf = malloc(100*sizeof(sa_buf));
   test_assert_critical(sa_buf != NULL);
   
   // seq_1
   test_assert(txt_append("TAGCNTCGACA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_get(2, sar) == 8);
   test_assert(sar_get_range(0, strlen(seq_1), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, sar_1, strlen(seq_1)*sizeof(int64_t)) == 0);
   sar_free(sar);
   txt_free(txt);

   // seq_2
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append("TGNGA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TCGAT", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   
   sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_get(4,sar) == 7);
   test_assert(sar_get_range(0, strlen(seq_2), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, sar_2, strlen(seq_2)*sizeof(int64_t)) == 0);
   sar_free(sar);
   txt_free(txt);

   // seq_3
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append("AAaaAAaaA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   
   sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_get(6,sar) == 3);
   test_assert(sar_get_range(0, strlen(seq_3), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, sar_3, strlen(seq_3)*sizeof(int64_t)) == 0);
   sar_free(sar);
   txt_free(txt);

   // seq_4
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append("acgTPGaKcnN", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   
   sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_get(0,sar) == 11);
   test_assert(sar_get_range(0, strlen(seq_4), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, sar_4, strlen(seq_4)*sizeof(int64_t)) == 0);
   sar_free(sar);
   txt_free(txt);

   // seq_5
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TgAgcatGGC", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("acTCGA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   
   sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_get(5,sar) == 3);
   test_assert(sar_get_range(0, strlen(seq_5), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, sar_5, strlen(seq_5)*sizeof(int64_t)) == 0);
   sar_free(sar);
   txt_free(txt);

   // long sequence
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("TGCAAAAAA",txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }
   sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   for (int i = 0; i < 1000; i++) {
      test_assert(sar_get(i, sar) == 9999-i*10);
      test_assert(sar_get(1000+i, sar) == 9998-i*10);
      test_assert(sar_get(2000+i, sar) == 9997-i*10);
      test_assert(sar_get(3000+i, sar) == 9996-i*10);
      test_assert(sar_get(4000+i, sar) == 9995-i*10);
      test_assert(sar_get(5000+i, sar) == 9994-i*10);
      test_assert(sar_get(6000+i, sar) == 9993-i*10);
      test_assert(sar_get(7000+i, sar) == 9992-i*10);
      test_assert(sar_get(8000+i, sar) == 9991-i*10);
      test_assert(sar_get(9000+i, sar) == 9990-i*10);
   }

   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   free(sa_buf);

   unredirect_stderr();
}


void
test_sar_get
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TgAgcatGGC", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("acTCGA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   sar_t * sar = sar_build(txt);
   test_assert(sar_get(10, NULL) == -1);
   test_assert(sar_get(-1, sar)  == -1);
   test_assert(sar_get(100, sar) == -1);

   test_assert(sar_get(0,sar) == 18);
   test_assert(sar_get(1,sar) == 11);
   test_assert(sar_get(2,sar) == 0);
   test_assert(sar_get(3,sar) == 17);
   test_assert(sar_get(4,sar) == 12);
   test_assert(sar_get(5,sar) == 3);
   test_assert(sar_get(6,sar) == 6);
   test_assert(sar_get(7,sar) == 10);
   test_assert(sar_get(8,sar) == 5);
   test_assert(sar_get(9,sar) == 15);
   test_assert(sar_get(10,sar) == 13);
   test_assert(sar_get(11,sar) == 16);
   test_assert(sar_get(12,sar) == 2);
   test_assert(sar_get(13,sar) == 9);
   test_assert(sar_get(14,sar) == 4);
   test_assert(sar_get(15,sar) == 8);
   test_assert(sar_get(16,sar) == 14);
   test_assert(sar_get(17,sar) == 1);
   test_assert(sar_get(18,sar) == 7);

   sar_free(sar);
   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


void
test_sar_get_range
(void)
{
   redirect_stderr();

   int64_t  ref_sa[] = {18,11,0,17,12,3,6,10,5,15,13,16,2,9,4,8,14,1,7};
   int64_t * sa_buf = malloc(30*sizeof(int64_t));
   test_assert_critical(sa_buf != NULL);

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TgAgcatGGC", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("acTCGA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_get_range(0,10, sa_buf, NULL) == -1);
   test_assert(sar_get_range(0,10, NULL, sar) == -1);
   test_assert(sar_get_range(-1, 10, sa_buf, sar)  == -1);
   test_assert(sar_get_range(14, 0, sa_buf, sar) == -1);
   test_assert(sar_get_range(14, -1, sa_buf, sar) == -1);
   test_assert(sar_get_range(14, 20, sa_buf, sar) == -1);
   
   test_assert(sar_get_range(0,19,sa_buf,sar) == 0);
   test_assert(memcmp(sa_buf, ref_sa, 19*sizeof(int64_t)) == 0);
   
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   free(sa_buf);

   unredirect_stderr();
}


void
test_sar_file
(void)
{
   redirect_stderr();

   int64_t  ref_sa[] = {18,11,0,17,12,3,6,10,5,15,13,16,2,9,4,8,14,1,7};
   int64_t * sa_buf = malloc(30*sizeof(int64_t));
   test_assert_critical(sa_buf != NULL);

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TgAgcatGGC", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("acTCGA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_file_write("test00.sar", sar) == 0);
   sar_free(sar);

   sar = sar_file_read("test00.sar");
   test_assert_critical(sar != NULL);
   test_assert(sar_get_range(0,19,sa_buf,sar) == 0);
   test_assert(memcmp(sa_buf, ref_sa, 19*sizeof(int64_t)) == 0);

   sar_free(sar);
   txt_free(txt);

   // long sequence
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("TGCAAAAAA",txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }
   sar = sar_build(txt);
   test_assert_critical(sar != NULL);
   test_assert(sar_file_write("test01.sar", sar) == 0);
   sar_free(sar);

   sar = sar_file_read("test01.sar");
   test_assert_critical(sar != NULL);
   for (int i = 0; i < 1000; i++) {
      test_assert(sar_get(i, sar) == 9999-i*10);
      test_assert(sar_get(1000+i, sar) == 9998-i*10);
      test_assert(sar_get(2000+i, sar) == 9997-i*10);
      test_assert(sar_get(3000+i, sar) == 9996-i*10);
      test_assert(sar_get(4000+i, sar) == 9995-i*10);
      test_assert(sar_get(5000+i, sar) == 9994-i*10);
      test_assert(sar_get(6000+i, sar) == 9993-i*10);
      test_assert(sar_get(7000+i, sar) == 9992-i*10);
      test_assert(sar_get(8000+i, sar) == 9991-i*10);
      test_assert(sar_get(9000+i, sar) == 9990-i*10);
   }

   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   free(sa_buf);

   unredirect_stderr();
}

// Define test cases to be run (for export).
const test_case_t test_cases_index_sar[] = {
   {"index_sar/sar_build",             test_sar_build},
   {"index_sar/sar_get",               test_sar_get},
   {"index_sar/sar_get_range",         test_sar_get_range},
   {"index_sar/sar_file",              test_sar_file},
   {NULL, NULL}, // Sentinel. //
};

