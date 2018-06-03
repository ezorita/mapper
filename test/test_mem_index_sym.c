#include "unittest.h"
#include "index_sym.h"

void
test_mem_sym_new
(void)
{
   sym_t * sym = NULL;
   char * alph[] = {"Aa","Bb","Cc","Dd",NULL};
   char * comp[] = {"AB", "BA", "CD", "KQ", NULL};
   // Test without expecting any output, the goal is not to have any SIGSEGV.
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym = sym_new(alph, comp, 0);
      sym_free(sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sym = sym_new(alph, comp, 0);
      sym_free(sym);
   }
   reset_alloc();
   unredirect_stderr();
}

void
test_mem_sym_set_complement
(void)
{
   char * alph[] = {"Aa","Bb","Cc","Dd",NULL};
   char * comp[] = {"AB","CD",NULL};

   redirect_stderr();
   sym_t * sym = sym_new(alph, NULL, 0);
   test_assert_critical(sym != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_set_complement(comp, sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sym_set_complement(comp, sym);
   }
   reset_alloc();

   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_sym_character
(void)
{
   char * alph[] = {"Aa","Bb","Cc","Dd","Ee",NULL};

   redirect_stderr();

   sym_t * sym = sym_new(alph, NULL, 0);
   test_assert_critical(sym != NULL);
   
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_character(0, sym);
      sym_character(1, sym);
      sym_character(2, sym);
      sym_character(3, sym);
      sym_character(4, sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sym_character(0, sym);
      sym_character(1, sym);
      sym_character(2, sym);
      sym_character(3, sym);
      sym_character(4, sym);
   }
   reset_alloc();

   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_sym_index
(void)
{
   sym_t * sym;
   char * alph[] = {"Aa","Bb","Cc","Dd","Ee",NULL};

   redirect_stderr();

   sym = sym_new(alph, NULL, 0);
   test_assert_critical(sym != NULL);
   
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_index('A', sym);
      sym_index('a', sym);
      sym_index('B', sym);
      sym_index('b', sym);
      sym_index('C', sym);
      sym_index('c', sym);
      sym_index('D', sym);
      sym_index('d', sym);
      sym_index('E', sym);
      sym_index('e', sym);
      sym_index('N', sym);
      sym_index('n', sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      sym_index('A', sym);
      sym_index('a', sym);
      sym_index('B', sym);
      sym_index('b', sym);
      sym_index('C', sym);
      sym_index('c', sym);
      sym_index('D', sym);
      sym_index('d', sym);
      sym_index('E', sym);
      sym_index('e', sym);
      sym_index('N', sym);
      sym_index('n', sym);
   }
   reset_alloc();

   sym_free(sym);
   unredirect_stderr();

}

void
test_mem_sym_str_index
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   char * text = "ATGCATCGTAGCA";

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      uint8_t * s = sym_str_index(text, sym);
      free(s);   
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      uint8_t * s = sym_str_index(text, sym);
      free(s);   
   }
   reset_alloc();

   sym_free(sym);
}


void
test_mem_sym_is_canonical
(void)
{
   sym_t * sym;
   char * alph[] = {"Aa","Bb","Cc","Dd","Ee",NULL};

   redirect_stderr();

   sym = sym_new(alph, NULL, 0);
   test_assert_critical(sym != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_is_canonical('A',sym);
      sym_is_canonical('a',sym);
      sym_is_canonical('B',sym);
      sym_is_canonical('b',sym);
      sym_is_canonical('C',sym);
      sym_is_canonical('c',sym);
      sym_is_canonical('D',sym);
      sym_is_canonical('d',sym);
      sym_is_canonical('E',sym);
      sym_is_canonical('e',sym);
      sym_is_canonical('N',sym);
      sym_is_canonical('n',sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      sym_is_canonical('A',sym);
      sym_is_canonical('a',sym);
      sym_is_canonical('B',sym);
      sym_is_canonical('b',sym);
      sym_is_canonical('C',sym);
      sym_is_canonical('c',sym);
      sym_is_canonical('D',sym);
      sym_is_canonical('d',sym);
      sym_is_canonical('E',sym);
      sym_is_canonical('e',sym);
      sym_is_canonical('N',sym);
      sym_is_canonical('n',sym);
   }
   reset_alloc();

   sym_free(sym);
   unredirect_stderr();

}

void
test_mem_sym_complement
(void)
{
   char * alph[] = {"Aa","Bb","Cc","Dd",NULL};
   char * comp[] = {"AB","CD",NULL};

   redirect_stderr();

   sym_t * sym;
   sym = sym_new(alph, comp, 0);
   test_assert_critical(sym != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_complement(0, sym);
      sym_complement(1, sym);
      sym_complement(2, sym);
      sym_complement(3, sym);
      sym_complement(4, sym);
      sym_complement(5, sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      sym_complement(0, sym);
      sym_complement(1, sym);
      sym_complement(2, sym);
      sym_complement(3, sym);
      sym_complement(4, sym);
      sym_complement(5, sym);

   }
   reset_alloc();

   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_sym_count
(void)
{
   char * alph[] = {"Aa","Bb","Cc","Dd",NULL};
   char * comp[] = {"AB","CD",NULL};

   redirect_stderr();

   sym_t * sym;
   sym = sym_new(alph, comp, 0);
   test_assert_critical(sym != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_count(sym);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sym_count(sym);
   }
   reset_alloc();

   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_sym_file
(void)
{
   char * alph[] = {"Aa","Bb","Cc","Dd",NULL};

   redirect_stderr();

   sym_t * sym_o, * sym_i;
   sym_o = sym_new(alph, NULL, 0);
   test_assert_critical(sym_o != NULL);
   
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sym_file_write("test20.sym", sym_o);
      sym_i = sym_file_read("test20.sym");
      sym_free(sym_i);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      sym_file_write("test20.sym", sym_o);
      sym_i = sym_file_read("test20.sym");
      sym_free(sym_i);
   }
   reset_alloc();

   sym_free(sym_o);
   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_mem_index_sym[] = {
   {"mem/index_sym/sym_new",             test_mem_sym_new},
   {"mem/index_sym/sym_set_complement",  test_mem_sym_set_complement},
   {"mem/index_sym/sym_character",       test_mem_sym_character},
   {"mem/index_sym/sym_complement",      test_mem_sym_complement},
   {"mem/index_sym/sym_index",           test_mem_sym_index},
   {"mem/index_sym/sym_str_index",       test_mem_sym_str_index},
   {"mem/index_sym/sym_is_canonical",    test_mem_sym_is_canonical},
   {"mem/index_sym/sym_count",           test_mem_sym_count},
   {"mem/index_sym/sym_file",            test_mem_sym_file},
   {NULL, NULL}, // Sentinel. //
};
