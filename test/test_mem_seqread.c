#include "unittest.h"
#include "seqread.h"

void
test_mem_seqread_new
(void)
{
   char * tag     = "sequence_name";
   char * seq0    = "ATCGANTANTAGCGN";
   char * qscore0 = "012391923823893";

   sym_t * sym = sym_new_dna();

   seqread_t * read = NULL;

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      read = seqread_new(tag, seq0, qscore0, sym);
      seqread_free(read);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      read = seqread_new(tag, seq0, qscore0, sym);
      seqread_free(read);
   }
   reset_alloc();
   unredirect_stderr();
   sym_free(sym);
}

void
test_mem_seqread_stack
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      gstack_t * stack = seqread_stack(1);
      gstack_free(stack);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      gstack_t * stack = seqread_stack(1);
      gstack_free(stack);
   }
   reset_alloc();

   unredirect_stderr();
   sym_free(sym);
}

void
test_mem_seqread_stack_pop
(void)
{
   char * tag     = "sequence_name";
   char * seq0    = "ATCGANTANTAGCGN";
   char * seq1    = "ACGTACJFDKACGGAACTGT";
   char * qscore0 = "012391923823893";
   char * qscore1 = "19239194385439334402";

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   gstack_t * stack = seqread_stack(1);
   test_assert_critical(stack != NULL);

   for (int i = 0; i < 50; i++) {
      seqread_t * read = seqread_new(tag, seq0, qscore0, sym);
      test_assert_critical(read != NULL);
      test_assert(gstack_push(read, stack) == 0);

      read = seqread_new(tag, seq1, qscore1, sym);
      test_assert_critical(read != NULL);
      test_assert(gstack_push(read, stack) == 0);
   }

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      seqread_t * read = seqread_pop(stack);
      seqread_free(read);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      seqread_t * read = seqread_pop(stack);
      seqread_free(read);
   }
   reset_alloc();

   unredirect_stderr();
   gstack_free(stack);
   sym_free(sym);
}


const test_case_t test_mem_seqread[] = {
   {"mem/seqread/seqread_new",         test_mem_seqread_new},
   {"mem/seqread/seqread_stack",       test_mem_seqread_stack},
   {"mem/seqread/seqread_stack_pop",   test_mem_seqread_stack_pop},
   {NULL, NULL}, // Sentinel. //
};
