#include "unittest.h"
#include "seqread.h"
#include "index_sym.h"

void
test_seqread_new
(void)
{
   char * tag     = "sequence_name";
   char * seq0    = "ATCGANTANTAGCGN";
   char * seq1    = "ACGTACJFDKACGGAACTGT";
   char * qscore0 = "012391923823893";
   char * qscore1 = "19239194385439334402";
   uint8_t sym0[15] = {0,3,1,2,0,4,3,0,4,3,0,2,1,2,4};
   uint8_t sym1[20] = {0,1,2,3,0,1,4,4,4,4,0,1,2,2,0,0,1,3,2,3};

   sym_t * sym = sym_new_dna();

   seqread_t * read = NULL;

   redirect_stderr();
   test_assert(seqread_new(tag, NULL, qscore0) == NULL);
   test_assert(seqread_new(tag, seq1, qscore0) == NULL);
   test_assert(seqread_new(tag, seq0, qscore1) == NULL);

   read = seqread_new(tag, seq0, qscore0);
   test_assert_critical(read != NULL);
   uint8_t * syms = sym_str_index(seqread_seq(read), sym);
   test_assert(memcmp(syms, sym0, 15) == 0);
   free(syms);
   seqread_free(read);

   read = seqread_new(tag, seq1, qscore1);
   test_assert_critical(read != NULL);
   syms = sym_str_index(seqread_seq(read), sym);
   test_assert(memcmp(syms, sym1, 20) == 0);
   free(syms);
   seqread_free(read);

   read = seqread_new(NULL, seq1, NULL);
   test_assert_critical(read != NULL);
   test_assert(seqread_tag(read) == NULL);
   test_assert(seqread_qscore(read) == NULL);
   seqread_free(read);

   unredirect_stderr();

   sym_free(sym);
}

void
test_seqread_stack
(void)
{
   char * tag     = "sequence_name";
   char * seq0    = "ATCGANTANTAGCGN";
   char * seq1    = "ACGTACJFDKACGGAACTGT";
   char * qscore0 = "012391923823893";
   char * qscore1 = "19239194385439334402";

   gstack_t * stack = seqread_stack(1);
   test_assert_critical(stack != NULL);

   seqread_t * read = seqread_new(tag, seq0, qscore0);
   test_assert_critical(read != NULL);
   test_assert(gstack_push(read, stack) == 0);
   test_assert(gstack_num_elm(stack) == 1);

   read = seqread_new(tag, seq1, qscore1);
   test_assert_critical(read != NULL);
   test_assert(gstack_push(read, stack) == 0);
   test_assert(gstack_num_elm(stack) == 2);

   redirect_stderr();
   test_assert(seqread_get(10, stack) == NULL);
   test_assert(seqread_get(0, NULL) == NULL);   
   unredirect_stderr();
   
   read = seqread_get(0, stack);
   test_assert_critical(read != NULL);
   char * s = seqread_seq(read);
   test_assert(s != NULL);
   test_assert(strcmp(seq0, s) == 0);

   read = seqread_pop(stack);
   test_assert_critical(read != NULL);
   test_assert(gstack_num_elm(stack) == 1);
   s = seqread_seq(read);
   test_assert(s != NULL);
   test_assert(strcmp(seq1, s) == 0);
   seqread_free(read);

   read = seqread_pop(stack);
   test_assert_critical(read != NULL);
   test_assert(gstack_num_elm(stack) == 0);
   s = seqread_seq(read);
   test_assert(s != NULL);
   test_assert(strcmp(seq0, s) == 0);
   seqread_free(read);

   read = seqread_pop(stack);
   test_assert(read == NULL);

   gstack_free(stack);
}

void
test_seqread_helpers
(void)
{
   char * tag     = "sequence_name";   
   char * seq1    = "ACGTACJFDKACGGAACTGT";
   char * qscore1 = "19239194385439334402";
   uint8_t sym1[20] = {0,1,2,3,0,1,4,4,4,4,0,1,2,2,0,0,1,3,2,3};

   sym_t * sym = sym_new_dna();
   seqread_t * read = seqread_new(tag, seq1, qscore1);
   test_assert_critical(read != NULL);
   test_assert(seqread_len(read) == 20);
   test_assert(strcmp(seqread_tag(read), tag) == 0);
   test_assert(strcmp(seqread_seq(read), seq1) == 0);
   test_assert(strcmp(seqread_qscore(read), qscore1) == 0);
   uint8_t * syms = sym_str_index(seqread_seq(read), sym);
   test_assert(memcmp(syms, sym1, 20) == 0);
   free(syms);
   seqread_free(read);
   sym_free(sym);
}


const test_case_t test_cases_seqread[] = {
   {"seqread/seqread_new",         test_seqread_new},
   {"seqread/seqread_stack",       test_seqread_stack},
   {"seqread/seqread_helpers",     test_seqread_helpers},

   {NULL, NULL}, // Sentinel. //
};
