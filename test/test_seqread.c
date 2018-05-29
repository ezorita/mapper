#include "unittest.h"
#include "seqread.h"

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
   test_assert(seqread_new(tag, NULL, qscore0, sym) == NULL);
   test_assert(seqread_new(tag, seq0, qscore0, NULL) == NULL);
   test_assert(seqread_new(tag, seq1, qscore0, sym) == NULL);
   test_assert(seqread_new(tag, seq0, qscore1, sym) == NULL);

   read = seqread_new(tag, seq0, qscore0, sym);
   test_assert_critical(read != NULL);
   test_assert(memcmp(seqread_sym(read), sym0, 15) == 0);
   seqread_free(read);

   read = seqread_new(tag, seq1, qscore1, sym);
   test_assert_critical(read != NULL);
   test_assert(memcmp(seqread_sym(read), sym1, 20) == 0);
   seqread_free(read);

   read = seqread_new(NULL, seq1, NULL, sym);
   test_assert_critical(read != NULL);
   test_assert(seqread_tag(read) == NULL);
   test_assert(seqread_qscore(read) == NULL);
   seqread_free(read);

   unredirect_stderr();

   sym_free(sym);
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
   seqread_t * read = seqread_new(tag, seq1, qscore1, sym);
   test_assert_critical(read != NULL);
   test_assert(seqread_len(read) == 20);
   test_assert(strcmp(seqread_tag(read), tag) == 0);
   test_assert(strcmp(seqread_seq(read), seq1) == 0);
   test_assert(strcmp(seqread_qscore(read), qscore1) == 0);
   test_assert(memcmp(seqread_sym(read), sym1,20) == 0);
   seqread_free(read);
   sym_free(sym);
}


const test_case_t test_cases_seqread[] = {
   {"seqread/seqread_new",         test_seqread_new},
   {"seqread/seqread_helpers",     test_seqread_helpers},
   {NULL, NULL}, // Sentinel. //
};
