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

const test_case_t test_mem_seqread[] = {
   {"mem/seqread/seqread_new",         test_mem_seqread_new},
   {NULL, NULL}, // Sentinel. //
};
