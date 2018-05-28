#include "unittest.h"
#include "index_ann.h"

void
test_mem_ann_build
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("ATCGATATCAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("ATCGATATCAGgCACTACGAGACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("ATCGATATCAGCCACTACGAtACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("cTCGATATCAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("ATCGATATCAGCCACTACGAGACAc", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("ATCGATATCAcCCACTACGAGACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("ATCGATATaAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("ATCGATATtAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);
   
   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      ann_t * ann = ann_build(25, 1, bwt, sar, 1);
      ann_free(ann);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      ann_t * ann = ann_build(25, 1, bwt, sar, 1);
      ann_free(ann);
   }
   reset_alloc();

   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_ann_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("ATCGATATCAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_append("ATCGATATCAGgCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq1",txt) == 0);
   test_assert(txt_append("ATCGATATCAGCCACTACGAtACAA", txt) == 0);
   test_assert(txt_commit_seq("seq2",txt) == 0);
   test_assert(txt_append("cTCGATATCAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq3",txt) == 0);
   test_assert(txt_append("ATCGATATCAGCCACTACGAGACAc", txt) == 0);
   test_assert(txt_commit_seq("seq4",txt) == 0);
   test_assert(txt_append("ATCGATATCAcCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq5",txt) == 0);
   test_assert(txt_append("ATCGATATaAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq6",txt) == 0);
   test_assert(txt_append("ATCGATATtAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq7",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);
   
   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   redirect_stderr();
   ann_t * ann = ann_build(25, 1, bwt, sar, 1);
   test_assert_critical(ann != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      free(ann_query(26,ann));
      free(ann_query(23,ann));
      free(ann_query(0,ann));
      free(ann_query(32,ann));
      free(ann_query(541,ann));
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      free(ann_query(26,ann));
      free(ann_query(23,ann));
      free(ann_query(0,ann));
      free(ann_query(32,ann));
      free(ann_query(541,ann));
   }
   reset_alloc();

   ann_free(ann);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}


void
test_mem_ann_file
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("ATCGATATCAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_append("ATCGATATCAGgCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq1",txt) == 0);
   test_assert(txt_append("ATCGATATCAGCCACTACGAtACAA", txt) == 0);
   test_assert(txt_commit_seq("seq2",txt) == 0);
   test_assert(txt_append("cTCGATATCAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq3",txt) == 0);
   test_assert(txt_append("ATCGATATCAGCCACTACGAGACAc", txt) == 0);
   test_assert(txt_commit_seq("seq4",txt) == 0);
   test_assert(txt_append("ATCGATATCAcCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq5",txt) == 0);
   test_assert(txt_append("ATCGATATaAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq6",txt) == 0);
   test_assert(txt_append("ATCGATATtAGCCACTACGAGACAA", txt) == 0);
   test_assert(txt_commit_seq("seq7",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);
   
   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   redirect_stderr();
   ann_t * ann = ann_build(25, 1, bwt, sar, 1);
   test_assert_critical(ann != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      ann_file_write("test30.ann", ann);
      ann_t * ann_i = ann_file_read("test30.ann");
      ann_free(ann_i);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 100; i++) {
      set_alloc_failure_countdown_to(i);
      ann_file_write("test30.ann", ann);
      ann_t * ann_i = ann_file_read("test30.ann");
      ann_free(ann_i);
   }
   reset_alloc();

   ann_free(ann);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}



// Define test cases to be run (for export).
const test_case_t test_mem_index_ann[] = {
   {"mem/index_ann/ann_build",            test_mem_ann_build},
   {"mem/index_ann/ann_query",            test_mem_ann_query},
   {"mem/index_ann/ann_file",             test_mem_ann_file},
   {NULL, NULL}, // Sentinel. //
};
