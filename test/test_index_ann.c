#include "unittest.h"
#include "index_ann.h"

void
test_ann_build
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
   test_assert(ann_build(1, 1, bwt, sar, 1) == NULL);
   test_assert(ann_build(-1, 1, bwt, sar, 1) == NULL);
   test_assert(ann_build(10, 0, bwt, sar, 1) == NULL);
   test_assert(ann_build(20, -1, bwt, sar, 1) == NULL);
   test_assert(ann_build(3, 4, bwt, sar, 1) == NULL);
   test_assert(ann_build(25, 8, bwt, sar, 1) == NULL);
   test_assert(ann_build(25, 1, NULL, sar, 1) == NULL);
   test_assert(ann_build(25, 1, bwt, NULL, 1) == NULL);
   test_assert(ann_build(25, 1, bwt, sar, 0) == NULL);
   test_assert(ann_build(25, 1, bwt, sar, -1) == NULL);

   ann_t * ann = ann_build(25, 1, bwt, sar, 1);
   unredirect_stderr();
   test_assert_critical(ann != NULL);

   ann_free(ann);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // The random chromosome.
   char * text = malloc(10001);
   char dna[4] = {'A','C','G','T'};
   for (int i = 0; i < 10000; i++) {
      int64_t v = random();
      if (v % 50 == 0) {
         text[i] = 'N';
      } else {
         text[i] = dna[v%4];
      }
   }
   text[10000] = 0;
   
   test_assert(txt_append(text, txt) == 0);
   test_assert(txt_commit_seq("seq0", txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   redirect_stderr();
   ann = ann_build(10, 1, bwt, sar, 1);
   test_assert_critical(sar != NULL);
   unredirect_stderr();

   ann_free(ann);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_ann_query
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
   unredirect_stderr();
   test_assert_critical(ann != NULL);

   /*
   ** 0        10        20   25
   ** ATCGATATCAGCCACTACGAGACAA$
   **    30        40        50
   ** ATCGATATCAGgCACTACGAGACAA$
   **        60        70     77
   ** ATCGATATCAGCCACTACGAtACAA$
   */

   locinfo_t * li = ann_query(0,ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 7);
   test_assert(li->align_cnt == 0);
   free(li);

   li = ann_query(26,ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 11);
   free(li);

   li = ann_query(23,ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 0);
   test_assert(li->neigh_cnt == 0);
   test_assert(li->align_cnt == 0);
   free(li);

   li = ann_query(txt_str_to_pos("seq1:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 11);
   free(li);

   li = ann_query(txt_str_to_pos("seq1:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24-11);
   free(li);

   li = ann_query(txt_str_to_pos("seq2:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 20);
   free(li);

   li = ann_query(txt_str_to_pos("seq2:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24-20);
   free(li);

   li = ann_query(txt_str_to_pos("seq3:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 0);
   free(li);

   li = ann_query(txt_str_to_pos("seq3:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24);
   free(li);

   li = ann_query(txt_str_to_pos("seq4:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24);
   free(li);

   li = ann_query(txt_str_to_pos("seq4:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 0);
   free(li);

   li = ann_query(txt_str_to_pos("seq5:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 10);
   free(li);

   li = ann_query(txt_str_to_pos("seq5:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24-10);
   free(li);

   li = ann_query(txt_str_to_pos("seq6:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 2);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 8);
   free(li);

   li = ann_query(txt_str_to_pos("seq6:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 2);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24-8);
   free(li);

   li = ann_query(txt_str_to_pos("seq7:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 2);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 8);
   free(li);

   li = ann_query(txt_str_to_pos("seq7:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 2);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24-8);
   free(li);

   ann_free(ann);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}


void
test_ann_file
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
   unredirect_stderr();
   test_assert_critical(ann != NULL);

   ann_file_write("test00.ann",ann);
   ann_free(ann);

   ann = ann_file_read("test00.ann");

   locinfo_t * li = ann_query(0,ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 7);
   test_assert(li->align_cnt == 0);
   free(li);

   li = ann_query(26,ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 11);
   free(li);

   li = ann_query(23,ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 0);
   test_assert(li->neigh_cnt == 0);
   test_assert(li->align_cnt == 0);
   free(li);

   li = ann_query(txt_str_to_pos("seq5:1:-",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 1);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 24-10);
   free(li);

   li = ann_query(txt_str_to_pos("seq6:1:+",txt), ann);
   test_assert_critical(li != NULL);
   test_assert(li->dist == 1);
   test_assert(li->neigh_cnt == 2);
   test_assert(li->align_cnt == 1);
   test_assert(li->align_pos[0] == 8);
   free(li);
   

   ann_free(ann);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}


// Define test cases to be run (for export).
const test_case_t test_cases_index_ann[] = {
   {"index_ann/ann_build",            test_ann_build},
   {"index_ann/ann_query",            test_ann_query},
   {"index_ann/ann_file",             test_ann_file},
   {NULL, NULL}, // Sentinel. //
};

