#include "unittest.h"
#include "index.h"

void
test_index_build
(void)
{
   redirect_stderr();
   test_assert(index_build("fakefile.fasta", "test_base") == NULL);
   test_assert(index_build(NULL, "test_base") == NULL);
   test_assert(index_build("examples/repeats.fa", NULL) == NULL);

   index_t * index = index_build("examples/repeats.fa", "test_base00");
   test_assert_critical(index != NULL);
   test_assert(strcmp(index->fname_base, "test_base00") == 0);
   test_assert(index->sym != NULL);
   test_assert(index->txt != NULL);
   test_assert(index->sar != NULL);
   test_assert(index->bwt != NULL);
   test_assert(index->ann == NULL);
   test_assert(index->ann_cnt == 0);

   // Test index content.
   uint8_t seq_one_suf[25] = {0,3,1,2,0,3,0,3,1,0,2,1,1,0,1,3,0,1,2,0,2,0,1,0,0};
   uint8_t seq_two_pre[25] = {0,0,1,0,2,0,2,1,0,3,1,0,1,2,2,0,1,3,0,3,0,2,1,3,0};

   bwtquery_t * q = bwt_new_query(index->bwt);
   for (int i = 0; i < 25; i ++) {
      test_assert(bwt_query(seq_one_suf[i], BWT_QUERY_SUFFIX, q, q) == 0);
   }
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 25);
   test_assert(strcmp(txt_pos_to_str(sar_get(bwt_start(q), index->sar), index->txt), "one:1:+") == 0);
   free(q);

   q = bwt_new_query(index->bwt);
   for (int i = 0; i < 25; i ++) {
      test_assert(bwt_query(seq_two_pre[i], BWT_QUERY_PREFIX, q, q) == 0);
   }
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 25);
   test_assert(strcmp(txt_pos_to_str(sar_get(bwt_start(q), index->sar), index->txt), "two:1:+") == 0);
   free(q);

   index_free(index);
   unredirect_stderr();
}

void
test_index_ann_new
(void)
{
   redirect_stderr();
   index_t * index = index_build("examples/repeats.fa", "test_base01");
   test_assert_critical(index != NULL);

   redirect_stderr();
   index_ann_new(25, 1, 1, index);
   unredirect_stderr();
   test_assert_critical(index->ann != NULL);
   test_assert(index->ann_cnt == 1);


   locinfo_t * lci = ann_query(txt_str_to_pos("one:1:+", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 7);
   test_assert(lci->align_cnt == 0);
   free(lci);

   lci = ann_query(txt_str_to_pos("five:1:+", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 1);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 24);
   free(lci);

   lci = ann_query(txt_str_to_pos("seven:1:-", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 2);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 16);
   free(lci);

   index_free(index);
   unredirect_stderr();
}


void
test_index_read
(void)
{
   redirect_stderr();
   index_t * index = index_build("examples/repeats.fa", "test_base02");
   test_assert_critical(index != NULL);
   test_assert(index->sym != NULL);
   test_assert(index->txt != NULL);
   test_assert(index->sar != NULL);
   test_assert(index->bwt != NULL);
   test_assert(index->ann == NULL);
   test_assert(index->ann_cnt == 0);


   redirect_stderr();
   index_ann_new(25, 1, 1, index);
   unredirect_stderr();
   test_assert_critical(index->ann != NULL);

   index_free(index);

   index = index_read("test_base02");
   test_assert(strcmp(index->fname_base, "test_base02") == 0);
   test_assert(index->sym != NULL);
   test_assert(index->txt != NULL);
   test_assert(index->sar != NULL);
   test_assert(index->bwt != NULL);
   test_assert(index->ann != NULL);
   test_assert(index->ann_cnt == 1);

   // Test index content.
   uint8_t seq_one_suf[25] = {0,3,1,2,0,3,0,3,1,0,2,1,1,0,1,3,0,1,2,0,2,0,1,0,0};
   uint8_t seq_two_pre[25] = {0,0,1,0,2,0,2,1,0,3,1,0,1,2,2,0,1,3,0,3,0,2,1,3,0};

   bwtquery_t * q = bwt_new_query(index->bwt);
   for (int i = 0; i < 25; i ++) {
      test_assert(bwt_query(seq_one_suf[i], BWT_QUERY_SUFFIX, q, q) == 0);
   }
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 25);
   test_assert(strcmp(txt_pos_to_str(sar_get(bwt_start(q), index->sar), index->txt), "one:1:+") == 0);
   free(q);

   q = bwt_new_query(index->bwt);
   for (int i = 0; i < 25; i ++) {
      test_assert(bwt_query(seq_two_pre[i], BWT_QUERY_PREFIX, q, q) == 0);
   }
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 25);
   test_assert(strcmp(txt_pos_to_str(sar_get(bwt_start(q), index->sar), index->txt), "two:1:+") == 0);
   free(q);

   // Test annotation.
   locinfo_t * lci = ann_query(txt_str_to_pos("six:1:+", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 1);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 10);
   free(lci);

   lci = ann_query(txt_str_to_pos("eight:1:+", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 2);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 8);
   free(lci);

   lci = ann_query(txt_str_to_pos("three:1:-", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 1);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 4);
   free(lci);

   index_free(index);

   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_cases_index[] = {
   {"index/index_build",         test_index_build},
   {"index/index_ann_new",       test_index_ann_new},
   {"index/index_read",          test_index_read},
   {NULL, NULL}, // Sentinel. //
};
