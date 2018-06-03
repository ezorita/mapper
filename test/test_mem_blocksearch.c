#include "unittest.h"
#include "index.h"

void
test_mem_blocksc_trail
(void)
{
   uint8_t * query = NULL;
   index_t * index = index_build("examples/repeats.fa", "test_bs00");
   test_assert_critical(index != NULL);
   char * three = "ATCGATATCAGCCACTACGAtACAA";

   pstree_t * pstree = alloc_stack_tree(1);
   test_assert_critical(pstree != NULL);

   // Prepare query.
   bwtquery_t ** qarray = malloc((strlen(three)+1)*sizeof(bwtquery_t *));
   test_assert_critical(qarray != NULL);
   qarray[0] = bwt_new_query(index->bwt);
   test_assert_critical(qarray[0] != NULL);

   query = sym_str_index(three, index->sym);
   test_assert_critical(query != NULL);
   for (int i = 0; i < strlen(three); i++) {
      qarray[i+1] = bwt_new_query(index->bwt);
      test_assert_critical(qarray[i+1] != NULL);
      bwt_query(query[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   redirect_stderr();
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      blocksc_trail(query, qarray, 25, 1, 0, pstree);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      blocksc_trail(query, qarray, 25, 1, 0, pstree);
   }
   reset_alloc();

   // Free memory.
   for (int i = 0; i < strlen(three)+1; i++) {
      free(qarray[i]);
   }
   free(qarray);
   free(query);
   free_stack_tree(pstree);
   index_free(index);
   unredirect_stderr();
}

void
test_mem_blocksearch_trail_rec
(void)
{
   uint8_t * query = NULL;
   index_t * index = index_build("examples/repeats.fa", "test_bs00");
   test_assert_critical(index != NULL);
   char * three = "ATCGATATCAGCCACTACGAtACAA";

   pstree_t * pstree = alloc_stack_tree(1);
   test_assert_critical(pstree != NULL);

   // Prepare query.
   query = sym_str_index(three, index->sym);
   test_assert_critical(query != NULL);

   redirect_stderr();
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      blocksearch_trail_rec(query, 0, 24, 2, 0, index->bwt, pstree);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      blocksearch_trail_rec(query, 0, 24, 2, 0, index->bwt, pstree);
   }
   reset_alloc();

   free(query);
   free_stack_tree(pstree);
   index_free(index);
   unredirect_stderr();
}


void
test_mem_seqsearch
(void)
{
   uint8_t * query = NULL;
   spath_t * path = calloc(1, sizeof(spath_t));
   test_assert_critical(path != NULL);
   pathstack_t * hits = pathstack_new(10);
   test_assert_critical(hits != NULL);

   index_t * index = index_build("examples/repeats.fa", "test_bs00");
   test_assert_critical(index != NULL);

   char * three = "ATCGATATCAGCCACTACGAtACAA";
   // Prepare query.
   path->bwtq = bwt_new_query(index->bwt);
   test_assert_critical(path->bwtq != NULL);

   query = sym_str_index(three, index->sym);
   test_assert_critical(query != NULL);

   redirect_stderr();
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      seqsearch_fw(*path, query, 0, 24, 1, 0, 0, &hits);
      seqsearch_bw(*path, query, 24, 0, 1, 0, 0, &hits);
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   hits->pos = 0;
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      seqsearch_fw(*path, query, 0, 24, 1, 0, 0, &hits);
      seqsearch_bw(*path, query, 24, 0, 1, 0, 0, &hits);
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   hits->pos = 0;
   reset_alloc();

   free(query);
   free(hits);
   free(path->bwtq);
   free(path);
   index_free(index);
   unredirect_stderr();
}


void
test_mem_scsearch
(void)
{
   uint8_t * query = NULL;
   spath_t * path = calloc(1, sizeof(spath_t));
   test_assert_critical(path != NULL);
   pathstack_t * hits = pathstack_new(10);
   test_assert_critical(hits != NULL);

   index_t * index = index_build("examples/repeats.fa", "test_bs00");
   test_assert_critical(index != NULL);

   char * three = "ATCGATATCAGCCACTACGAtACAA";
   // Prepare query.
   path->bwtq = bwt_new_query(index->bwt);
   test_assert_critical(path->bwtq != NULL);

   query = sym_str_index(three, index->sym);
   test_assert_critical(query != NULL);

   redirect_stderr();
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      scsearch_fw(*path, query, 0, 24, 1, 0, 0, 1, &hits);
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   hits->pos = 0;
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      scsearch_fw(*path, query, 0, 24, 1, 0, 0, 1, &hits);
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   hits->pos = 0;
   reset_alloc();

   free(query);
   free(hits);
   free(path->bwtq);
   free(path);
   index_free(index);
   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_mem_blocksearch[] = {
   {"mem/blocksearch/blocksc_trail",         test_mem_blocksc_trail},
   {"mem/blocksearch/blocksearch_trail_rec", test_mem_blocksearch_trail_rec},
   {"mem/blocksearch/seqsearch",             test_mem_seqsearch},
   {"mem/blocksearch/scsearch",              test_mem_scsearch},
   {NULL, NULL}, // Sentinel. //
};
