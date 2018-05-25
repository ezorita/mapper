#include <unistd.h>
#include <getopt.h>
#include "unittest.h"
#include "index.h"
#include "blocksearch.h"

void
test_blocksc_trail
(void)
{
   uint8_t * query8 = malloc(50);
   test_assert_critical(query8 != NULL);
   index_t * index = index_build("examples/repeats.fa", "test_bs00");
   test_assert_critical(index != NULL);

   /*
     >one
     ATCGATATCAGCCACTACGAGACAA
     >two
     ATCGATATCAGgCACTACGAGACAA
     >three
     ATCGATATCAGCCACTACGAtACAA
     >four
     cTCGATATCAGCCACTACGAGACAA
     >five
     ATCGATATCAGCCACTACGAGACAc
     >six
     ATCGATATCAcCCACTACGAGACAA
     >seven
     ATCGATATaAGCCACTACGAGACAA
     >eight
     ATCGATATtAGCCACTACGAGACAA
   */
   char * three = "ATCGATATCAGCCACTACGAtACAA";
   char * five =  "ATCGATATCAGCCACTACGAGACAc";
   char * none =  "NNNNATATCAGCCACTACGAGACAA";

   pstree_t * pstree = alloc_stack_tree(1);
   test_assert_critical(pstree != NULL);

   // Prepare query.
   bwtquery_t ** qarray = malloc((strlen(five)+1)*sizeof(bwtquery_t *));
   test_assert_critical(qarray != NULL);
   qarray[0] = bwt_new_query(index->bwt);
   test_assert_critical(qarray[0] != NULL);

   int32_t * query32 = sym_str_index(five, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(five); i++) {
      query8[i] = (uint8_t)query32[i];
      qarray[i+1] = bwt_new_query(index->bwt);
      test_assert_critical(qarray[i+1] != NULL);
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   // Invalid arguments.
   redirect_stderr();
   test_assert(blocksc_trail(NULL, qarray, 25, 1, 0, pstree) == -1);
   test_assert(blocksc_trail(query8, NULL, 25, 1, 0, pstree) == -1);
   test_assert(blocksc_trail(query8, qarray, 0, 1, 0, pstree) == -1);
   test_assert(blocksc_trail(query8, qarray, 25, -1, 0, pstree) == -1);
   test_assert(blocksc_trail(query8, qarray, 25, 1, -1, pstree) == -1);
   test_assert(blocksc_trail(query8, qarray, 25, 1, 0, NULL) == -1);

   // Query sequence FIVE
   test_assert(blocksc_trail(query8, qarray, 25, 1, 0, pstree) == 0);
   test_assert(pstree->stack->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (pstree->stack->path[i].score == 1) {
         test_assert(pstree->stack->path[i].align[0] == 1<<24);
         char * locus = txt_pos_to_str(sar_get(bwt_start(pstree->stack->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   free(query32);

   // Query sequence THREE
   query32 = sym_str_index(three, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(three); i++) {
      query8[i] = (uint8_t)query32[i];
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   test_assert(blocksc_trail(query8, qarray, 25, 1, 15, pstree) == 0);
   test_assert(pstree->stack->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (pstree->stack->path[i].score == 1) {
         test_assert(pstree->stack->path[i].align[0] == 1<<20);
         char * locus = txt_pos_to_str(sar_get(bwt_start(pstree->stack->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   free(query32);

   // Query sequence NONE
   query32 = sym_str_index(none, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(none); i++) {
      query8[i] = (uint8_t)query32[i];
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   test_assert(blocksc_trail(query8, qarray, 25, 1, 0, pstree) == 0);
   test_assert(pstree->stack->pos == 0);
   free(query32);

   // Free memory.   
   for (int i = 0; i < strlen(five)+1; i++) {
      free(qarray[i]);
   }
   free(qarray);
   free(query8);
   free_stack_tree(pstree);
   index_free(index);
   unredirect_stderr();
}

void
test_blocksearch_trail_rec
(void)
{

}

// Define test cases to be run (for export).
const test_case_t test_cases_blocksearch[] = {
   {"blocksearch/blocksc_trail",      test_blocksc_trail},
   {NULL, NULL}, // Sentinel. //
};
