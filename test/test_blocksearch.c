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
   char * seven = "ATCGATATaAGCCACTACGAGACAA";
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

   // Query sequence SEVEN (mistmatch in first half, should not find it).
   query32 = sym_str_index(seven, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(seven); i++) {
      query8[i] = (uint8_t)query32[i];
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   test_assert(blocksc_trail(query8, qarray, 25, 1, 0, pstree) == 0);
   test_assert(pstree->stack->pos == 1);
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
   char * one =   "ATCGATATCAGCCACTACGAGACAA";
   char * three = "ATCGATATCAGCCACTACGAtACAA";
   char * five =  "ATCGATATCAGCCACTACGAGACAc";
   char * seven = "ATCGATATaAGCCACTACGAGACAA";
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
   test_assert(blocksearch_trail_rec(NULL, 0, 24, 2, 0, index->bwt, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, -1, 24, 2, 0, index->bwt, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, 0, -1, 2, 0, index->bwt, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, 12, 9, 2, 0, index->bwt, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, 0, 24, -1, 0, index->bwt, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, -1, index->bwt, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 0, NULL, pstree) == -1);
   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 0, index->bwt, NULL) == -1);

   // Query sequence FIVE
   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 0, index->bwt, pstree) == 0);
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

   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 15, index->bwt, pstree) == 0);
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

   // Query sequence SEVEN.
   query32 = sym_str_index(seven, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(seven); i++) {
      query8[i] = (uint8_t)query32[i];
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 0, index->bwt, pstree) == 0);
   test_assert(pstree->stack->pos == 3);
   for (int i = 0; i < 3; i++) {
      if (pstree->stack->path[i].score == 1) {
         test_assert(pstree->stack->path[i].align[0] == 1<<8);
         char * locus = txt_pos_to_str(sar_get(bwt_start(pstree->stack->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0 || strcmp("eight:1:+", locus) == 0);
         free(locus);
      }
   }
   free(query32);

   // Query sequence ONE.
   query32 = sym_str_index(one, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(one); i++) {
      query8[i] = (uint8_t)query32[i];
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 0, index->bwt, pstree) == 0);
   test_assert(pstree->stack->pos == 8);
   free(query32);

   // Query sequence NONE (in seqsearch_trail_rec 'N' are treated as wildcards)
   query32 = sym_str_index(none, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(none); i++) {
      query8[i] = (uint8_t)query32[i];
      bwt_query(query8[i], BWT_QUERY_SUFFIX, qarray[i], qarray[i+1]);
   }

   test_assert(blocksearch_trail_rec(query8, 0, 24, 2, 0, index->bwt, pstree) == 0);
   test_assert(pstree->stack->pos == 8);
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
test_seqsearch
(void)
{
   uint8_t * query8 = malloc(50);
   test_assert_critical(query8 != NULL);
   spath_t * path = calloc(1, sizeof(spath_t));
   test_assert_critical(path != NULL);
   pathstack_t * hits = pathstack_new(10);
   test_assert_critical(hits != NULL);

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
   char * one =   "ATCGATATCAGCCACTACGAGACAA";
   char * three = "ATCGATATCAGCCACTACGAtACAA";
   char * five =  "ATCGATATCAGCCACTACGAGACAc";
   char * seven = "ATCGATATaAGCCACTACGAGACAA";
   char * none =  "NNNNATATCAGCCACTACGAGACAA";

   // Prepare query.
   path->bwtq = bwt_new_query(index->bwt);
   test_assert_critical(path->bwtq != NULL);

   int32_t * query32 = sym_str_index(five, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(five); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   redirect_stderr();
   // Query sequence FIVE (fw).
   hits->pos = 0;
   test_assert(seqsearch_fw(*path, query8, 0, 24, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<24);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }

   // Query sequence FIVE (bw).
   hits->pos = 0;
   test_assert(seqsearch_bw(*path, query8, 24, 0, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<24);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }

   free(query32);

   // Query sequence THREE (fw).
   query32 = sym_str_index(three, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(three); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   hits->pos = 0;
   test_assert(seqsearch_fw(*path, query8, 0, 24, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<20);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }

   // Query sequence THREE (bw).
   hits->pos = 0;
   test_assert(seqsearch_bw(*path, query8, 24, 0, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<20);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence SEVEN (fw).
   query32 = sym_str_index(seven, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(seven); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   hits->pos = 0;
   test_assert(seqsearch_fw(*path, query8, 0, 24, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 3);
   for (int i = 0; i < 3; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<8);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0 || strcmp("eight:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }

   // Query sequence SEVEN (bw).
   hits->pos = 0;
   test_assert(seqsearch_bw(*path, query8, 24, 0, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 3);
   for (int i = 0; i < 3; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<8);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0 || strcmp("eight:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence ONE.
   query32 = sym_str_index(one, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(one); i++) {
      query8[i] = (uint8_t)query32[i];
   }
   // fw
   hits->pos = 0;
   test_assert(seqsearch_fw(*path, query8, 0, 24, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 8);
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   // bw
   hits->pos = 0;
   test_assert(seqsearch_bw(*path, query8, 24, 0, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 8);
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence NONE (in seqsearch_* 'N' are treated as wildcards)
   query32 = sym_str_index(none, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(none); i++) {
      query8[i] = (uint8_t)query32[i];
   }
   // fw
   hits->pos = 0;
   test_assert(seqsearch_fw(*path, query8, 0, 24, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 8);
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   // bw
   hits->pos = 0;
   test_assert(seqsearch_bw(*path, query8, 24, 0, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 8);
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Free memory.   
   free(query8);
   free(hits);
   free(path->bwtq);
   free(path);
   index_free(index);
   unredirect_stderr();
   
}


void
test_scsearch
(void)
{
   // SC search is Seeq&Construct search. Only searches the left space of the trie,
   // where the boundary is the query sequence.
   uint8_t * query8 = malloc(50);
   test_assert_critical(query8 != NULL);
   spath_t * path = calloc(1, sizeof(spath_t));
   test_assert_critical(path != NULL);
   pathstack_t * hits = pathstack_new(10);
   test_assert_critical(hits != NULL);

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
   char * one =   "ATCGATATCAGCCACTACGAGACAA";
   char * three = "ATCGATATCAGCCACTACGAtACAA";
   char * five =  "ATCGATATCAGCCACTACGAGACAc";
   char * seven = "ATCGATATaAGCCACTACGAGACAA";
   char * eight = "ATCGATATtAGCCACTACGAGACAA";
   char * none =  "NNNNATATCAGCCACTACGAGACAA";

   // Prepare query.
   path->bwtq = bwt_new_query(index->bwt);
   test_assert_critical(path->bwtq != NULL);

   int32_t * query32 = sym_str_index(five, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(five); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   redirect_stderr();
   // Query sequence FIVE (sc).
   hits->pos = 0;
   test_assert(scsearch_fw(*path, query8, 0, 24, 1, 0, 0, 1, &hits) == 0);
   test_assert(hits->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<24);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence THREE (sc).
   query32 = sym_str_index(three, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(three); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   hits->pos = 0;
   test_assert(scsearch_fw(*path, query8, 0, 24, 1, 0, 0, 1, &hits) == 0);
   test_assert(hits->pos == 2);
   for (int i = 0; i < 2; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<20);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence SEVEN (sc).
   query32 = sym_str_index(seven, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(seven); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   hits->pos = 0;
   test_assert(scsearch_fw(*path, query8, 0, 24, 1, 0, 0, 1, &hits) == 0);
   test_assert(hits->pos == 1);
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence EIGHT (sc).
   query32 = sym_str_index(eight, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(eight); i++) {
      query8[i] = (uint8_t)query32[i];
   }

   hits->pos = 0;
   test_assert(scsearch_fw(*path, query8, 0, 24, 1, 0, 0, 1, &hits) == 0);
   test_assert(hits->pos == 3);
   for (int i = 0; i < 3; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<8);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("one:1:+", locus) == 0 || strcmp("seven:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence ONE.
   query32 = sym_str_index(one, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(one); i++) {
      query8[i] = (uint8_t)query32[i];
   }
   // sc
   hits->pos = 0;
   test_assert(scsearch_fw(*path, query8, 0, 24, 1, 0, 0, 1, &hits) == 0);
   test_assert(hits->pos == 3);
   for (int i = 0; i < 3; i++) {
      if (hits->path[i].score == 1) {
         test_assert(hits->path[i].align[0] == 1<<8 || hits->path[i].align[0] == 1<<10);
         char * locus = txt_pos_to_str(sar_get(bwt_start(hits->path[i].bwtq), index->sar), index->txt);
         test_assert(strcmp("six:1:+", locus) == 0 || strcmp("seven:1:+", locus) == 0);
         free(locus);
      }
   }
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Query sequence NONE (in seqsearch_* 'N' are treated as wildcards)
   query32 = sym_str_index(none, index->sym);
   test_assert_critical(query32 != NULL);
   for (int i = 0; i < strlen(none); i++) {
      query8[i] = (uint8_t)query32[i];
   }
   // sc
   hits->pos = 0;
   test_assert(seqsearch_fw(*path, query8, 0, 24, 1, 0, 0, &hits) == 0);
   test_assert(hits->pos == 8);
   for (int i = 0; i < hits->pos; i++) {
      free(hits->path[i].bwtq);
   }
   free(query32);

   // Free memory.   
   free(query8);
   free(hits);
   free(path->bwtq);
   free(path);
   index_free(index);
   unredirect_stderr();
   
}



// Define test cases to be run (for export).
const test_case_t test_cases_blocksearch[] = {
   {"blocksearch/blocksc_trail",         test_blocksc_trail},
   {"blocksearch/blocksearch_trail_rec", test_blocksearch_trail_rec},
   {"blocksearch/seqsearch",             test_seqsearch},
   {"blocksearch/scsearch",              test_scsearch},
   {NULL, NULL}, // Sentinel. //
};
