#include "unittest.h"
#include "seed.h"

void
test_seed_new
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   
   bwtquery_t * q = bwt_new_query(index->bwt);
   test_assert_critical(q != NULL);
   
   seqread_t * read = seqread_new("one-exact", "ATCGATATCAGCCACTACGAGACAA", NULL);
   test_assert_critical(read != NULL);
   
   redirect_stderr();
   test_assert(seed_new(-1, 10, q, read) == NULL);
   test_assert(seed_new(0, -10, q, read) == NULL);
   test_assert(seed_new(11, 10, q, read) == NULL);
   test_assert(seed_new(0, 25, NULL, read) == NULL);
   test_assert(seed_new(0, 25, q, NULL) == NULL);

   seed_t * seed = seed_new(0, 25, q, read);
   test_assert_critical(seed != NULL);

   unredirect_stderr();

   seed_free(seed);
   seqread_free(read);
   index_free(index);
}

void
test_seed_stack
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   
   bwtquery_t * q0 = bwt_new_query(index->bwt);
   test_assert_critical(q0 != NULL);
   bwtquery_t * q1 = bwt_new_query(index->bwt);
   test_assert_critical(q1 != NULL);

   
   seqread_t * read = seqread_new("one-exact", "ATCGATATCAGCCACTACGAGACAA", NULL);
   test_assert_critical(read != NULL);
   
   gstack_t * stack = seed_stack(1);
   test_assert_critical(stack != NULL);
   
   seed_t * seed0 = seed_new(0, 25, q0, read);
   test_assert_critical(seed0 != NULL);
   seed_t * seed1 = seed_new(0, 25, q1, read);
   test_assert_critical(seed1 != NULL);

   redirect_stderr();
   test_assert(gstack_push(seed0, stack) == 0);
   test_assert(gstack_push(seed1, stack) == 0);

   seed_t * get0 = seed_get(0, stack);
   test_assert(get0 == seed0);
   seed_t * get1 = seed_get(1, stack);
   test_assert(get1 == seed1);

   seed_t * pop1 = seed_pop(stack);
   test_assert(pop1 == seed1);
   seed_free(pop1);
   unredirect_stderr();
   
   gstack_free(stack);
   seqread_free(read);
   index_free(index);
}


void
test_seed_next_mem
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);

   seqread_t * read0 = seqread_new("one-exact", "ATCGATATCAGCCACTACGAGACAA", NULL);
   test_assert_critical(read0 != NULL);
   seqread_t * read1 = seqread_new("one-exact", "ATCGATATCAGCTACTACGAGACAA", NULL);
   test_assert_critical(read1 != NULL);

   seed_t * s0 = seed_next_mem(NULL, read0, index);
   test_assert_critical(s0 != NULL);
   bwtquery_t * q0 = seed_bwtq(s0);
   test_assert_critical(q0 != NULL);
   test_assert(bwt_size(q0) == 1);
   test_assert(bwt_depth(q0) == 25);
   char * pos = txt_pos_to_str(sar_get(bwt_start(q0),index->sar),index->txt);
   test_assert(strcmp("one:1:+", pos) == 0);
   free(pos);
   test_assert(seed_beg(s0) == 0);
   test_assert(seed_end(s0) == 25);
   test_assert(seed_next_mem(s0, read0, index) == NULL);
   seed_free(s0);

   redirect_stderr();
   test_assert(seed_next_mem(NULL, NULL, index) == NULL);
   test_assert(seed_next_mem(NULL, read1, NULL) == NULL);
   test_assert(seed_next_mem(NULL, NULL, NULL) == NULL);
   unredirect_stderr();

   s0 = seed_next_mem(NULL, read1, index);
   test_assert_critical(s0 != NULL);
   q0 = seed_bwtq(s0);
   test_assert_critical(q0 != NULL);
   test_assert(bwt_size(q0) == 3);
   test_assert(bwt_depth(q0) == 12);
   test_assert(seed_beg(s0) == 0);
   test_assert(seed_end(s0) == 12);
   
   seed_t * s1 = seed_next_mem(s0, read1, index);
   test_assert_critical(s1 != NULL);
   // finds GCTA in reverse complement of eight
   bwtquery_t * q1 = seed_bwtq(s1);
   test_assert_critical(q1 != NULL);
   test_assert(bwt_size(q1) == 1);
   test_assert(bwt_depth(q1) == 4);
   test_assert(seed_beg(s1) == 10);
   test_assert(seed_end(s1) == 14);
   
   seed_t * s2 = seed_next_mem(s1, read1, index);
   test_assert_critical(s2 != NULL);
   bwtquery_t * q2 = seed_bwtq(s2);
   test_assert_critical(q2 != NULL);
   test_assert(bwt_size(q2) == 8);
   test_assert(bwt_depth(q2) == 4);
   test_assert(seed_beg(s2) == 11);
   test_assert(seed_end(s2) == 15);

   seed_t * s3 = seed_next_mem(s2, read1, index);
   test_assert_critical(s3 != NULL);
   bwtquery_t * q3 = seed_bwtq(s3);
   test_assert_critical(q3 != NULL);
   test_assert(bwt_size(q3) == 6);
   test_assert(bwt_depth(q3) == 12);
   test_assert(seed_beg(s3) == 13);
   test_assert(seed_end(s3) == 25);
   
   test_assert(seed_next_mem(s3, read1, index) == NULL);

   seed_free(s0);
   seed_free(s1);
   seed_free(s2);
   seed_free(s3);
   seqread_free(read0);
   seqread_free(read1);
   index_free(index);
}


void
test_seed_mems
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   
   seqread_t * read0 = seqread_new("one-exact", "ATCGATATCAGCCACTACGAGACAA", NULL);
   test_assert_critical(read0 != NULL);
   seqread_t * read1 = seqread_new("one-exact", "ATCGATATCAGCTACTACGAGACAA", NULL);
   test_assert_critical(read1 != NULL);
   
   redirect_stderr();
   test_assert(seed_mems(NULL, index) == NULL);
   test_assert(seed_mems(read0, NULL) == NULL);
   test_assert(seed_mems(NULL, NULL) == NULL);
   unredirect_stderr();

   gstack_t * stack0 = seed_mems(read0, index);
   test_assert_critical(stack0 != NULL);
   test_assert(gstack_num_elm(stack0) == 1);
   seed_t * s0 = seed_pop(stack0);
   test_assert_critical(s0 != NULL);
   bwtquery_t * q0 = seed_bwtq(s0);
   test_assert_critical(q0 != NULL);
   test_assert(bwt_size(q0) == 1);
   test_assert(bwt_depth(q0) == 25);
   char * pos = txt_pos_to_str(sar_get(bwt_start(q0),index->sar),index->txt);
   test_assert(strcmp("one:1:+", pos) == 0);
   free(pos);
   test_assert(seed_beg(s0) == 0);
   test_assert(seed_end(s0) == 25);
   seed_free(s0);
   gstack_free(stack0);

   gstack_t * stack1 = seed_mems(read1, index);
   test_assert_critical(stack1 != NULL);
   test_assert(gstack_num_elm(stack1) == 4);

   seed_t * s3 = seed_pop(stack1);
   test_assert_critical(s3 != NULL);
   bwtquery_t * q3 = seed_bwtq(s3);
   test_assert_critical(q3 != NULL);
   test_assert(bwt_size(q3) == 6);
   test_assert(bwt_depth(q3) == 12);
   test_assert(seed_beg(s3) == 13);
   test_assert(seed_end(s3) == 25);

   seed_t * s2 = seed_pop(stack1);
   test_assert_critical(s2 != NULL);
   bwtquery_t * q2 = seed_bwtq(s2);
   test_assert_critical(q2 != NULL);
   test_assert(bwt_size(q2) == 8);
   test_assert(bwt_depth(q2) == 4);
   test_assert(seed_beg(s2) == 11);
   test_assert(seed_end(s2) == 15);

   seed_t * s1 = seed_pop(stack1);
   test_assert_critical(s1 != NULL);
   bwtquery_t * q1 = seed_bwtq(s1);
   test_assert_critical(q1 != NULL);
   test_assert(bwt_size(q1) == 1);
   test_assert(bwt_depth(q1) == 4);
   test_assert(seed_beg(s1) == 10);
   test_assert(seed_end(s1) == 14);

   s0 = seed_pop(stack1);
   test_assert_critical(s0 != NULL);
   q0 = seed_bwtq(s0);
   test_assert_critical(q0 != NULL);
   test_assert(bwt_size(q0) == 3);
   test_assert(bwt_depth(q0) == 12);
   test_assert(seed_beg(s0) == 0);
   test_assert(seed_end(s0) == 12);
   
   seed_free(s0);
   seed_free(s1);
   seed_free(s2);
   seed_free(s3);
   gstack_free(stack1);
   seqread_free(read0);
   seqread_free(read1);
   index_free(index);
}

void
test_seed_helpers
(void)
{
   redirect_stderr();
   test_assert(seed_beg(NULL) == -1);
   test_assert(seed_end(NULL) == -1);
   test_assert(seed_bwtq(NULL) == NULL);
   test_assert(seed_read(NULL) == NULL);
   unredirect_stderr();
}

// Define test cases to be run (for export).
const test_case_t test_cases_seed[] = {
   {"seed/seed_new",            test_seed_new},
   {"seed/seed_stack",          test_seed_stack},
   {"seed/seed_next_mem",       test_seed_next_mem},
   {"seed/seed_mems",           test_seed_mems},
   {"seed/seed_helpers",           test_seed_helpers},
   {NULL, NULL}, // Sentinel. //
};
