#include "unittest.h"
#include "seed.h"

void
test_mem_seed_new
(void)
{   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   
   seqread_t * read = seqread_new("one-exact", "ATCGATATCAGCCACTACGAGACAA", NULL);
   test_assert_critical(read != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   for (int i = 0; i < 1000; i++) {
      bwtquery_t * q = bwt_new_query(index->bwt);
      test_assert_critical(q != NULL);
      set_alloc_failure_rate_to(0.1);
      seed_t * seed = seed_new(0, 25, q, read);
      reset_alloc();
      if (seed != NULL) {
	 seed_free(seed);
      } else {
	 free(q);
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      bwtquery_t * q = bwt_new_query(index->bwt);
      test_assert_critical(q != NULL);
      set_alloc_failure_countdown_to(i);
      seed_t * seed = seed_new(0, 25, q, read);
      reset_alloc();
      if (seed != NULL) {
	 seed_free(seed);
      } else {
	 free(q);
      }
   }
   reset_alloc();
   unredirect_stderr();

   seqread_free(read);
   index_free(index);
}

void
test_mem_seed_stack
(void)
{
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      gstack_t * stack = seed_stack(100);
      gstack_free(stack);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      gstack_t * stack = seed_stack(100);
      gstack_free(stack);
   }
   reset_alloc();
   unredirect_stderr();
}
void
test_mem_seed_get_pop
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   
   seqread_t * read = seqread_new("one-exact", "ATCGATATCAGCCACTACGAGACAA", NULL);
   test_assert_critical(read != NULL);

   gstack_t * stack = seed_stack(100);
   test_assert_critical(stack != NULL);


   for (int i = 0; i < 200; i++) {
      bwtquery_t * q = bwt_new_query(index->bwt);
      test_assert_critical(q != NULL);
      seed_t * seed = seed_new(0, 25, q, read);
      test_assert_critical(seed != NULL);
      test_assert(gstack_push(seed, stack) == 0);
   }


   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      seed_t * s = seed_get(0,stack);
      s = seed_pop(stack);
      if (s != NULL) {
	 seed_free(s);
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i < 100; i++) {
      set_alloc_failure_countdown_to(i);
      seed_t * s = seed_get(0,stack);
      s = seed_pop(stack);
      if (s != NULL) {
	 seed_free(s);
      }
   }
   
   reset_alloc();
   unredirect_stderr();

   gstack_free(stack);
   seqread_free(read);
   index_free(index);
}

void
test_mem_seed_next_mem
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   seqread_t * read = seqread_new("one-exact", "ATCGATATCAGCTACTACGAGACAA", NULL);
   test_assert_critical(read != NULL);

   seed_t * s0 = seed_next_mem(NULL, read, index);
   test_assert_critical(s0 != NULL);
   seed_t * s1 = seed_next_mem(s0, read, index);
   test_assert_critical(s1 != NULL);
   seed_t * s2 = seed_next_mem(s1, read, index);
   test_assert_critical(s2 != NULL);
   
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      seed_t * s = seed_next_mem(NULL, read, index);
      if (s) seed_free(s);
      s = seed_next_mem(s0, read, index);
      if (s) seed_free(s);
      s = seed_next_mem(s1, read, index);
      if (s) seed_free(s);
      s = seed_next_mem(s2, read, index);
      if (s) seed_free(s);
   }
   reset_alloc();

      // Set alloc countdown 0->10.
   for (int i = 0; i < 200; i++) {
      set_alloc_failure_countdown_to(i);
      seed_t * s = seed_next_mem(NULL, read, index);
      if (s) seed_free(s);
      s = seed_next_mem(s0, read, index);
      if (s) seed_free(s);
      s = seed_next_mem(s1, read, index);
      if (s) seed_free(s);
      s = seed_next_mem(s2, read, index);
      if (s) seed_free(s);
   }
   reset_alloc();
   unredirect_stderr();

   seed_free(s0);
   seed_free(s1);
   seed_free(s2);
   seqread_free(read);
   index_free(index);
}

void
test_mem_seed_mems
(void)
{
   index_t * index = index_read("test_base00");
   test_assert_critical(index != NULL);
   seqread_t * read = seqread_new("one-exact", "ATCGATATCAGCTACTACGAGACAA", NULL);
   test_assert_critical(read != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      gstack_t * stack = seed_mems(read, index);
      if (stack != NULL) {
	 gstack_free(stack);
      }
   }
   reset_alloc();

      // Set alloc countdown 0->10.
   for (int i = 0; i < 200; i++) {
      set_alloc_failure_countdown_to(i);
            gstack_t * stack = seed_mems(read, index);
      if (stack != NULL) {
	 gstack_free(stack);
      }
   }
   reset_alloc();

   unredirect_stderr();
   seqread_free(read);
   index_free(index);
}



// Define test cases to be run (for export).
const test_case_t test_mem_seed[] = {
   {"mem/seed/seed_new",            test_mem_seed_new},
   {"mem/seed/seed_stack",          test_mem_seed_stack},
   {"mem/seed/seed_next_mem",       test_mem_seed_next_mem},
   {"mem/seed/seed_mems",           test_mem_seed_mems},
   {NULL, NULL}, // Sentinel. //
};
