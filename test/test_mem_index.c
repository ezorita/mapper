#include "unittest.h"
#include "index.h"

void
test_mem_index_build
(void)
{
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      index_t * index = index_build("examples/repeats.fa", "test_base20");
      index_free(index);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      index_t * index = index_build("examples/repeats.fa", "test_base20");
      index_free(index);
   }
   reset_alloc();
   unredirect_stderr();
}

void
test_mem_index_ann_new
(void)
{
   redirect_stderr();
   
   index_t * index = index_build("examples/repeats.fa", "test_base20");
   test_assert_critical(index != NULL);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      if (index_ann_new(25, 1, 1, index) == 0) {
         unlink("test_base20.ann.25.1");
         ann_free(index->ann[0]);
         free(index->ann);
         index->ann = NULL;
         index->ann_cnt = 0;
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 1000; i++) {
      set_alloc_failure_countdown_to(i);
      if (index_ann_new(25, 1, 1, index) == 0) {
         unlink("test_base20.ann.25.1");
         ann_free(index->ann[0]);
         free(index->ann);
         index->ann = NULL;
         index->ann_cnt = 0;
      }
   }
   reset_alloc();

   index_free(index);
   unredirect_stderr();
}

void
test_mem_index_read
(void)
{
   redirect_stderr();
   
   index_t * index = index_build("examples/repeats.fa", "test_base30");
   test_assert_critical(index != NULL);
   test_assert(index_ann_new(25, 1, 1, index) == 0);

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      index_t * index_i = index_read("test_base30");
      index_free(index_i);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      index_t * index_i = index_read("test_base30");
      index_free(index_i);
   }
   reset_alloc();

   index_free(index);
   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_mem_index[] = {
   {"mem/index/index_build",         test_mem_index_build},
   {"mem/index/index_ann_new",       test_mem_index_ann_new},
   {"mem/index/index_read",          test_mem_index_read},
   {NULL, NULL}, // Sentinel. //
};
