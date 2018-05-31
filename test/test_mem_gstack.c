#include "unittest.h"
#include "gstack.h"

void
test_mem_gstack_new
(void)
{
   redirect_stderr();

   gstack_t * gstack = NULL;

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      gstack = gstack_new(100, NULL);
      gstack_free(gstack);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      gstack = gstack_new(100, NULL);
      gstack_free(gstack);
   }
   reset_alloc();

   unredirect_stderr();
}

void
test_mem_gstack_push
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   int64_t * value;

   // Set alloc failure rate to 0.1.
   for (int i = 0; i < 1000; i++) {
      reset_alloc();
      value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      set_alloc_failure_rate_to(0.1);
      if (gstack_push(value, gstack) == -1) {
	 free(value);
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      reset_alloc();
      value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      set_alloc_failure_countdown_to(i);
      if (gstack_push(value, gstack) == -1) {
	 free(value);
      }
   }
   reset_alloc();

   gstack_free(gstack);

   unredirect_stderr();
}   

void
test_mem_gstack_push_array
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   int64_t ** values = malloc(10*sizeof(int64_t*));
   test_assert_critical(values != NULL);

   // Set alloc failure rate to 0.1.
   for (int i = 0; i < 100; i++) {
      reset_alloc();
      for (int i = 0 ; i < 10; i++) {
	 values[i] = malloc(sizeof(int64_t));
	 test_assert_critical(values[i] != NULL);
      }
      set_alloc_failure_rate_to(0.1);
      gstack_push_array((void **)values, 10, gstack);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      reset_alloc();
      for (int i = 0 ; i < 10; i++) {
	 values[i] = malloc(sizeof(int64_t));
	 test_assert_critical(values[i] != NULL);
      }
      set_alloc_failure_countdown_to(i);
      gstack_push_array((void **)values, 10, gstack);
   }
   reset_alloc();

   gstack_free(gstack);
   free(values);
   unredirect_stderr();
}

void
test_mem_gstack_pop
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   for (int i = 0; i < 900; i++) {
      int64_t * value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      gstack_push(value, gstack);	 
   }

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      int64_t * ptr = gstack_pop(gstack);
      free(ptr);
   }
   reset_alloc();

   for (int i = 0; i < 150; i++) {
      int64_t * value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      gstack_push(value, gstack);	 
   }

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      int64_t * ptr = gstack_pop(gstack);
      free(ptr);
   }
   reset_alloc();

   gstack_free(gstack);

   unredirect_stderr();
}

void
test_mem_gstack_get
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   for (int i = 0; i < 900; i++) {
      int64_t * value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      gstack_push(value, gstack);	 
   }
 
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      gstack_get(i,gstack);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      gstack_get(i,gstack);
   }
   reset_alloc();

   gstack_free(gstack);

   unredirect_stderr();
}

void
test_mem_gstack_popr
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
  // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   for (int i = 0; i < 1200; i++) {
      int64_t * value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      gstack_push(value, gstack);	 
   }

    // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      gstack_popr(gstack);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      gstack_popr(gstack);
   }
   reset_alloc();

   gstack_free(gstack);

   unredirect_stderr();
}

void
test_mem_gstack_clear
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   // fill stack
   for (int i = 0; i < 900; i++) {
      int64_t * value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      gstack_push(value, gstack);	 
   }

    // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      gstack_clear(gstack);
   }
   reset_alloc();

   //refill
   for (int i = 0; i < 900; i++) {
      int64_t * value = malloc(sizeof(int64_t));
      test_assert_critical(value != NULL);
      gstack_push(value, gstack);	 
   }

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 200; i++) {
      set_alloc_failure_countdown_to(i);
      gstack_clear(gstack);
   }
   reset_alloc();

   gstack_free(gstack);

   unredirect_stderr();
}



// Define test cases to be run (for export).
const test_case_t test_mem_gstack[] = {
   {"mem/gstack/gstack_new",           test_mem_gstack_new},
   {"mem/gstack/gstack_push",          test_mem_gstack_push},
   {"mem/gstack/gstack_push_array",    test_mem_gstack_push_array},
   {"mem/gstack/gstack_pop",           test_mem_gstack_pop},
   {"mem/gstack/gstack_get",           test_mem_gstack_get},
   {"mem/gstack/gstack_popr",          test_mem_gstack_popr},
   {"mem/gstack/gstack_clear",         test_mem_gstack_clear},
   {NULL, NULL}, // Sentinel. //
};
