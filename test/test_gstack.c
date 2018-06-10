#include "unittest.h"
#include "gstack.h"

void
test_gstack_new
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Invalid parameters.
   test_assert(gstack_new(0x0FFFFFFFFFFFFFFF, free) == NULL);

   // Create some gstacks.
   gstack = gstack_new(10,free);
   test_assert(gstack != NULL);
   gstack_free(gstack);
   gstack = gstack_new(100,free);
   test_assert(gstack != NULL);
   gstack_free(gstack);
   gstack = gstack_new(1000,free);
   test_assert(gstack != NULL);
   gstack_free(gstack);
   gstack = gstack_new(10000,free);
   test_assert(gstack != NULL);
   gstack_free(gstack);

   unredirect_stderr();
}

void
test_gstack_push
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);
   
   int64_t * value = NULL;
   
   // Invalid parameters.
   test_assert(gstack_push(NULL, gstack) == -1);
   
   // Fill gstack.
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   test_assert(gstack_push(value, NULL) == -1);
   
   *value = 1;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 49;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 193;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 257;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 3324;
   test_assert(gstack_push(value, gstack) == 0);

   gstack_free(gstack);
   unredirect_stderr();
}

void
test_gstack_push_array
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   int64_t fibo[13] = {1,2,3,5,8,13,21,34,55,89,144,233,377};

   int64_t ** values = malloc(sizeof(void *)*13);
   test_assert_critical(values != NULL);
   for (int i = 0; i < 13; i++) {
      values[i] = malloc(sizeof(int64_t));
      test_assert_critical(values[i] != NULL);
      *values[i] = fibo[i];
   }
   
   // Invalid parameters.
   test_assert(gstack_push_array((void **)values, 1, NULL) == -1);
   test_assert(gstack_push_array(NULL, 1, gstack) == -1);
   test_assert(gstack_push_array((void **)values, 0, gstack) == 0);
   
   // Fill gstack.
   test_assert(gstack_push_array((void **)values, 10, gstack) == 0);

   // Check values.
   int64_t * get_val;
   get_val = gstack_get(0, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 1);

   get_val = gstack_get(3, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 5);

   get_val = gstack_get(6, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 21);

   get_val = gstack_get(9, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 89);

   // Fill further.
   test_assert(gstack_push_array((void **)values + 10, 3, gstack) == 0);

   // Check values.
   get_val = gstack_get(10, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 144);

   get_val = gstack_get(11, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 233);

   get_val = gstack_get(12, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 377);

   gstack_free(gstack);
   free(values);
   unredirect_stderr();
}

void
test_gstack_transfer_all
(void)
{
   redirect_stderr();
   gstack_t * stack_src = NULL;
   gstack_t * stack_dst = NULL;
   
   // Create a gstack.
   stack_src = gstack_new(1, free);
   test_assert_critical(stack_src != NULL);
   stack_dst = gstack_new(1, free);
   test_assert_critical(stack_dst != NULL);

   int64_t fibo[13] = {1,2,3,5,8,13,21,34,55,89,144,233,377};

   int64_t ** values = malloc(sizeof(void *)*13);
   test_assert_critical(values != NULL);
   for (int i = 0; i < 13; i++) {
      values[i] = malloc(sizeof(int64_t));
      test_assert_critical(values[i] != NULL);
      *values[i] = fibo[i];
   }
 
   // Fill gstack.
   test_assert(gstack_push_array((void **)values, 10, stack_src) == 0);
   test_assert(gstack_push(values[10], stack_dst) == 0);

   // Check contents.
   test_assert(gstack_num_elm(stack_src) == 10);
   test_assert(gstack_max_elm(stack_src) == 16);
   test_assert(gstack_num_elm(stack_dst) == 1);
   test_assert(gstack_max_elm(stack_dst) == 1);

   // Transfer all.
   test_assert(gstack_transfer_all(stack_dst, stack_src) == 0);

   // Check final contents.
   test_assert(gstack_num_elm(stack_dst) == 11);
   test_assert(*(int64_t *)gstack_get(0, stack_dst) == 144);
   test_assert(*(int64_t *)gstack_get(1, stack_dst) == 1);
   test_assert(*(int64_t *)gstack_get(10, stack_dst) == 89);
   test_assert(gstack_num_elm(stack_src) == 0);
   test_assert(gstack_max_elm(stack_src) == 1);

   // Push remaining values.
   test_assert(gstack_push_array((void **)values+11, 2, stack_dst) == 0);
   test_assert(gstack_num_elm(stack_dst) == 13);
   test_assert(*(int64_t *)gstack_get(11, stack_dst) == 233);
   test_assert(*(int64_t *)gstack_get(12, stack_dst) == 377);

   gstack_free(stack_src);
   gstack_free(stack_dst);
   free(values);
}


void
test_gstack_pop
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(100, free);
   test_assert_critical(gstack != NULL);

   // Invalid parameters.
   test_assert(gstack_pop(NULL) == NULL);
   
   // Pop from empty gstack.
   test_assert(gstack_pop(gstack) == NULL);

   // Fill gstack.
   int64_t * value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 10;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 49;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 193;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 257;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 3324;
   test_assert(gstack_push(value, gstack) == 0);

   // Test pop.
   int64_t * pop_val;
   pop_val = gstack_pop(gstack);
   test_assert_critical(pop_val != NULL);
   test_assert(*pop_val == 3324);
   free(pop_val);

   pop_val = gstack_pop(gstack);
   test_assert_critical(pop_val != NULL);
   test_assert(*pop_val == 257);
   free(pop_val);

   pop_val = gstack_pop(gstack);
   test_assert_critical(pop_val != NULL);
   test_assert(*pop_val == 193);
   free(pop_val);

   pop_val = gstack_pop(gstack);
   test_assert_critical(pop_val != NULL);
   test_assert(*pop_val == 49);
   free(pop_val);

   pop_val = gstack_pop(gstack);
   test_assert_critical(pop_val != NULL);
   test_assert(*pop_val == 10);
   free(pop_val);

   pop_val = gstack_pop(gstack);
   test_assert(pop_val == NULL);

   gstack_free(gstack);
   unredirect_stderr();   
}

void
test_gstack_get
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(100, free);
   test_assert_critical(gstack != NULL);

   // Invalid parameters.
   test_assert(gstack_get(0, NULL) == NULL);
   test_assert(gstack_get(0, gstack) == NULL);

   // Fill gstack.
   int64_t * value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 10;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 49;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 193;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 257;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 3324;
   test_assert(gstack_push(value, gstack) == 0);

   // Get stack elements.
   int64_t * get_val;
   get_val = gstack_get(0, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 10);
    
   get_val = gstack_get(1, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 49);

   get_val = gstack_get(2, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 193);

   get_val = gstack_get(3, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 257);

   get_val = gstack_get(4, gstack);
   test_assert_critical(get_val != NULL);
   test_assert(*get_val == 3324);

   gstack_free(gstack);
   unredirect_stderr();
}

void
test_gstack_popr
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(100, free);
   test_assert_critical(gstack != NULL);

   // Fill gstack.
   int64_t * value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 10;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 49;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 193;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 257;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 3324;
   test_assert(gstack_push(value, gstack) == 0);

   // popr
   test_assert(gstack_popr(NULL) == -1);
   test_assert(gstack_num_elm(gstack) == 5);
   test_assert(gstack_popr(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 4);
   test_assert(gstack_popr(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 3);
   test_assert(gstack_popr(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 2);
   test_assert(gstack_popr(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 1);
   test_assert(gstack_popr(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 0);
   test_assert(gstack_popr(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 0);

   gstack_free(gstack);
   unredirect_stderr();
}

void
test_gstack_clear
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(100, free);
   test_assert_critical(gstack != NULL);

   // Fill gstack.
   int64_t * value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 10;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 49;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 193;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 257;
   test_assert(gstack_push(value, gstack) == 0);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 3324;
   test_assert(gstack_push(value, gstack) == 0);

   // popr
   test_assert(gstack_clear(NULL) == -1);
   test_assert(gstack_num_elm(gstack) == 5);
   test_assert(gstack_clear(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 0);
   test_assert(gstack_clear(gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 0);

   gstack_free(gstack);
   unredirect_stderr();
}

void
test_gstack_helpers
(void)
{
   redirect_stderr();
   gstack_t * gstack = NULL;
   
   // Create a gstack.
   gstack = gstack_new(1, free);
   test_assert_critical(gstack != NULL);

   // Fill gstack.
   test_assert(gstack_num_elm(gstack) == 0);
   test_assert(gstack_max_elm(gstack) == 1);

   int64_t * value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 10;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 1);
   test_assert(gstack_max_elm(gstack) == 1);

   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 49;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 2);
   test_assert(gstack_max_elm(gstack) == 2);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 193;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 3);
   test_assert(gstack_max_elm(gstack) == 4);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 257;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 4);
   test_assert(gstack_max_elm(gstack) == 4);
   
   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 3324;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 5);
   test_assert(gstack_max_elm(gstack) == 8);

   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 24194;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 6);
   test_assert(gstack_max_elm(gstack) == 8);

   value = malloc(sizeof(int64_t));
   test_assert_critical(value != NULL);
   *value = 38113;
   test_assert(gstack_push(value, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 7);
   test_assert(gstack_max_elm(gstack) == 8);

   // Push array.
   int64_t fibo[13] = {1,2,3,5,8,13,21,34,55,89,144,233,377};

   int64_t ** values = malloc(sizeof(void *)*13);
   test_assert_critical(values != NULL);
   for (int i = 0; i < 13; i++) {
      values[i] = malloc(sizeof(int64_t));
      test_assert_critical(values[i] != NULL);
      *values[i] = fibo[i];
   }

   test_assert(gstack_push_array((void **)values, 13, gstack) == 0);
   test_assert(gstack_num_elm(gstack) == 20);
   test_assert(gstack_max_elm(gstack) == 32);

   gstack_free(gstack);
   free(values);
   unredirect_stderr();
}

// Define test cases to be run (for export).
const test_case_t test_cases_gstack[] = {
   {"gstack/gstack_new",           test_gstack_new},
   {"gstack/gstack_push",          test_gstack_push},
   {"gstack/gstack_push_array",    test_gstack_push_array},
   {"gstack/gstack_transfer_all",    test_gstack_transfer_all},
   {"gstack/gstack_pop",           test_gstack_pop},
   {"gstack/gstack_get",           test_gstack_get},
   {"gstack/gstack_popr",          test_gstack_popr},
   {"gstack/gstack_clear",         test_gstack_clear},
   {"gstack/gstack_helpers",       test_gstack_helpers},
   {NULL, NULL}, // Sentinel. //
};
