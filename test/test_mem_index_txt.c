#include "unittest.h"
#include "index_txt.h"

void
test_mem_txt_new
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   redirect_stderr();

   txt_t * txt;

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      txt = txt_new(sym);
      txt_free(txt);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      txt = txt_new(sym);
      txt_free(txt);
   }
   reset_alloc();
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_sym
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("ACTAGQJKCNATGCAGTCGAANNTG", txt) == 0);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      txt_sym(0, txt);
      txt_sym(10, txt);
      txt_sym(25, txt);
      txt_sym(-1, txt);
      txt_sym(0, NULL);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      txt_sym(0, txt);
      txt_sym(10, txt);
      txt_sym(25, txt);
      txt_sym(-1, txt);
      txt_sym(0, NULL);
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_sym_range
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("ACTAGQJKCNATGCAGTCGAANNTG", txt) == 0);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      txt_sym_range(0, 10, txt);
      txt_sym_range(10, 5, txt);
      txt_sym_range(25, 30, txt);
      txt_sym_range(-1, 5, txt);
      txt_sym_range(0, 10, NULL);
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      txt_sym_range(0, 10, txt);
      txt_sym_range(10, 5, txt);
      txt_sym_range(25, 30, txt);
      txt_sym_range(-1, 5, txt);
      txt_sym_range(0, 10, NULL);
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}


void
test_mem_txt_append
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 1000; i++) {
      txt_append("ACTAGQJKCNATGCAGTCGAANNTG", txt);
   }
   reset_alloc();

   txt_free(txt);
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      for (int j = 0; j < 1000; j++) {
         txt_append("ACTAGQJKCNATGCAGTCGAANNTG", txt);
      }
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_append_wildcard
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 10000; i++) {
      txt_append_wildcard(txt);
   }
   reset_alloc();

   txt_free(txt);
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      for (int j = 0; j < 10000; j++) {
         txt_append_wildcard(txt);
      }
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_commit_seq
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   
   char * seqname = malloc(100);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 100; i++) {
      sprintf(seqname, "%d", i);
      txt_commit_seq(seqname, txt);
   }
   reset_alloc();

   txt_free(txt);
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      for (int j = 0; j < 100; j++) {
         sprintf(seqname, "%d-%d", i,j);
         txt_commit_seq(seqname, txt);
      }
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   free(seqname);
   unredirect_stderr();
}

void
test_mem_txt_commit_rc
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   for (int i = 0; i < TXT_BUFFER_SIZE/20-2; i++) {
      txt_append("AAAACCCCGGGGTTTTNNNN", txt);
   }
   
   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 10; i++) {
      txt_commit_rc(txt);
   }
   reset_alloc();

   txt_free(txt);
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   for (int i = 0; i < TXT_BUFFER_SIZE/20-2; i++) {
      txt_append("AAAACCCCGGGGTTTTNNNN", txt);
   }

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      txt_commit_rc(txt);
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_pos_to_str
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_commit_seq("seq0", txt) == 0);
   test_assert(txt_append("TGATCGATCNTAGCT", txt) == 0);
   test_assert(txt_commit_seq("seq1", txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 10; i++) {
      free(txt_pos_to_str(9,txt));
      free(txt_pos_to_str(12,txt));
      free(txt_pos_to_str(27,txt));
      free(txt_pos_to_str(52,txt));
   }

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      free(txt_pos_to_str(9,txt));
      free(txt_pos_to_str(12,txt));
      free(txt_pos_to_str(27,txt));
      free(txt_pos_to_str(52,txt));
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_str_to_pos
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_commit_seq("seq0", txt) == 0);
   test_assert(txt_append("TGATCGATCNTAGCT", txt) == 0);
   test_assert(txt_commit_seq("seq1", txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 10; i++) {
      txt_str_to_pos("seq0:10:+", txt);
      txt_str_to_pos("seq1:2:+", txt);
      txt_str_to_pos("seq1:15:-", txt);
      txt_str_to_pos("seq0:1:-", txt);
      txt_str_to_pos("seq1:15:+", txt);
   }

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      txt_str_to_pos("seq0:10:+", txt);
      txt_str_to_pos("seq1:2:+", txt);
      txt_str_to_pos("seq1:15:-", txt);
      txt_str_to_pos("seq0:1:-", txt);
      txt_str_to_pos("seq1:15:+", txt);
   }
   reset_alloc();
   
   txt_free(txt);
   sym_free(sym);
   unredirect_stderr();
}

void
test_mem_txt_file
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt_o, * txt_i;
   txt_o = txt_new(sym);
   test_assert_critical(txt_o != NULL);

   test_assert(txt_append("AACCGGTTNN", txt_o) == 0);
   test_assert(txt_commit_seq("seq0", txt_o) == 0);
   test_assert(txt_append("TGATCGATCNTAGCT", txt_o) == 0);
   test_assert(txt_commit_seq("seq1", txt_o) == 0);
   test_assert(txt_commit_rc(txt_o) == 0);

   redirect_stderr();
   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.1);
   for (int i = 0; i < 10; i++) {
      txt_file_write("test30.txt", txt_o);
      txt_i = txt_file_read("test30.txt", sym);
      txt_free(txt_i);
   }

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 10; i++) {
      set_alloc_failure_countdown_to(i);
      txt_file_write("test30.txt", txt_o);
      txt_i = txt_file_read("test30.txt", sym);
      txt_free(txt_i);
   }
   reset_alloc();
   
   txt_free(txt_o);
   sym_free(sym);
   unredirect_stderr();
}



// Define test cases to be run (for export).
const test_case_t test_mem_index_txt[] = {
   {"mem/index_txt/txt_new",             test_mem_txt_new},
   {"mem/index_txt/txt_sym",             test_mem_txt_sym},
   {"mem/index_txt/txt_sym_range",       test_mem_txt_sym_range},
   {"mem/index_txt/txt_append",          test_mem_txt_append},
   {"mem/index_txt/txt_append_wildcard", test_mem_txt_append_wildcard},
   {"mem/index_txt/txt_commit_seq",      test_mem_txt_commit_seq},
   {"mem/index_txt/txt_commit_rc",       test_mem_txt_commit_rc},
   {"mem/index_txt/txt_seq_pos_to_str",  test_mem_txt_pos_to_str},
   {"mem/index_txt/txt_seq_str_to_pos",  test_mem_txt_str_to_pos},
   {"mem/index_txt/txt_file",            test_mem_txt_file},
   {NULL, NULL}, // Sentinel. //
};
