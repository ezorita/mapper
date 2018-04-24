#include "unittest.h"
#include "index_txt.h"

void
test_txt_new
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;

   txt = txt_new(NULL);
   test_assert(txt == NULL);

   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_length(txt) == 0);
   test_assert(txt_wildcard_count(txt) == 0);
   test_assert(txt_get_symbols(txt) == sym);
   
   test_assert(txt_append("ACTAGC", txt) == 0);
   test_assert(txt_length(txt) == 6);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 7);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_append("GGG", txt) == 0);
   test_assert(txt_length(txt) == 10);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 11);
   test_assert(txt_wildcard_count(txt) == 2);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_sym
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_sym(0, txt) == -1);

   test_assert(txt_append("ACTAGQJKCNATGCAGTCGAANNTG", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 26);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_sym(0, NULL) == -1);
   
   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(1, txt) == 1);
   test_assert(txt_sym(2, txt) == 3);
   test_assert(txt_sym(3, txt) == 0);
   test_assert(txt_sym(4, txt) == 2);
   test_assert(txt_sym(5, txt) == 4);
   test_assert(txt_sym(6, txt) == 4);
   test_assert(txt_sym(7, txt) == 4);
   test_assert(txt_sym(8, txt) == 1);
   test_assert(txt_sym(9, txt) == 4);
   test_assert(txt_sym(25, txt) == 5);
   test_assert(txt_sym(-1, txt) == -1);
   test_assert(txt_sym(132, txt) == -1);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_sym_range
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);

   test_assert_critical(txt != NULL);

   test_assert(txt_sym_range(0,1,txt) == NULL);

   test_assert(txt_append("ACTAGQJKCNATGCAGTCGAANNTG", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 26);
   test_assert(txt_wildcard_count(txt) == 1);

   uint8_t * sr;
   sr = txt_sym_range(0,100,txt);
   test_assert(sr == NULL);

   sr = txt_sym_range(-1,10,txt);
   test_assert(sr == NULL);
   
   sr = txt_sym_range(30,10,txt);
   test_assert(sr == NULL);

   sr = txt_sym_range(0,10, NULL);
   test_assert(sr == NULL);
   
   sr = txt_sym_range(0,10,txt);
   test_assert_critical(sr != NULL);
   test_assert(memcmp(sr,"\0\1\3\0\2\4\4\4\1\4",10) == 0);

   sr = txt_sym_range(10,10,txt);
   test_assert_critical(sr != NULL);
   test_assert(memcmp(sr,"\0\3\2\1\0\2\3\1\2\0",10) == 0);

   sr = txt_sym_range(15,10,txt);
   test_assert_critical(sr != NULL);
   test_assert(memcmp(sr,"\2\3\1\2\0\0\4\4\3\2",10) == 0);

   sr = txt_sym_range(20,5,txt);
   test_assert_critical(sr != NULL);
   test_assert(memcmp(sr,"\0\4\4\3\2\5",5) == 0);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_append
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append(NULL, NULL) == -1);
   test_assert(txt_append(NULL, txt) == -1);
   test_assert(txt_append("AGTACGTA", NULL) == -1);

   test_assert(txt_append("TAGCTGATCGGACGACATCN", txt) == 0);
   test_assert(txt_length(txt) == 20);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("AAAAAAAAAAAAAAAAAAAA", txt) == 0);
   test_assert(txt_length(txt) == 40);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("CCCCCCCCCCCCCCCCCCCC", txt) == 0);
   test_assert(txt_length(txt) == 60);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("GGGGGGGGGGGGGGGGGGGG", txt) == 0);
   test_assert(txt_length(txt) == 80);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("TTTTTTTTTTTTTTTTTTTT", txt) == 0);
   test_assert(txt_length(txt) == 100);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("NNNNNNNNNNNNNNNNNNNN", txt) == 0);
   test_assert(txt_length(txt) == 120);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("DEYSHDYSEDHYSEYDRRYD", txt) == 0);
   test_assert(txt_length(txt) == 140);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append("TC", txt) == 0);
   test_assert(txt_length(txt) == 142);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 143);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_sym(0,txt) == 3);
   test_assert(txt_sym(19,txt) == 4);
   test_assert(txt_sym(20,txt) == 0);
   test_assert(txt_sym(47,txt) == 1);
   test_assert(txt_sym(68,txt) == 2);
   test_assert(txt_sym(93,txt) == 3);
   test_assert(txt_sym(110,txt) == 4);
   test_assert(txt_sym(133,txt) == 4);

   uint8_t * text = txt_sym_range(0,143,txt);
   test_assert_critical(text != NULL);

   char * ref = "\3\0\2\1\3\2\0\3\1\2\2\0\1\2\0\1\0\3\1\4\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\1\1\1\1\1\1\1\1\1\1\1\1\1\1\1\1\1\1\1\1\2\2\2\2\2\2\2\2\2\2\2\2\2\2\2\2\2\2\2\2\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\3\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\4\3\1\5";
   test_assert(memcmp(text, ref, 143) == 0);

   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("0123456789",txt) == 0);
   }
   test_assert(txt_length(txt) == 10143);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_append_wildcard
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append_wildcard(NULL) == -1);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 1);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 4);
   test_assert(txt_wildcard_count(txt) == 4);

   test_assert(txt_append("TCGAC", txt) == 0);
   test_assert(txt_length(txt) == 9);
   test_assert(txt_wildcard_count(txt) == 4);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 10);
   test_assert(txt_wildcard_count(txt) == 5);

   for (int i = 0; i < 1200; i++) {
      test_assert(txt_append_wildcard(txt) == 0);
   }
   test_assert(txt_length(txt) == 1210);
   test_assert(txt_wildcard_count(txt) == 1205);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_length
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   test_assert(txt_length(NULL) == -1);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_length(txt) == 0);

   test_assert(txt_append("012345678901234567890123456789012345678901234567890123456789", txt) == 0);
   test_assert(txt_length(txt) == 60);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_length(txt) == 61);
   test_assert(txt_append("012345678901234567890123456789012345678901234567890123456789", txt) == 0);
   test_assert(txt_length(txt) == 121);

   for (int i = 0; i < 100; i++) {
      test_assert(txt_append("01234567890123456789012345678901234567890123456789", txt) == 0);
   }
   test_assert(txt_length(txt) == 5121);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_wildcard_count
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   test_assert(txt_length(NULL) == -1);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_length(txt) == 0);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_wildcard_count(txt) == 1);

   for (int i = 0; i < 1200; i++) {
      test_assert(txt_append_wildcard(txt) == 0);
   }
   test_assert(txt_wildcard_count(txt) == 1201);

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_get_symbols
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   test_assert(txt_get_symbols(NULL) == NULL);
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_get_symbols(txt) == sym);   

   txt_free(txt);
   sym_free(sym);
}


void
test_txt_file
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   test_assert(txt_get_symbols(NULL) == NULL);
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   for (int i = 0; i < 100; i++) {
      test_assert(txt_append("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTNNNNNNNNNN", txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }
   test_assert(txt_length(txt) == 5100);
   test_assert(txt_wildcard_count(txt) == 100);
   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(10, txt) == 1);
   test_assert(txt_sym(20, txt) == 2);
   test_assert(txt_sym(30, txt) == 3);
   test_assert(txt_sym(40, txt) == 4);
   test_assert(txt_sym(50, txt) == 5);
   test_assert(txt_sym(5099, txt) == 5);
   test_assert(txt_sym(5089, txt) == 4);
   test_assert(txt_sym(5079, txt) == 3);
   test_assert(txt_sym(5069, txt) == 2);
   test_assert(txt_sym(5059, txt) == 1);
   test_assert(txt_sym(5049, txt) == 0);

   test_assert(txt_file_write(NULL, txt) == -1); 
   test_assert(txt_file_write("", txt) == -1); 
   test_assert(txt_file_write("test00.txt", NULL) == -1); 

   test_assert(txt_file_write("test01.txt", txt) == 0);

   txt_free(txt);

   // Force wrong magic number.
   test_assert(txt_file_read("test01.sym", sym) == NULL);
   test_assert(txt_file_read("", sym) == NULL);
   test_assert(txt_file_read(NULL, sym) == NULL);
   test_assert(txt_file_read("test01.txt", NULL) == NULL);

   txt = txt_file_read("test01.txt", sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_length(txt) == 5100);
   test_assert(txt_wildcard_count(txt) == 100);
   test_assert(txt_get_symbols(txt) == sym);
   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(10, txt) == 1);
   test_assert(txt_sym(20, txt) == 2);
   test_assert(txt_sym(30, txt) == 3);
   test_assert(txt_sym(40, txt) == 4);
   test_assert(txt_sym(50, txt) == 5);
   test_assert(txt_sym(5099, txt) == 5);
   test_assert(txt_sym(5089, txt) == 4);
   test_assert(txt_sym(5079, txt) == 3);
   test_assert(txt_sym(5069, txt) == 2);
   test_assert(txt_sym(5059, txt) == 1);
   test_assert(txt_sym(5049, txt) == 0);

   txt_free(txt);
   sym_free(sym);
}


// Define test cases to be run (for export).
const test_case_t test_cases_index_txt[] = {
   {"index_txt/txt_new",             test_txt_new},
   {"index_txt/txt_sym",             test_txt_sym},
   {"index_txt/txt_sym_range",       test_txt_sym_range},
   {"index_txt/txt_append",          test_txt_append},
   {"index_txt/txt_append_wildcard", test_txt_append_wildcard},
   {"index_txt/txt_length",          test_txt_length},
   {"index_txt/txt_wildcard_count",  test_txt_wildcard_count},
   {"index_txt/txt_get_symbols",     test_txt_get_symbols},
   {"index_txt/txt_file",      test_txt_file},
   {NULL, NULL}, // Sentinel. //
};

