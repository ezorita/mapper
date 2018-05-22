#include "unittest.h"
#include "index_txt.h"

void
test_txt_new
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   redirect_stderr();
   
   txt_t * txt;

   txt = txt_new(NULL);
   test_assert(txt == NULL);

   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_length(txt) == 0);
   test_assert(txt_wildcard_count(txt) == 0);
   test_assert(txt_seq_count(txt) == 0);
   test_assert(txt_get_symbols(txt) == sym);
   
   test_assert(txt_append("ACTAGC", txt) == 0);
   test_assert(txt_length(txt) == 6);
   test_assert(txt_seq_count(txt) == 0);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_commit_seq("seq1",txt) == 0);
   test_assert(txt_length(txt) == 7);
   test_assert(txt_seq_count(txt) == 1);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_append("GGG", txt) == 0);
   test_assert(txt_length(txt) == 10);
   test_assert(txt_seq_count(txt) == 1);
   test_assert(txt_wildcard_count(txt) == 1);

   test_assert(txt_commit_seq("seq2",txt) == 0);
   test_assert(txt_length(txt) == 11);
   test_assert(txt_seq_count(txt) == 2);
   test_assert(txt_wildcard_count(txt) == 2);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


void
test_txt_sym
(void)
{
   redirect_stderr();

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

   unredirect_stderr();
}


void
test_txt_sym_range
(void)
{
   redirect_stderr();

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

   unredirect_stderr();
}


void
test_txt_append
(void)
{
   redirect_stderr();

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

   unredirect_stderr();
}


void
test_txt_append_wildcard
(void)
{
   redirect_stderr();

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

   unredirect_stderr();
}


void
test_txt_commit_seq
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_length(txt) == 10);

   test_assert(txt_commit_seq(NULL,txt) == -1);
   test_assert(txt_commit_seq("",txt) == -1);
   test_assert(txt_commit_seq("seq0",NULL) == -1);

   test_assert(txt_commit_seq("seq0",txt) == 0);

   test_assert(txt_wildcard_count(txt) == 1);
   test_assert(txt_length(txt) == 11);
   test_assert(txt_seq_count(txt) == 1);
   test_assert(txt_seq_start(0,txt) == 0);
   test_assert(txt_seq_length(0,txt) == 11);
   test_assert(strcmp(txt_seq_name(0,txt), "seq0") == 0);

   test_assert(txt_append("TGATCGATCNTAGCT", txt) == 0);
   test_assert(txt_length(txt) == 26);
   // Seq name must be unique
   test_assert(txt_commit_seq("seq0",txt) == -1);
   test_assert(txt_commit_seq("seq1",txt) == 0);

   test_assert(txt_wildcard_count(txt) == 2);
   test_assert(txt_length(txt) == 27);
   test_assert(txt_seq_count(txt) == 2);
   test_assert(txt_seq_start(1,txt) == 11);
   test_assert(txt_seq_length(1,txt) == 16);
   test_assert(strcmp(txt_seq_name(1,txt), "seq1") == 0);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


void
test_txt_commit_rc
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_length(txt) == 10);
   test_assert(txt_wildcard_count(txt) == 0);

   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(1, txt) == 0);
   test_assert(txt_sym(2, txt) == 1);
   test_assert(txt_sym(3, txt) == 1);
   test_assert(txt_sym(4, txt) == 2);
   test_assert(txt_sym(5, txt) == 2);
   test_assert(txt_sym(6, txt) == 3);
   test_assert(txt_sym(7, txt) == 3);
   test_assert(txt_sym(8, txt) == 4);
   test_assert(txt_sym(9, txt) == 4);

   test_assert(txt_commit_rc(txt) == 0);
   test_assert(txt_length(txt) == 22);
   test_assert(txt_wildcard_count(txt) == 2);

   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(1, txt) == 0);
   test_assert(txt_sym(2, txt) == 1);
   test_assert(txt_sym(3, txt) == 1);
   test_assert(txt_sym(4, txt) == 2);
   test_assert(txt_sym(5, txt) == 2);
   test_assert(txt_sym(6, txt) == 3);
   test_assert(txt_sym(7, txt) == 3);
   test_assert(txt_sym(8, txt) == 4);
   test_assert(txt_sym(9, txt) == 4);
   test_assert(txt_sym(10, txt) == 5);
   test_assert(txt_sym(11, txt) == 4);
   test_assert(txt_sym(12, txt) == 4);
   test_assert(txt_sym(13, txt) == 0);
   test_assert(txt_sym(14, txt) == 0);
   test_assert(txt_sym(15, txt) == 1);
   test_assert(txt_sym(16, txt) == 1);
   test_assert(txt_sym(17, txt) == 2);
   test_assert(txt_sym(18, txt) == 2);
   test_assert(txt_sym(19, txt) == 3);
   test_assert(txt_sym(20, txt) == 3);
   test_assert(txt_sym(21, txt) == 5);

   txt_free(txt);

   // Another example with more wildcards: AC$GT$NN$
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AC", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("GT", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("NN", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   test_assert(txt_length(txt) == 9);
   test_assert(txt_wildcard_count(txt) == 3);

   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(1, txt) == 1);
   test_assert(txt_sym(2, txt) == 5);
   test_assert(txt_sym(3, txt) == 2);
   test_assert(txt_sym(4, txt) == 3);
   test_assert(txt_sym(5, txt) == 5);
   test_assert(txt_sym(6, txt) == 4);
   test_assert(txt_sym(7, txt) == 4);
   test_assert(txt_sym(8, txt) == 5);

   test_assert(txt_commit_rc(txt) == 0);
   test_assert(txt_length(txt) == 18);
   test_assert(txt_wildcard_count(txt) == 6);

   test_assert(txt_sym(0, txt) == 0);
   test_assert(txt_sym(1, txt) == 1);
   test_assert(txt_sym(2, txt) == 5);
   test_assert(txt_sym(3, txt) == 2);
   test_assert(txt_sym(4, txt) == 3);
   test_assert(txt_sym(5, txt) == 5);
   test_assert(txt_sym(6, txt) == 4);
   test_assert(txt_sym(7, txt) == 4);
   test_assert(txt_sym(8, txt) == 5);
   test_assert(txt_sym(9, txt) == 4);
   test_assert(txt_sym(10, txt) == 4);
   test_assert(txt_sym(11, txt) == 5);
   test_assert(txt_sym(12, txt) == 0);
   test_assert(txt_sym(13, txt) == 1);
   test_assert(txt_sym(14, txt) == 5);
   test_assert(txt_sym(15, txt) == 2);
   test_assert(txt_sym(16, txt) == 3);
   test_assert(txt_sym(17, txt) == 5);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


void
test_txt_length
(void)
{
   redirect_stderr();

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

   unredirect_stderr();
}


void
test_txt_wildcard_count
(void)
{
   redirect_stderr();

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

   unredirect_stderr();
}


void
test_txt_get_symbols
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   test_assert(txt_get_symbols(NULL) == NULL);
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   test_assert(txt_get_symbols(txt) == sym);   

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


void
test_txt_seq_helpers
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_append("TGATCGATCNTAGCT", txt) == 0);
   test_assert(txt_commit_seq("seq1",txt) == 0);

   test_assert(txt_seq_count(NULL) == -1);
   test_assert(txt_seq_count(txt) == 2);

   test_assert(txt_seq_start(0,NULL) == -1);
   test_assert(txt_seq_start(-1,txt) == -1);
   test_assert(txt_seq_start(2,txt) == -1);
   test_assert(txt_seq_start(0,txt) == 0);
   test_assert(txt_seq_start(1,txt) == 11);

   test_assert(txt_seq_length(0,NULL) == -1);
   test_assert(txt_seq_length(-1,txt) == -1);
   test_assert(txt_seq_length(2,txt) == -1);
   test_assert(txt_seq_length(0,txt) == 11);
   test_assert(txt_seq_length(1,txt) == 16);

   test_assert(txt_seq_name(0,NULL) == NULL);
   test_assert(txt_seq_name(-1,txt) == NULL);
   test_assert(txt_seq_name(2,txt) == NULL);
   test_assert(strcmp(txt_seq_name(0,txt), "seq0") == 0);
   test_assert(strcmp(txt_seq_name(1,txt), "seq1") == 0);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


void
test_txt_pos_to_str
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_append("TGATCGATCNTAGCT", txt) == 0);
   test_assert(txt_commit_seq("seq1",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   //0        90        90     6   0         0         0
   //AACCGGTTNN$TGATCGATCNTAGCT$TGATCGATCNTAGCT$AACCGGTTNN$

   test_assert(txt_pos_to_str(-1,txt) == NULL);
   test_assert(txt_pos_to_str(0,NULL) == NULL);
   test_assert(strcmp(txt_pos_to_str(9,txt), "seq0:10:+") == 0);
   test_assert(strcmp(txt_pos_to_str(12,txt), "seq1:2:+") == 0);
   test_assert(strcmp(txt_pos_to_str(27,txt), "seq1:15:-") == 0);
   test_assert(strcmp(txt_pos_to_str(52,txt), "seq0:1:-") == 0);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}

void
test_txt_str_to_pos
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   test_assert(txt_append("AACCGGTTNN", txt) == 0);
   test_assert(txt_commit_seq("seq0",txt) == 0);
   test_assert(txt_append("TGATCGATCNTAGCT", txt) == 0);
   test_assert(txt_commit_seq("seq1",txt) == 0);
   test_assert(txt_commit_rc(txt) == 0);

   //0        90        90     6   0         0         0
   //AACCGGTTNN$TGATCGATCNTAGCT$TGATCGATCNTAGCT$AACCGGTTNN$

   test_assert(txt_str_to_pos(NULL,txt) == -1);
   test_assert(txt_str_to_pos("seq0:8:+",NULL) == -1);
   test_assert(txt_str_to_pos("seqNULL:10:+",txt) == -1);
   test_assert(txt_str_to_pos("seqNULL:+",txt) == -1);
   test_assert(txt_str_to_pos("seq0,10,+",txt) == -1);
   test_assert(txt_str_to_pos("seq0:-1:+",txt) == -1);
   test_assert(txt_str_to_pos("seq0:100:+",txt) == -1);
   test_assert(txt_str_to_pos("seq0:0:+",txt) == -1);

   test_assert(txt_str_to_pos("seq0:10:+",txt) == 9);
   test_assert(txt_str_to_pos("seq1:2:+",txt) == 12);
   test_assert(txt_str_to_pos("seq1:15:-",txt) == 27);
   test_assert(txt_str_to_pos("seq0:1:-",txt) == 52);
   test_assert(txt_str_to_pos("seq1:15:+",txt) == 25);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}



void
test_txt_file
(void)
{
   redirect_stderr();

   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);
   
   test_assert(txt_get_symbols(NULL) == NULL);
   txt_t * txt;
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);

   char * sn = malloc(10);
   for (int i = 0; i < 100; i++) {
      sprintf(sn, "seq%d", i);
      test_assert(txt_append("AAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTNNNNNNNNNN", txt) == 0);
      test_assert(txt_commit_seq(sn,txt) == 0);
   }
   free(sn);

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
   
   test_assert(strcmp(txt_seq_name(0,txt), "seq0") == 0);
   test_assert(strcmp(txt_seq_name(71,txt), "seq71") == 0);
   test_assert(strcmp(txt_seq_name(99,txt), "seq99") == 0);
   test_assert(txt_seq_name(100,txt) == NULL);

   test_assert(txt_seq_length(0,txt) == 51);
   test_assert(txt_seq_length(34,txt) == 51);
   test_assert(txt_seq_length(99,txt) == 51);
   test_assert(txt_seq_length(100,txt) == -1);

   test_assert(txt_seq_start(0,txt) == 0);
   test_assert(txt_seq_start(1,txt) == 51);
   test_assert(txt_seq_start(26,txt) == 1326);
   test_assert(txt_seq_start(99,txt) == 5049);
   test_assert(txt_seq_start(100,txt) == -1);

   test_assert(strcmp(txt_pos_to_str(0,txt), "seq0:1:+") == 0);

   txt_free(txt);
   sym_free(sym);

   unredirect_stderr();
}


// Define test cases to be run (for export).
const test_case_t test_cases_index_txt[] = {
   {"index_txt/txt_new",             test_txt_new},
   {"index_txt/txt_sym",             test_txt_sym},
   {"index_txt/txt_sym_range",       test_txt_sym_range},
   {"index_txt/txt_append",          test_txt_append},
   {"index_txt/txt_append_wildcard", test_txt_append_wildcard},
   {"index_txt/txt_commit_seq",      test_txt_commit_seq},
   {"index_txt/txt_commit_rc",       test_txt_commit_rc},
   {"index_txt/txt_length",          test_txt_length},
   {"index_txt/txt_wildcard_count",  test_txt_wildcard_count},
   {"index_txt/txt_get_symbols",     test_txt_get_symbols},
   {"index_txt/txt_seq_helpers",     test_txt_seq_helpers},
   {"index_txt/txt_seq_pos_to_str",  test_txt_pos_to_str},
   {"index_txt/txt_seq_str_to_pos",  test_txt_str_to_pos},
   {"index_txt/txt_file",      test_txt_file},
   {NULL, NULL}, // Sentinel. //
};

