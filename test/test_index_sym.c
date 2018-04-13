#include "unittest.h"
#include "index_sym.h"

// EDUARD ELS TESTS NOMES HAURIEN DE TESTEJAR LA INTERFICIE!!!!
// PENSA QUE SI CANVIES ALGUNES MINUCIES INTERNES DELS TEST HAS DE CANVIAR-LOS TOTS.
// LA IDEA ES QUE EL CONTINGUT SIGUI INDIFERENT, PERO QUE ELS RESULTATS SEMPRE SIGUIN ELS ESPERATS! COMPRENDE?

void
test_sym_new
(void)
{
   char * alph00[] = {NULL};
   char * alph01[] = {"Aa",NULL};
   char * alph02[] = {"A","A",NULL};
   char * alph03[] = {"Aa","Bbk","Cjc","Da",NULL};
   char * alph04[] = {"Aa","Bb","Cc","Dd",NULL};
   char * alph05[] = {"Aa","Bb","Cc","Dd","Ee","F","Ggf","HhIi",NULL};
   char * comp01[] = {"AB", "BA", "CD", "KQ", NULL};
   char * comp02[] = {"AB", "bA", "CD", "DC", NULL};
   char * comp03[] = {"aB", "BA", "CD", "DC", NULL};
   char * comp04[] = {"AB", "BA", "CD", "dc", NULL};
   char * comp05[] = {"AB", "BA", "CD", "DC", NULL};
   // NULL alphabet.
   sym_t * sym0 = sym_new(NULL, NULL, 0);
   test_assert(sym0 == NULL);

   // Empty alphabet.
   sym0 = sym_new(alph00, NULL, 0);
   test_assert(sym0 == NULL);

   // Alphabet with less than 2 symbols.
   sym0 = sym_new(alph01, NULL, 0);
   test_assert(sym0 == NULL);

   // Ambiguous symbols.
   sym0 = sym_new(alph02, NULL, 0);
   test_assert(sym0 == NULL);

   sym0 = sym_new(alph03, NULL, 0);
   test_assert(sym0 == NULL);

   // Wrong complements.
   sym0 = sym_new(alph04, NULL, 0);
   test_assert_critical(sym0 != NULL);
   test_assert(sym_count(sym0) == 4);
   sym_free(sym0);

   sym0 = sym_new(alph04, comp01, 0);
   test_assert(sym0 == NULL);

   sym0 = sym_new(alph04, comp02, 0);
   test_assert(sym0 == NULL);

   sym0 = sym_new(alph04, comp03, 0);
   test_assert(sym0 == NULL);

   sym0 = sym_new(alph04, comp04, 0);
   test_assert(sym0 == NULL);

   // Correct complement.
   sym0 = sym_new(alph04, comp05, 0);
   test_assert(sym0 != NULL);
   sym_free(sym0);

   sym0 = sym_new(alph05, comp05, 0);
   test_assert_critical(sym0 != NULL);
   test_assert(sym_count(sym0) == 8);
   sym_free(sym0);

   // Test output consistency.
   char * alph1[] = {"Aa","Cc","Gg","Tt","Nn",NULL};
   sym_t * sym1 = sym_new(alph1, NULL, 0);
   test_assert_critical(sym1 != NULL);
   test_assert(sym_count(sym1) == 5);
   test_assert(sym_character(0, sym1) == 'A');
   test_assert(sym_character(1, sym1) == 'C');
   test_assert(sym_character(2, sym1) == 'G');
   test_assert(sym_character(3, sym1) == 'T');
   test_assert(sym_character(4, sym1) == 'N');
   test_assert(sym_complement(0, sym1) == 0);
   test_assert(sym_complement(1, sym1) == 1);
   test_assert(sym_complement(2, sym1) == 2);
   test_assert(sym_complement(3, sym1) == 3);
   test_assert(sym_complement(4, sym1) == 4);
   test_assert(sym_index('A', sym1) == 0);
   test_assert(sym_index('a', sym1) == 0);
   test_assert(sym_index('C', sym1) == 1);
   test_assert(sym_index('c', sym1) == 1);
   test_assert(sym_index('G', sym1) == 2);
   test_assert(sym_index('g', sym1) == 2);
   test_assert(sym_index('T', sym1) == 3);
   test_assert(sym_index('t', sym1) == 3);
   test_assert(sym_index('N', sym1) == 4);
   test_assert(sym_index('n', sym1) == 4);
   test_assert(sym_index('x', sym1) == 0);
   test_assert(sym_index('K', sym1) == 0);
   test_assert(sym_index('$', sym1) == 0);
   sym_free(sym1);

   sym1 = sym_new(alph1, NULL, 4);
   test_assert_critical(sym1 != NULL);
   test_assert(sym_count(sym1) == 5);
   test_assert(sym_character(0, sym1) == 'A');
   test_assert(sym_character(1, sym1) == 'C');
   test_assert(sym_character(2, sym1) == 'G');
   test_assert(sym_character(3, sym1) == 'T');
   test_assert(sym_character(4, sym1) == 'N');
   test_assert(sym_complement(0, sym1) == 0);
   test_assert(sym_complement(1, sym1) == 1);
   test_assert(sym_complement(2, sym1) == 2);
   test_assert(sym_complement(3, sym1) == 3);
   test_assert(sym_complement(4, sym1) == 4);
   test_assert(sym_index('A', sym1) == 0);
   test_assert(sym_index('a', sym1) == 0);
   test_assert(sym_index('C', sym1) == 1);
   test_assert(sym_index('c', sym1) == 1);
   test_assert(sym_index('G', sym1) == 2);
   test_assert(sym_index('g', sym1) == 2);
   test_assert(sym_index('T', sym1) == 3);
   test_assert(sym_index('t', sym1) == 3);
   test_assert(sym_index('N', sym1) == 4);
   test_assert(sym_index('n', sym1) == 4);
   test_assert(sym_index('x', sym1) == 4);
   test_assert(sym_index('K', sym1) == 4);
   test_assert(sym_index('$', sym1) == 4);
   sym_free(sym1);

   char * comp12[] = {"AT", "TA", "CG", "GC", NULL};
   sym1 = sym_new(alph1, NULL, 4);
   test_assert_critical(sym1 != NULL);
   test_assert(sym_character(0, sym1) == 'A');
   test_assert(sym_character(1, sym1) == 'C');
   test_assert(sym_character(2, sym1) == 'G');
   test_assert(sym_character(3, sym1) == 'T');
   test_assert(sym_character(4, sym1) == 'N');
   test_assert(sym_complement(0, sym1) == 0);
   test_assert(sym_complement(1, sym1) == 1);
   test_assert(sym_complement(2, sym1) == 2);
   test_assert(sym_complement(3, sym1) == 3);
   test_assert(sym_complement(4, sym1) == 4);
   sym_set_complement(comp12, sym1);
   test_assert(sym_complement(0, sym1) == 3);
   test_assert(sym_complement(1, sym1) == 2);
   test_assert(sym_complement(2, sym1) == 1);
   test_assert(sym_complement(3, sym1) == 0);
   test_assert(sym_complement(4, sym1) == 4);
   sym_free(sym1);

   sym1 = sym_new(alph1, comp12, 4);
   test_assert_critical(sym1 != NULL);
   test_assert(sym_character(0, sym1) == 'A');
   test_assert(sym_character(1, sym1) == 'C');
   test_assert(sym_character(2, sym1) == 'G');
   test_assert(sym_character(3, sym1) == 'T');
   test_assert(sym_character(4, sym1) == 'N');
   test_assert(sym_complement(0, sym1) == 3);
   test_assert(sym_complement(1, sym1) == 2);
   test_assert(sym_complement(2, sym1) == 1);
   test_assert(sym_complement(3, sym1) == 0);
   test_assert(sym_complement(4, sym1) == 4);
   sym_free(sym1);

   char * alph2[] = {"Tt","Nn","Cc","Aa","Gg",NULL};
   sym_t * sym2 = sym_new(alph2, NULL, 2);
   test_assert_critical(sym2 != NULL);
   test_assert(sym_count(sym2) == 5);
   test_assert(sym_character(0, sym2) == 'T');
   test_assert(sym_character(1, sym2) == 'N');
   test_assert(sym_character(2, sym2) == 'C');
   test_assert(sym_character(3, sym2) == 'A');
   test_assert(sym_character(4, sym2) == 'G');
   test_assert(sym_index('A', sym2) == 3);
   test_assert(sym_index('a', sym2) == 3);
   test_assert(sym_index('C', sym2) == 2);
   test_assert(sym_index('c', sym2) == 2);
   test_assert(sym_index('G', sym2) == 4);
   test_assert(sym_index('g', sym2) == 4);
   test_assert(sym_index('T', sym2) == 0);
   test_assert(sym_index('t', sym2) == 0);
   test_assert(sym_index('N', sym2) == 1);
   test_assert(sym_index('n', sym2) == 1);
   test_assert(sym_index('x', sym2) == 2);
   test_assert(sym_index('K', sym2) == 2);
   test_assert(sym_index('$', sym2) == 2);
   test_assert(sym_complement(0, sym2) == 0);
   test_assert(sym_complement(1, sym2) == 1);
   test_assert(sym_complement(2, sym2) == 2);
   test_assert(sym_complement(3, sym2) == 3);
   test_assert(sym_complement(4, sym2) == 4);
   sym_free(sym2);

   char * comp22[] = {"AT", "TA", "CG", "GC", NULL};
   sym2 = sym_new(alph2, comp22, 1);
   test_assert_critical(sym2 != NULL);
   test_assert(sym_count(sym2) == 5);
   test_assert(sym_character(0, sym2) == 'T');
   test_assert(sym_character(1, sym2) == 'N');
   test_assert(sym_character(2, sym2) == 'C');
   test_assert(sym_character(3, sym2) == 'A');
   test_assert(sym_character(4, sym2) == 'G');
   test_assert(sym_index('A', sym2) == 3);
   test_assert(sym_index('a', sym2) == 3);
   test_assert(sym_index('C', sym2) == 2);
   test_assert(sym_index('c', sym2) == 2);
   test_assert(sym_index('G', sym2) == 4);
   test_assert(sym_index('g', sym2) == 4);
   test_assert(sym_index('T', sym2) == 0);
   test_assert(sym_index('t', sym2) == 0);
   test_assert(sym_index('N', sym2) == 1);
   test_assert(sym_index('n', sym2) == 1);
   test_assert(sym_index('x', sym2) == 1);
   test_assert(sym_index('K', sym2) == 1);
   test_assert(sym_index('$', sym2) == 1);
   test_assert(sym_complement(0, sym2) == 3);
   test_assert(sym_complement(1, sym2) == 1);
   test_assert(sym_complement(2, sym2) == 4);
   test_assert(sym_complement(3, sym2) == 0);
   test_assert(sym_complement(4, sym2) == 2);
   sym_free(sym2);

}




// Define test cases to be run (for export).
const test_case_t test_cases_from_file_1[] = {
   {"index_sym/sym_new",  test_sym_new},
   {NULL, NULL}, // Sentinel. //
};
