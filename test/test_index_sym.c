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

void
test_sym_set_complement
(void)
{
   char * alph0[] = {"Aa","Bb","Cc","Dd",NULL};
   char * comp00[] = {"AB","CD",NULL};
   char * comp01[] = {"AC","BD",NULL};
   char * comp02[] = {"AB","BA","CD","DC",NULL};
   char * comp03[] = {NULL};
   char * comp04[] = {"AA","BA","CD","DD",NULL};

   test_assert(sym_set_complement(NULL, NULL) == -1);
   test_assert(sym_set_complement(comp00, NULL) == -1);

   sym_t * sym;
   sym = sym_new(alph0, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_complement(0,sym) == 0);
   test_assert(sym_complement(1,sym) == 1);
   test_assert(sym_complement(2,sym) == 2);
   test_assert(sym_complement(3,sym) == 3);

   test_assert(sym_set_complement(comp00, sym) == 0);
   test_assert(sym_complement(0,sym) == 1);
   test_assert(sym_complement(1,sym) == 1);
   test_assert(sym_complement(2,sym) == 3);
   test_assert(sym_complement(3,sym) == 3);

   test_assert(sym_set_complement(comp01, sym) == 0);
   test_assert(sym_complement(0,sym) == 2);
   test_assert(sym_complement(1,sym) == 3);
   test_assert(sym_complement(2,sym) == 2);
   test_assert(sym_complement(3,sym) == 3);

   test_assert(sym_set_complement(comp02, sym) == 0);
   test_assert(sym_complement(0,sym) == 1);
   test_assert(sym_complement(1,sym) == 0);
   test_assert(sym_complement(2,sym) == 3);
   test_assert(sym_complement(3,sym) == 2);

   test_assert(sym_set_complement(comp03, sym) == 0);
   test_assert(sym_complement(0,sym) == 0);
   test_assert(sym_complement(1,sym) == 1);
   test_assert(sym_complement(2,sym) == 2);
   test_assert(sym_complement(3,sym) == 3);

   test_assert(sym_set_complement(comp04, sym) == 0);
   test_assert(sym_complement(0,sym) == 0);
   test_assert(sym_complement(1,sym) == 0);
   test_assert(sym_complement(2,sym) == 3);
   test_assert(sym_complement(3,sym) == 3);

   test_assert(sym_set_complement(NULL, sym) == 0);
   test_assert(sym_complement(0,sym) == 0);
   test_assert(sym_complement(1,sym) == 1);
   test_assert(sym_complement(2,sym) == 2);
   test_assert(sym_complement(3,sym) == 3);

   char * wrong00[] = {"AB","BA","CD","D",NULL};
   char * wrong01[] = {"Ab","BA","CD","DC",NULL};
   char * wrong02[] = {"AB","BA","CD","DF",NULL};
   char * wrong03[] = {"AB","BA","CD","Dc",NULL};
   
   test_assert(sym_set_complement(wrong00, sym) == -1);
   test_assert(sym_set_complement(wrong01, sym) == -1);
   test_assert(sym_set_complement(wrong02, sym) == -1);
   test_assert(sym_set_complement(wrong03, sym) == -1);

   sym_free(sym);
}

void
test_sym_character
(void)
{
   sym_t * sym;
   char * alph0[] = {"Aa","Bb","Cc","Dd","Ee",NULL};
   char * alph1[] = {"aA","Bb","cC","Dd","eE",NULL};

   test_assert(sym_character(0, NULL) == -1);
   
   sym = sym_new(alph0, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_character(0, sym) == 'A');
   test_assert(sym_character(1, sym) == 'B');
   test_assert(sym_character(2, sym) == 'C');
   test_assert(sym_character(3, sym) == 'D');
   test_assert(sym_character(4, sym) == 'E');

   test_assert(sym_character(5, sym) == -1);
   test_assert(sym_character(6, sym) == -1);
   test_assert(sym_character(-1, sym) == -1);
   sym_free(sym);

   sym = sym_new(alph1, NULL, 2);
   test_assert_critical(sym != NULL);
   test_assert(sym_character(0, sym) == 'a');
   test_assert(sym_character(1, sym) == 'B');
   test_assert(sym_character(2, sym) == 'c');
   test_assert(sym_character(3, sym) == 'D');
   test_assert(sym_character(4, sym) == 'e');

   test_assert(sym_character(5, sym) == -1);
   test_assert(sym_character(6, sym) == -1);
   test_assert(sym_character(-1, sym) == -1);
   sym_free(sym);   
}

void
test_sym_index
(void)
{
   sym_t * sym;
   char * alph0[] = {"Aa","Bb","Cc","Dd","Ee",NULL};
   char * alph1[] = {"aA","Bb","cC","Dd","eE",NULL};

   test_assert(sym_index(0, NULL) == (uint8_t)-1);
   
   sym = sym_new(alph0, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_character(sym_index('A', sym), sym) == 'A');
   test_assert(sym_character(sym_index('a', sym), sym) == 'A');
   test_assert(sym_character(sym_index('B', sym), sym) == 'B');
   test_assert(sym_character(sym_index('b', sym), sym) == 'B');
   test_assert(sym_character(sym_index('C', sym), sym) == 'C');
   test_assert(sym_character(sym_index('c', sym), sym) == 'C');
   test_assert(sym_character(sym_index('D', sym), sym) == 'D');
   test_assert(sym_character(sym_index('d', sym), sym) == 'D');
   test_assert(sym_character(sym_index('E', sym), sym) == 'E');
   test_assert(sym_character(sym_index('e', sym), sym) == 'E');
   test_assert(sym_character(sym_index('F', sym), sym) == 'A');
   test_assert(sym_character(sym_index('f', sym), sym) == 'A');
   test_assert(sym_character(sym_index('Q', sym), sym) == 'A');
   test_assert(sym_character(sym_index('$', sym), sym) == 'A');
   test_assert(sym_character(sym_index('9', sym), sym) == 'A');
   sym_free(sym);

   sym = sym_new(alph1, NULL, 2);
   test_assert_critical(sym != NULL);
   test_assert(sym_character(sym_index('A', sym), sym) == 'a');
   test_assert(sym_character(sym_index('a', sym), sym) == 'a');
   test_assert(sym_character(sym_index('B', sym), sym) == 'B');
   test_assert(sym_character(sym_index('b', sym), sym) == 'B');
   test_assert(sym_character(sym_index('C', sym), sym) == 'c');
   test_assert(sym_character(sym_index('c', sym), sym) == 'c');
   test_assert(sym_character(sym_index('D', sym), sym) == 'D');
   test_assert(sym_character(sym_index('d', sym), sym) == 'D');
   test_assert(sym_character(sym_index('E', sym), sym) == 'e');
   test_assert(sym_character(sym_index('e', sym), sym) == 'e');
   test_assert(sym_character(sym_index('F', sym), sym) == 'c');
   test_assert(sym_character(sym_index('f', sym), sym) == 'c');
   test_assert(sym_character(sym_index('Q', sym), sym) == 'c');
   test_assert(sym_character(sym_index('$', sym), sym) == 'c');
   test_assert(sym_character(sym_index('9', sym), sym) == 'c');
   sym_free(sym);

}

void
test_sym_is_canonical
(void)
{
   sym_t * sym;
   char * alph0[] = {"Aa","Bb","Cc","Dd","Ee",NULL};
   char * alph1[] = {"aA","Bb","cC","Dd","eE",NULL};

   test_assert(sym_is_canonical(0, NULL) == -1);
   
   sym = sym_new(alph0, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_is_canonical('A',sym) == 1);
   test_assert(sym_is_canonical('B',sym) == 1);
   test_assert(sym_is_canonical('C',sym) == 1);
   test_assert(sym_is_canonical('D',sym) == 1);
   test_assert(sym_is_canonical('E',sym) == 1);
   test_assert(sym_is_canonical('a',sym) == 0);
   test_assert(sym_is_canonical('b',sym) == 0);
   test_assert(sym_is_canonical('c',sym) == 0);
   test_assert(sym_is_canonical('d',sym) == 0);
   test_assert(sym_is_canonical('e',sym) == 0);
   test_assert(sym_is_canonical('#',sym) == 0);
   test_assert(sym_is_canonical('!',sym) == 0);
   test_assert(sym_is_canonical('8',sym) == 0);
   test_assert(sym_is_canonical('^',sym) == 0);
   test_assert(sym_is_canonical('%',sym) == 0);
   test_assert(sym_is_canonical('$',sym) == 0);
   sym_free(sym);

   sym = sym_new(alph1, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_is_canonical('A',sym) == 0);
   test_assert(sym_is_canonical('B',sym) == 1);
   test_assert(sym_is_canonical('C',sym) == 0);
   test_assert(sym_is_canonical('D',sym) == 1);
   test_assert(sym_is_canonical('E',sym) == 0);
   test_assert(sym_is_canonical('a',sym) == 1);
   test_assert(sym_is_canonical('b',sym) == 0);
   test_assert(sym_is_canonical('c',sym) == 1);
   test_assert(sym_is_canonical('d',sym) == 0);
   test_assert(sym_is_canonical('e',sym) == 1);
   test_assert(sym_is_canonical('#',sym) == 0);
   test_assert(sym_is_canonical('!',sym) == 0);
   test_assert(sym_is_canonical('8',sym) == 0);
   test_assert(sym_is_canonical('^',sym) == 0);
   test_assert(sym_is_canonical('%',sym) == 0);
   test_assert(sym_is_canonical('$',sym) == 0);
   sym_free(sym);   
}

void
test_sym_complement
(void)
{
   char * alph0[] = {"Aa","Bb","Cc","Dd",NULL};
   char * comp00[] = {"AB","CD",NULL};
   char * comp01[] = {"AC","BD",NULL};
   char * comp02[] = {"AB","BA","CD","DC",NULL};
   char * comp03[] = {NULL};

   sym_t * sym;
   sym = sym_new(alph0, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_complement(0, sym) == 0);
   test_assert(sym_complement(1, sym) == 1);
   test_assert(sym_complement(2, sym) == 2);
   test_assert(sym_complement(3, sym) == 3);

   test_assert(sym_complement(4, sym) == (uint8_t)-1);
   test_assert(sym_complement(-1, sym) == (uint8_t)-1);

   test_assert(sym_set_complement(comp00, sym) == 0);
   test_assert(sym_complement(0, sym) == 1);
   test_assert(sym_complement(1, sym) == 1);
   test_assert(sym_complement(2, sym) == 3);
   test_assert(sym_complement(3, sym) == 3);

   test_assert(sym_complement(4, sym) == (uint8_t)-1);
   test_assert(sym_complement(-1, sym) == (uint8_t)-1);

   test_assert(sym_set_complement(comp01, sym) == 0);
   test_assert(sym_complement(0, sym) == 2);
   test_assert(sym_complement(1, sym) == 3);
   test_assert(sym_complement(2, sym) == 2);
   test_assert(sym_complement(3, sym) == 3);

   test_assert(sym_complement(4, sym) == (uint8_t)-1);
   test_assert(sym_complement(-1, sym) == (uint8_t)-1);

   test_assert(sym_set_complement(comp02, sym) == 0);
   test_assert(sym_complement(0, sym) == 1);
   test_assert(sym_complement(1, sym) == 0);
   test_assert(sym_complement(2, sym) == 3);
   test_assert(sym_complement(3, sym) == 2);

   test_assert(sym_complement(4, sym) == (uint8_t)-1);
   test_assert(sym_complement(-1, sym) == (uint8_t)-1);

   test_assert(sym_set_complement(comp03, sym) == 0);
   test_assert(sym_complement(0, sym) == 0);
   test_assert(sym_complement(1, sym) == 1);
   test_assert(sym_complement(2, sym) == 2);
   test_assert(sym_complement(3, sym) == 3);

   test_assert(sym_complement(4, sym) == (uint8_t)-1);
   test_assert(sym_complement(-1, sym) == (uint8_t)-1);

   sym_free(sym);
}

void
test_sym_count
(void)
{
   char * alph0[] = {"a","b","c","d",NULL};
   char * alph1[] = {"a","b","c","d","e",NULL};
   char * alph2[] = {"a","b","c","d","e","f",NULL};
   char * alph3[] = {"a","b","c","d","e","fF","gG","h","iI","Jj",NULL};
   char * alph4[] = {"a","b",NULL};

   sym_t * sym;
   sym = sym_new(alph0, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_count(sym) == 4);
   sym_free(sym);

   sym = sym_new(alph1, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_count(sym) == 5);
   sym_free(sym);

   sym = sym_new(alph2, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_count(sym) == 6);
   sym_free(sym);

   sym = sym_new(alph3, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_count(sym) == 10);
   sym_free(sym);

   sym = sym_new(alph4, NULL, 0);
   test_assert_critical(sym != NULL);
   test_assert(sym_count(sym) == 2);
   sym_free(sym);
}

// Define test cases to be run (for export).
const test_case_t test_cases_from_file_1[] = {
   {"index_sym/sym_new",             test_sym_new},
   {"index_sym/sym_set_complement",  test_sym_set_complement},
   {"index_sym/sym_character",       test_sym_character},
   {"index_sym/sym_complement",      test_sym_complement},
   {"index_sym/sym_index",           test_sym_index},
   {"index_sym/sym_is_canonical",    test_sym_is_canonical},
   {"index_sym/sym_count",           test_sym_count},
   {NULL, NULL}, // Sentinel. //
};
