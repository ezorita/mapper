#include "index_sym.h"

// Interface data types.
struct sym_t {
   uint8_t   sym_count;
   char    * sym_canon;
   uint8_t * sym_table;
   uint8_t * com_table;
};


// Interface function sources.
sym_t *
sym_new
(
  char    ** alphabet,
  char    ** complement,
  uint8_t    sym_default
)
/*
** Creates a new symbol alphabet.
**
**  alphabet:
**    must be an array of strings terminated with a NULL pointer. 
**    Each position of the array points to a string and represents a different symbol.
**    The string contains the valid characters that represent each symbol.
**    The canonical representation of each symbol is the first character of the string.
**
**    Example: 
**      // DNA symbols (A,C,T,G,N):
**      char * alphabet[6] = {"Aa","Cc","Gg","Tt","Nn",NULL};
**      // Canonicals are "ACGTN".
**
**  complement (optional, set to NULL):
**    An array of strings which defines the complement relationships between symbols,
**    terminated with a NULL pointer. If complement is set to NULL, an identity
**    relationship will be used, i.e. each symbol is its own complement.
**    Each string must contain two characters representing the complementary relationship.
**    Only canonical representations of symbols may be used.
**
**    Note: The relationship is unidirectional, so "AB" means that 'B' is the
**    complementary of 'A', but this has no implications on the complementary of 'B'.
**    For bidirectional complementarity, both complemens must be defined.
**
**      // DNA symbols (A,C,T,G,N):
**      char * alphabet[6] = {"Aa","Cc","Gg","Tt","Nn",NULL};
**      // Complementary relationships for DNA (A <-> T, C <-> G).
**      char * complement[5] = {"AT","TA","CG","GC",NULL};
**
**  sym_default:
**    The cardinality of the default symbol. Characters that do not represent any
**    symbol will be converted into the default symbol.
**
**    Example: 
**      // For DNA symbols (A,C,T,G,N) we want the default to be 'N'.
**      char * alphabet[6] = {"Aa","Cc","Gg","Tt","Nn",NULL};
**      // Complementary relationships for DNA (A <-> T, C <-> G).
**      char * complement[5] = {"AT","TA","CG","GC",NULL};
**      // 'N' is at position 4 of alphabet.
**      int sym_default = 4;
**      // All characters, except "ACGTacgt" will be considered 'N'.
**      sym_set_alphabet(alphabet, complement, sym_default, sym);
*/
{
   // Check arguments.
   if (alphabet == NULL)
      return NULL;

   // Count number of symbols.
   int sym_count = 0;
   while (alphabet[sym_count] != NULL) sym_count++;
   
   // There must be at least two symbols to store info.
   if (sym_count < 2 || sym_count >= SYM_MAX_ALPHABET_SIZE)
      return NULL;

   // Check default symbol.
   if (sym_default >= sym_count)
      return NULL;

   // Alloc sym.
   sym_t * sym = malloc(sizeof(sym_t));
   if (sym == NULL)
      return NULL;
   // Set symbol count.
   sym->sym_count = sym_count;

   // (Re)alloc canonical symbol array.
   if (sym->sym_canon == NULL) {
      sym->sym_canon = malloc(sym_count+1);
   } else {
      sym->sym_canon = realloc(sym->sym_canon, sym_count+1);
   }
   
   if (sym->sym_canon == NULL)
      return NULL;

   // Fill canonicals (The first character representing each symbol).
   for (int i = 0; i < sym_count; i++) {
      sym->sym_canon[i] = alphabet[i][0];
   }
   sym->sym_canon[sym->sym_count] = 0;

   // Fill sym table.
   if (sym->sym_table == NULL) {
      sym->sym_table = malloc(SYM_TABLE_SIZE);
   }

   // Fill all characters with default symbol.
   memset(sym->sym_table, sym_default, SYM_TABLE_SIZE);
   
   // Set symbols.
   for (int i = 0; i < sym_count; i++) {
      int j = 0;
      while (alphabet[i][j] != '\0') sym->sym_table[(uint8_t)alphabet[i][j++]] = i;
   }

   // (Re)alloc complement table.
   if (sym->com_table == NULL) {
      sym->com_table = malloc(sym_count + 1);
   } else {
      sym->com_table = realloc(sym->com_table, sym_count + 1);
   }

   if (sym->com_table == NULL)
      return NULL;
   
   // Set complement, identity if not set.
   if (complement == NULL) {
      for (int i = 0; i < sym_count; i++) {
         sym->com_table[i] = i;
      }
   } else {
      if (sym_set_complement(complement, sym))
         return NULL;
   }

   return sym;
}

void
sym_free
(
  sym_t  * sym
)
{
   if (sym != NULL) {
      free(sym->sym_canon);
      free(sym->sym_table);
      free(sym->com_table);
      free(sym);
   }

   return;
}


int
sym_set_complement
(
  char  ** complement,
  sym_t  * sym
)
/*
**
** Defines the complement table.
**
**  complement:
**    An array of strings which defines the complement relationships between symbols,
**    terminated with a NULL pointer.
**    Each string must contain two characters representing the complementary relationship.
**    Only canonical representations of symbols may be used.
**
**    Note: The relationship is unidirectional, so "AB" means that 'B' is the
**    complementary of 'A', but this has no implications on the complementary of 'B'.
**    For bidirectional complementarity, both complemens must be defined.
**
**    Example:
**      // Complementary relationships for DNA (A <-> T, C <-> G).
**      char * complement[5] = {"AT","TA","CG","GC",NULL};
**      sym_set_complement(complement, sym);
*/
{
   // Check arguments.
   if (sym == NULL || complement == NULL)
      return -1;

   if (sym->com_table == NULL)
      return -1;

   // Check all relationships before applying them.
   int i = 0;
   char * rel;
   while ((rel = complement[i++]) != NULL) {
      if (strlen(rel) < 2)
         return -1;
      if (!sym_is_canonical(rel[0], sym) || !sym_is_canonical(rel[1], sym))
         return -1;
   }

   // Apply relationships.
   i = 0;
   while ((rel = complement[i++]) != NULL) {
      sym->com_table[sym->sym_table[(uint8_t)rel[0]]] = sym->sym_table[(uint8_t)rel[1]];
   }

   return 0;
}

char
sym_character
(
  uint8_t   s,
  sym_t   * sym
)
{
   // Check arguments.
   if (sym == NULL)
      return -1;

   if (s >= sym->sym_count)
      return -1;

   // Return symbol representation.
   return sym->sym_canon[s];
}

uint8_t
sym_index
(
  char     c,
  sym_t  * sym
)
{
   // Check arguments.
   if (sym == NULL)
      return -1;

   // Return symbol index.
   return sym->sym_table[(uint8_t)c];
}

int
sym_is_canonical
(
  char     c,
  sym_t  * sym
)
{
   // Check arguments.
   if (sym == NULL)
      return -1;

   // Check whether symbol is in canonicals.
   if (c == sym->sym_canon[sym->sym_table[(uint8_t)c]])
         return 1;
   
   return 0;
}
