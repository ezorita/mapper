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
**    Alphabets must have at least 2 different symbols.
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

   // Check symbol consistency.
   uint8_t * used = calloc(SYM_TABLE_SIZE, 1);
   for (int i = 0; i < sym_count; i++) {
      int j = 0;
      while (alphabet[i][j] != '\0') {
         if (used[(uint8_t)alphabet[i][j]]) {
            free(used);
            return NULL;
         }
         used[(uint8_t)alphabet[i][j++]] = 1;
      }
   }
   free(used);

   // Alloc sym.
   sym_t * sym = malloc(sizeof(sym_t));
   if (sym == NULL)
      return NULL;
   
   // Initialize sym_t.
   *sym = (sym_t){sym_count, NULL, NULL, NULL};

   // (Re)alloc canonical symbol array.
   sym->sym_canon = malloc(sym_count+1);
   
   if (sym->sym_canon == NULL) {
      sym_free(sym);
      return NULL;
   }

   // Fill canonicals (The first character representing each symbol).
   for (int i = 0; i < sym_count; i++) {
      sym->sym_canon[i] = alphabet[i][0];
   }
   sym->sym_canon[sym->sym_count] = 0;

   // Fill sym table.
   sym->sym_table = malloc(SYM_TABLE_SIZE);

   if (sym->sym_table == NULL) {
      sym_free(sym);
      return NULL;
   }

   // Fill all characters with default symbol.
   memset(sym->sym_table, sym_default, SYM_TABLE_SIZE);
   
   // Set symbols.
   for (int i = 0; i < sym_count; i++) {
      int j = 0;
      while (alphabet[i][j] != '\0') sym->sym_table[(uint8_t)alphabet[i][j++]] = i;
   }

   // Alloc complements table.
   sym->com_table = malloc(sym_count + 1);

   if (sym->com_table == NULL) {
      sym_free(sym);
      return NULL;
   }
   
   // Set complements, identity if not set.
   if (sym_set_complement(complement, sym)) {
      sym_free(sym);
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
**    If complement is set to NULL, an identity relationship will be used, i.e. each 
**    symbol is its own complement.
**    Relationships of symbols not defined in complement are set to indentity.
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
   if (sym == NULL)
      return -1;

   if (sym->com_table == NULL)
      return -1;

   // Set identity base.
   for (int i = 0; i < sym->sym_count; i++) {
      sym->com_table[i] = i;
   }
   
   if (complement == NULL)
      return 0;

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

   if (sym->sym_canon == NULL || s >= sym->sym_count)
      return -1;

   // Return symbol representation.
   return sym->sym_canon[s];
}

uint8_t
sym_complement
(
  uint8_t   s,
  sym_t   * sym
)
{
   // Check arguments.
   if (sym == NULL)
      return -1;

   if (sym->com_table == NULL || s >= sym->sym_count)
      return -1;

   // Return symbol complement.
   return sym->com_table[s];
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


int
sym_count
(
  sym_t * sym
)
{
   // Check arguments.
   if (sym == NULL)
      return 0;
   
   // Return symbol count.
   return sym->sym_count;
}


// I/O functions
/* File structure:
** uint64_t x 1			magic_no  
** uint8_t  x 1			sym_count
** char     x sym_count+1	sym_canonicals
** uint8_t  x 256		sym_table
** uint8_t  x sym_count+1	com_table
*/

int
sym_file_write
(
  char   * filename,
  sym_t  * sym
)
{
   // Check arguments.
   if (filename == NULL || sym == NULL)
      return -1;
   
   // Open file.
   int fd = creat(filename, 0644);
   if (fd == -1)
      return -1;

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = SYM_FILE_MAGICNO;

   // Write magic.
   if (write(fd, &magic, sizeof(uint64_t)) == -1)
      goto close_and_error;

   // Write sym_count.
   if (write(fd, (uint8_t *)&(sym->sym_count), sizeof(uint8_t)) == -1)
      goto close_and_error;

   // Write sym_canon array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (char *)sym->sym_canon + e_cnt, (sym->sym_count + 1 - e_cnt)*sizeof(char));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(char);
   } while (e_cnt < sym->sym_count + 1);

   // Write sym_table array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint8_t *)sym->sym_table + e_cnt, (SYM_TABLE_SIZE - e_cnt)*sizeof(uint8_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < SYM_TABLE_SIZE);

   // Write com_table array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint8_t *)sym->com_table + e_cnt, (sym->sym_count + 1 - e_cnt)*sizeof(uint8_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < sym->sym_count + 1);


   close(fd);
   return 0;

 close_and_error:
   close(fd);
   return -1;
}

sym_t * 
sym_file_read
(
  char * filename
)
{
   // Check arguments.
   if (filename == NULL)
      return NULL;

   // Open file.
   int fd = open(filename, O_RDONLY);
   if (fd == -1)
      return NULL;

   // Alloc memory.
   sym_t * sym = malloc(sizeof(sym_t));
   if (sym == NULL)
      return NULL;
   // Set NULL pointers.
   sym->sym_canon = NULL;
   sym->sym_table = sym->com_table = NULL;
   
   // Read file.
   uint64_t magic;
   ssize_t b_cnt;
   ssize_t e_cnt;

   if (read(fd, &magic, sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;
   if (magic != SYM_FILE_MAGICNO)
      goto free_and_return;

   if (read(fd, &(sym->sym_count), sizeof(uint8_t)) < sizeof(uint8_t))
      goto free_and_return;
   if (sym->sym_count < 2)
      goto free_and_return;

   // Alloc canonical symbol array.
   sym->sym_canon = malloc((sym->sym_count + 1) * sizeof(uint8_t));
   if (sym->sym_canon == NULL)
      goto free_and_return;

   // Read symbol canonicals.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (char *)sym->sym_canon + e_cnt, (sym->sym_count + 1 - e_cnt) * sizeof(char));
      if (b_cnt == -1)
         goto free_and_return;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < sym->sym_count + 1);

   // Alloc symbol table.
   sym->sym_table = malloc(SYM_TABLE_SIZE * sizeof(uint8_t));
   if (sym->sym_table == NULL)
      goto free_and_return;

   // Read symbol table.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (uint8_t *)sym->sym_table + e_cnt, (SYM_TABLE_SIZE - e_cnt) * sizeof(uint8_t));
      if (b_cnt == -1)
         goto free_and_return;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < SYM_TABLE_SIZE);

   // Alloc symbol complement table.
   sym->com_table = malloc((sym->sym_count + 1) * sizeof(uint8_t));
   if (sym->com_table == NULL)
      goto free_and_return;

   // Read symbol table.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (uint8_t *)sym->com_table + e_cnt, (sym->sym_count + 1 - e_cnt) * sizeof(uint8_t));
      if (b_cnt == -1)
         goto free_and_return;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < sym->sym_count + 1);

   close(fd);
   return sym;

 free_and_return:
   close(fd);
   sym_free(sym);
   return NULL;
}

