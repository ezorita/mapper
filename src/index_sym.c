#include "index_sym.h"

// Interface data types.
struct sym_t {
   uint8_t   sym_count;
   char    * sym_canon;
   uint8_t * sym_table;
   uint8_t * com_table;
};

char * SYM_ALPHABET_DNA[] = {"Aa","Cc","Gg","Tt","Nn",NULL};
char * SYM_COMPLEMENT_DNA[] = {"AT","CG","GC","TA",NULL};
int    SYM_DEFAULT_DNA = 4;


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
   // Alloc names.
   sym_t   * sym  = NULL;
   uint8_t * used = NULL;

   // Check arguments.
   error_test_msg(alphabet == NULL, "argument 'alphabet' is NULL.");

   // Count number of symbols.
   int sym_count = 0;
   while (alphabet[sym_count] != NULL) sym_count++;
   
   // There must be at least two symbols to store info.
   error_test_msg(sym_count < 2, "less than 2 symbols defined.");
   error_test_msg(sym_count >= SYM_MAX_ALPHABET_SIZE, "alphabet size is greater than SYM_MAX_ALPHABET_SIZE.");

   // Check default symbol.
   error_test_msg(sym_default >= sym_count, "'sym_default' does not belong to alphabet.");

   // Check symbol consistency.
   used = calloc(SYM_TABLE_SIZE, 1);
   error_test_mem(used);

   for (int i = 0; i < sym_count; i++) {
      int j = 0;
      while (alphabet[i][j] != '\0') {
         if (used[(uint8_t)alphabet[i][j]]) {
            error_throw_msg("defined symbols must be unique.");
         }
         used[(uint8_t)alphabet[i][j++]] = 1;
      }
   }
   free(used);
   used = NULL;

   // Alloc sym.
   sym = malloc(sizeof(sym_t));
   error_test_mem(sym);
   
   // Initialize sym_t.
   *sym = (sym_t){sym_count, NULL, NULL, NULL};

   // (Re)alloc canonical symbol array.
   sym->sym_canon = malloc(sym_count+1);
   error_test_mem(sym->sym_canon);

   // Fill canonicals (The first character representing each symbol).
   for (int i = 0; i < sym_count; i++) {
      sym->sym_canon[i] = alphabet[i][0];
   }
   sym->sym_canon[sym->sym_count] = 0;

   // Fill sym table.
   sym->sym_table = malloc(SYM_TABLE_SIZE);
   error_test_mem(sym->sym_table);

   // Fill all characters with default symbol.
   memset(sym->sym_table, sym_default, SYM_TABLE_SIZE);
   
   // Set symbols.
   for (int i = 0; i < sym_count; i++) {
      int j = 0;
      while (alphabet[i][j] != '\0') sym->sym_table[(uint8_t)alphabet[i][j++]] = i;
   }

   // Alloc complements table.
   sym->com_table = calloc(sym_count + 1,sizeof(uint8_t));
   error_test_mem(sym->com_table);
   
   // Set complements, identity if not set.
   error_test(sym_set_complement(complement, sym) == -1);

   return sym;

 failure_return:
   free(used);
   sym_free(sym);
   return NULL;
}

sym_t *
sym_new_dna
(
  void
)
{
   return sym_new(SYM_ALPHABET_DNA,SYM_COMPLEMENT_DNA,SYM_DEFAULT_DNA);
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
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");
   error_test_msg(sym->com_table == NULL, "argument 'sym->com_table' is NULL.");

   // Set identity base.
   for (int i = 0; i < sym->sym_count + 1; i++) {
      sym->com_table[i] = i;
   }
   
   if (complement == NULL)
      return 0;

   // Check all relationships before applying them.
   int i = 0;
   char * rel;
   while ((rel = complement[i++]) != NULL) {
      error_test_msg(strlen(rel) < 2, "incorrect complement format.");
      error_test_msg(sym_is_canonical(rel[0], sym) == 0, "symbols must be canonicals.");
      error_test_msg(sym_is_canonical(rel[1], sym) == 0, "symbols must be canonicals.");
   }

   // Apply relationships.
   i = 0;
   while ((rel = complement[i++]) != NULL) {
      sym->com_table[sym->sym_table[(uint8_t)rel[0]]] = sym->sym_table[(uint8_t)rel[1]];
   }

   return 0;

 failure_return:
   return -1;
}

char
sym_character
(
  uint8_t   s,
  sym_t   * sym
)
{
   // Check arguments.
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");
   error_test_msg(sym->sym_canon == NULL, "argument 'sym->sym_canon' is NULL.");
   error_test_msg(s >= sym->sym_count, "symbol index out of bounds.");

   // Return symbol representation.
   return sym->sym_canon[s];
   
 failure_return:
   return -1;
}

int32_t
sym_complement
(
  uint8_t   s,
  sym_t   * sym
)
{
   // Check arguments.
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");
   error_test_msg(sym->com_table == NULL, "argument 'sym->com_table' is NULL.");
   error_test_msg(s >= sym->sym_count, "symbol index out of bounds.");

   // Return symbol complement.
   return sym->com_table[s];

 failure_return:
   return -1;
}

int32_t
sym_index
(
  char     c,
  sym_t  * sym
)
{
   // Check arguments.
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");

   // Return symbol index.
   return sym->sym_table[(uint8_t)c];

 failure_return:
   return -1;
}

int
sym_is_canonical
(
  char     c,
  sym_t  * sym
)
{
   // Check arguments.
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");

   // Check whether symbol is in canonicals.
   if (c == sym->sym_canon[sym->sym_table[(uint8_t)c]]) {
      return 1;
   } else {
      return 0;
   }

 failure_return:
   return -1;
}


int
sym_count
(
  sym_t * sym
)
{
   // Check arguments.
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");
   
   // Return symbol count.
   return sym->sym_count;

 failure_return:
   return -1;
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
   int fd = -1;
   // Check arguments.
   error_test_msg(filename == NULL, "argument 'filename' is NULL.");
   error_test_msg(sym == NULL, "argument 'sym' is NULL.");
   
   // Open file.
   fd = creat(filename, 0644);
   error_test_def(fd == -1);

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = SYM_FILE_MAGICNO;

   // Write magic.
   b_cnt = write(fd, &magic, sizeof(uint64_t));
   error_test_def(b_cnt == -1);

   // Write sym_count.
   b_cnt = write(fd, (uint8_t *)&(sym->sym_count), sizeof(uint8_t));
   error_test_def(b_cnt == -1);

   // Write sym_canon array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (char *)sym->sym_canon + e_cnt, (sym->sym_count + 1 - e_cnt)*sizeof(char));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(char);
   } while (e_cnt < sym->sym_count + 1);

   // Write sym_table array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint8_t *)sym->sym_table + e_cnt, (SYM_TABLE_SIZE - e_cnt)*sizeof(uint8_t));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < SYM_TABLE_SIZE);

   // Write com_table array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint8_t *)sym->com_table + e_cnt, (sym->sym_count + 1 - e_cnt)*sizeof(uint8_t));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < sym->sym_count + 1);

   close(fd);
   return 0;

 failure_return:
   if (fd != -1)
      close(fd);
   return -1;
}

sym_t * 
sym_file_read
(
  char * filename
)
{
   // Declare variables.
   int fd = -1;
   sym_t * sym = NULL;

   // Check arguments.
   error_test_msg (filename == NULL, "argument 'filename' is NULL.");

   // Open file.
   fd = open(filename, O_RDONLY);
   error_test_def(fd == -1);

   // Alloc memory.
   sym = malloc(sizeof(sym_t));
   error_test_mem(sym);

   // Set NULL pointers.
   sym->sym_canon = NULL;
   sym->sym_table = sym->com_table = NULL;
   
   // Read file.
   uint64_t magic;
   ssize_t b_cnt;
   ssize_t e_cnt;

   b_cnt = read(fd, &magic, sizeof(uint64_t));
   error_test_def(b_cnt == -1);
   error_test_msg(magic != SYM_FILE_MAGICNO, "unrecognized 'sym' file format (magicno).");

   b_cnt = read(fd, &(sym->sym_count), sizeof(uint8_t));
   error_test_def(b_cnt < -1);
   error_test_msg(sym->sym_count < 2, "incorrect 'sym' file format (sym_count < 2).");

   // Alloc canonical symbol array.
   sym->sym_canon = malloc((sym->sym_count + 1) * sizeof(uint8_t));
   error_test_mem(sym->sym_canon);

   // Read symbol canonicals.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (char *)sym->sym_canon + e_cnt, (sym->sym_count + 1 - e_cnt) * sizeof(char));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(char);
   } while (e_cnt < sym->sym_count + 1);

   // Alloc symbol table.
   sym->sym_table = malloc(SYM_TABLE_SIZE * sizeof(uint8_t));
   error_test_mem(sym->sym_table);

   // Read symbol table.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (uint8_t *)sym->sym_table + e_cnt, (SYM_TABLE_SIZE - e_cnt) * sizeof(uint8_t));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < SYM_TABLE_SIZE);

   // Alloc symbol complement table.
   sym->com_table = malloc((sym->sym_count + 1) * sizeof(uint8_t));
   error_test_mem(sym->com_table);

   // Read symbol table.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (uint8_t *)sym->com_table + e_cnt, (sym->sym_count + 1 - e_cnt) * sizeof(uint8_t));
      error_test_def(b_cnt == -1);
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < sym->sym_count + 1);

   close(fd);
   return sym;

 failure_return:
   if (fd != -1)
      close(fd);
   sym_free(sym);
   return NULL;
}

