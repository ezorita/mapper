#include "index_txt.h"

// Interface data types.
struct txt_t {
   int64_t   txt_len;
   size_t    mem_size;
   int64_t   wil_cnt;
   sym_t   * sym;
   uint8_t * text;
};


// Interface functions source.
txt_t *
txt_new
(
  sym_t  * sym
)
{
   // Check arguments.
   if (sym == NULL)
      return NULL;
   
   // Alloc txt_t and text buffer.
   txt_t * txt = malloc(sizeof(txt_t));
   if (txt == NULL)
      return NULL;

   uint8_t * text = malloc(TXT_BUFFER_SIZE*sizeof(uint8_t));
   if (text == NULL) {
      free(txt);
      return NULL;
   }

   // Initialize and return.
   txt->txt_len = 0;
   txt->mem_size = TXT_BUFFER_SIZE;
   txt->wil_cnt = 0;
   txt->sym = sym;
   txt->text = text;

   return txt;
}


void
txt_free
(
  txt_t * txt
)
{
   if (txt != NULL) {
      free(txt->text);
      free(txt);
   }

   return;
}


int32_t
txt_sym
(
  int64_t pos,
  txt_t * txt
)
{
   // Check arguments.
   if (txt == NULL) 
      return -1;

   if (pos < 0 || pos >= txt->txt_len)
      return -1;

   return txt->text[pos];
}


uint8_t *
txt_sym_range
(
  int64_t     beg,
  int64_t     len,
  txt_t     * txt
)
{
   // Check arguments.
   if (txt == NULL)
      return NULL;
   
   if (beg < 0 || beg + len > txt->txt_len)
      return NULL;

   // Since we are not compressing, return pointer to text.
   return txt->text + beg;
}


int
txt_append
(
 char   * text,
 txt_t  * txt
)
{
   if (text == NULL || txt == NULL)
      return -1;

   // Realloc text structure if necessary.
   size_t tlen = strlen(text);
   if (txt->txt_len + tlen >= txt->mem_size) {
      size_t new_mem_size = txt->mem_size * 2;
      txt->text = realloc(txt->text, new_mem_size * sizeof(uint8_t));
      if (txt->text == NULL)
         return -1;
      txt->mem_size = new_mem_size;
   }

   // Append text symbols.
   for (int i = 0; i < tlen; i++) {
      txt->text[txt->txt_len + i] = sym_index(text[i], txt->sym);
   }

   txt->txt_len += tlen;

   return 0;
}


int
txt_append_wildcard
(
  txt_t  * txt
)
/*
** Wildcards are represented with the (n_symbols+1)-th symbol.
*/
{
   // Check arguments.
   if (txt == NULL)
      return -1;

   // Realloc text structure if necessary.
   if (txt->txt_len >= txt->mem_size) {
      size_t new_mem_size = txt->mem_size + 1;
      txt->text = realloc(txt->text, new_mem_size * sizeof(uint8_t));
      if (txt->text == NULL)
         return -1;
      txt->mem_size = new_mem_size;
   }

   // Append wildcard.
   txt->text[txt->txt_len++] = sym_count(txt->sym);
   txt->wil_cnt++;
   
   return 0;
}


// Helper functions.
int64_t
txt_length
(
  txt_t * txt
)
{
   // Check arguments.
   if (txt == NULL)
      return -1;

   return txt->txt_len;
}


sym_t *
txt_get_symbols
(
 txt_t * txt
)
{
   // Check arguments.
   if (txt == NULL)
      return NULL;
   
   return txt->sym;
}


int64_t
txt_wildcard_count
(
  txt_t * txt
)
{
   // Check arguments.
   if (txt == NULL) 
      return -1;

   return txt->wil_cnt;
}

// I/O functions.
int
txt_file_write
(
  char   * filename,
  txt_t  * txt
)
{
   // Check arguments.
   if (filename == NULL || txt == NULL)
      return -1;
   
   // Open file.
   int fd = creat(filename, 0644);
   if (fd == -1)
      return -1;

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = TXT_FILE_MAGICNO;

   // Write magic.
   if (write(fd, &magic, sizeof(uint64_t)) == -1)
      goto close_and_error;

   // Write txt_len.
   if (write(fd, (int64_t *)&(txt->txt_len), sizeof(int64_t)) == -1)
      goto close_and_error;

   // write wil_cnt.
   if (write(fd, (int64_t *)&(txt->wil_cnt), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write text array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (uint8_t *)txt->text + e_cnt, (txt->txt_len - e_cnt)*sizeof(uint8_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < txt->txt_len);

   close(fd);
   return 0;

 close_and_error:
   close(fd);
   return -1;
}


txt_t *
txt_file_read
(
  char   * filename,
  sym_t  * sym
)
{
      // Check arguments.
   if (filename == NULL || sym == NULL)
      return NULL;

   // Open file.
   int fd = open(filename, O_RDONLY);
   if (fd == -1)
      return NULL;

   // Alloc memory.
   txt_t * txt = malloc(sizeof(txt_t));
   if (txt == NULL)
      return NULL;
   // Set NULL pointers.
   txt->text = NULL;
   
   // Read file.
   uint64_t magic;
   ssize_t b_cnt;
   ssize_t e_cnt;

   // Read magic number.
   if (read(fd, &magic, sizeof(uint64_t)) < sizeof(uint64_t))
      goto free_and_return;
   if (magic != TXT_FILE_MAGICNO)
      goto free_and_return;

   // Read 'txt_len'.
   if (read(fd, &(txt->txt_len), sizeof(int64_t)) < sizeof(int64_t))
      goto free_and_return;
   if (txt->txt_len < 1)
      goto free_and_return;

   // Read 'wil_cnt'.
   if (read(fd, &(txt->wil_cnt), sizeof(int64_t)) < sizeof(int64_t))
      goto free_and_return;
   if (txt->txt_len < 1)
      goto free_and_return;

   // Alloc text.
   txt->text = malloc(txt->txt_len * sizeof(uint8_t));
   if (txt->text == NULL)
      goto free_and_return;

   // Set 'mem_size'.
   txt->mem_size = txt->txt_len;

   // Read 'text'.
   e_cnt = 0;
   do {
      b_cnt = read(fd, (char *)txt->text + e_cnt, (txt->txt_len - e_cnt) * sizeof(uint8_t));
      if (b_cnt == -1)
         goto free_and_return;
      e_cnt += b_cnt / sizeof(uint8_t);
   } while (e_cnt < txt->txt_len);

   // Set symbol alphabet.
   txt->sym = sym;

   return txt;
   
 free_and_return:
   close(fd);
   txt_free(txt);
   return NULL;
}
