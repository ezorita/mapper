#include "index_txt.h"

// Interface data types.
struct txt_t {
   size_t      mmap_len;
   void      * mmap_ptr;
   size_t      mem_txt;
   size_t      mem_seq;
   int64_t     txt_len;
   int64_t     seq_cnt;
   int64_t     wil_cnt;
   int64_t     rc_flag;
   int64_t   * seq_len;
   int64_t   * seq_beg;
   char     ** seq_name;
   sym_t     * sym;
   uint8_t   * text;
};

// Private function headers.
int64_t   seq_bisect  (int64_t beg, int64_t end, int64_t target, int64_t * vals);


// Interface functions source.
txt_t *
txt_new
(
  sym_t  * sym
)
{
   // Variable names carrying alloc pointers.
   int64_t  * seq_len = NULL, * seq_beg = NULL;
   char    ** seq_name = NULL;
   uint8_t  * text = NULL;
   txt_t    * txt = NULL;

   // Check arguments.
   if (sym == NULL)
      goto failure_return;
   
   // Alloc txt_t.
   txt = malloc(sizeof(txt_t));

   if (txt == NULL)
      goto failure_return;

   // Alloc sequence memory.
   seq_len  = malloc(TXT_SEQ_SIZE*sizeof(uint64_t));
   seq_beg  = malloc(TXT_SEQ_SIZE*sizeof(uint64_t));
   seq_name = malloc(TXT_SEQ_SIZE*sizeof(char *));

   if (!seq_len || !seq_beg || !seq_name)
      goto failure_return;

   // Alloc text memory.
   text = malloc(TXT_BUFFER_SIZE*sizeof(uint8_t));

   if (text == NULL)
      goto failure_return;

   // Initialize and return.
   txt->mmap_len = 0;
   txt->mmap_ptr = NULL;
   txt->mem_txt  = TXT_BUFFER_SIZE;
   txt->mem_seq  = TXT_SEQ_SIZE;
   txt->txt_len  = 0;
   txt->seq_cnt  = 0;
   txt->wil_cnt  = 0;
   txt->rc_flag  = 0;
   txt->seq_len  = seq_len;
   txt->seq_beg  = seq_beg;
   txt->seq_name = seq_name;
   txt->sym      = sym;
   txt->text     = text;

   return txt;

 failure_return:
   free(seq_len);
   free(seq_beg);
   free(seq_name);
   free(text);
   free(txt);
   return NULL;
}


void
txt_free
(
  txt_t * txt
)
{
   if (txt != NULL) {
      if (txt->mmap_ptr == NULL) {
         for (int64_t i = 0; i < txt->seq_cnt; i++) {
            free(txt->seq_name[i]);
         }
         free(txt->seq_len);
         free(txt->seq_beg);
         free(txt->seq_name);
         free(txt->text);
      } else {
         free(txt->seq_name);
         munmap(txt->mmap_ptr, txt->mmap_len);
      }
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
   if (txt->txt_len + tlen >= txt->mem_txt) {
      size_t new_mem_size = txt->mem_txt * 2;
      txt->text = realloc(txt->text, new_mem_size * sizeof(uint8_t));
      if (txt->text == NULL)
         return -1;
      txt->mem_txt = new_mem_size;
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
   if (txt->txt_len >= txt->mem_txt) {
      size_t new_mem_size = txt->mem_txt + 1;
      txt->text = realloc(txt->text, new_mem_size * sizeof(uint8_t));
      if (txt->text == NULL)
         return -1;
      txt->mem_txt = new_mem_size;
   }

   // Append wildcard.
   txt->text[txt->txt_len++] = sym_count(txt->sym);
   txt->wil_cnt++;
   
   return 0;
}


int
txt_commit_seq
(
  char   * seqname,
  txt_t  * txt
)
{
   // Check arguments.
   if (txt == NULL || seqname == NULL)
      return -1;
   
   // Check seqname length.
   if (strlen(seqname) < 1)
      return -1;

   // Check seqname uniqueness.
   for (int64_t i = 0; i < txt->seq_cnt; i++) {
      if (strcmp(seqname, txt->seq_name[i]) == 0)
         return -1;
   }

   // Realloc structure if necessary.
   if (txt->txt_len >= txt->mem_txt) {
      size_t new_mem_size = txt->mem_txt + 1;
      txt->text = realloc(txt->text, new_mem_size * sizeof(uint8_t));
      if (txt->text == NULL)
         return -1;
      txt->mem_txt = new_mem_size;
   }

   if (txt->seq_cnt >= txt->mem_seq) {
      size_t new_mem_size = 2*txt->mem_seq;
      txt->seq_beg  = realloc(txt->seq_beg , new_mem_size * sizeof(uint64_t));
      txt->seq_len  = realloc(txt->seq_len , new_mem_size * sizeof(uint64_t));
      txt->seq_name = realloc(txt->seq_name, new_mem_size * sizeof(char *));
      if (!txt->seq_beg || !txt->seq_len || !txt->seq_name)
         return -1;
      txt->mem_seq = new_mem_size;
   }

   // Get current sequence beginning and length.
   uint64_t beg = 0;
   if (txt->seq_cnt > 0) {
      beg = txt->seq_beg[txt->seq_cnt-1] + txt->seq_len[txt->seq_cnt-1];
   }

   // Append wildcard.
   txt_append_wildcard(txt);

   // Store sequence info.
   txt->seq_beg[txt->seq_cnt]  = beg;
   txt->seq_len[txt->seq_cnt]  = txt->txt_len - beg;
   txt->seq_name[txt->seq_cnt] = strdup(seqname);

   if (txt->seq_name[txt->seq_cnt] == NULL)
      return -1;

   // Update sequence count.
   txt->seq_cnt++;
   
   return 0;
}


int
txt_commit_rc
(
  txt_t  * txt
)
{
   // Check arguments.
   if (txt == NULL)
      return -1;

   if (txt->txt_len < 1)
      return 0;

   // Realloc structure if necessary.
   if (2*txt->txt_len >= txt->mem_txt) {
      size_t new_mem_size = 2*(txt->mem_txt + 1);
      txt->text = realloc(txt->text, new_mem_size * sizeof(uint8_t));
      if (txt->text == NULL)
         return -1;
      txt->mem_txt = new_mem_size;
   }
   
   // Append reverse complement of text, including wildcards.
   sym_t  * sym = txt_get_symbols(txt);
   int64_t tlen = txt_length(txt);
   int     nsym = sym_count(sym);

   // Add last wildcard if not present.
   if (txt->text[tlen-1] < nsym) {
      txt_append_wildcard(txt);
   } else {
      tlen--;
   }
   
   // Append reverse complement of text.
   for (int i = 1; i <= tlen; i++) {
      if (__builtin_expect(txt->text[tlen-i] < nsym, 1)) {
         txt->text[tlen+i] = sym_complement(txt_sym(tlen-i,txt), sym);
         txt->txt_len++;
      } else {
         txt_append_wildcard(txt);
      }
   }

   // Append last wildcard.
   if (txt->text[txt_length(txt)-1] < nsym) {
      txt_append_wildcard(txt);
   }

   // Set RevComp flag.
   txt->rc_flag = 1;

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


int64_t
txt_seq_count
(
  txt_t * txt
)
{
   // Check arguments.
   if (txt == NULL) 
      return -1;

   return txt->seq_cnt;
}


int64_t
txt_seq_start
(
  int32_t    seq,
  txt_t    * txt
)
{
   if (txt == NULL)
      return -1;
   if (seq < 0 || seq >= txt->seq_cnt)
      return -1;

   return txt->seq_beg[seq];
}


int64_t
txt_seq_length
(
  int32_t    seq,
  txt_t    * txt
)
{
   if (txt == NULL)
      return -1;
   if (seq < 0 || seq >= txt->seq_cnt)
      return -1;

   return txt->seq_len[seq];
}


char *
txt_seq_name
(
  int32_t    seq,
  txt_t    * txt
)
{
   if (txt == NULL)
      return NULL;
   if (seq < 0 || seq >= txt->seq_cnt)
      return NULL;

   return txt->seq_name[seq];
}


char *
txt_pos_to_str
(
  int64_t    pos,
  txt_t    * txt
)
{
   if (txt == NULL)
      return NULL;
   if (pos < 0 || pos >= txt->txt_len)
      return NULL;

   // Check whether text has RC.
   int strand = 0;
   if (txt->rc_flag) {
      strand = (pos >= txt->txt_len/2 ? 1 : 0);
      pos    = (pos >= txt->txt_len/2 ? txt->txt_len - 2 - pos : pos);
   }

   // Find sequence name.
   int seq_id = seq_bisect(0, txt->seq_cnt, pos, txt->seq_beg);

   // Generate string.
   char * pos_str = malloc(strlen(txt->seq_name[seq_id])+30);
   if (pos_str == NULL)
      return NULL;

   sprintf(pos_str, "%s:%ld:%c", txt->seq_name[seq_id], pos - txt->seq_beg[seq_id] + 1, (strand ? '-' : '+'));

   return pos_str;
}


int64_t
txt_str_to_pos
(
  char   * str,
  txt_t  * txt
)
{
   if (str == NULL || txt == NULL)
      return -1;

   char * pos_str = strdup(str);
   if (pos_str == NULL)
      return -1;

   // Tokenize string.
   char * seq_name = strtok(pos_str, ":");
   char * seq_pos  = strtok(NULL, ":");
   char * seq_fw   = strtok(NULL, ":");

   // Incorrect format.
   if (seq_pos == NULL)
      return -1;

   // Read strand (default is '+').
   int strand = 1;
   if (seq_fw == NULL || seq_fw[0] == '+')
      strand = 0;

   // Find seq index by name.
   int seq_id = -1;
   for (int64_t i = 0; i < txt->seq_cnt; i++) {
      if (strcmp(seq_name, txt->seq_name[i]) == 0) {
         seq_id = i;
         break;
      }
   }
   if (seq_id < 0)
      return -1;

   // Convert position to int.
   int64_t pos = atol(seq_pos);
   if (pos < 1 || pos > txt->seq_len[seq_id])
      return -1;

   // Find position in sequence.
   pos = txt->seq_beg[seq_id] + pos - 1;
   if (strand) {
      pos = txt->txt_len - 2 - pos;
   }         

   free(pos_str);
   
   // Return absolute position.
   return pos;
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

   // write seq_cnt.
   if (write(fd, (int64_t *)&(txt->seq_cnt), sizeof(int64_t)) == -1)
      goto close_and_error;

   // write wil_cnt.
   if (write(fd, (int64_t *)&(txt->wil_cnt), sizeof(int64_t)) == -1)
      goto close_and_error;

   // write rc_flag.
   if (write(fd, (int64_t *)&(txt->rc_flag), sizeof(int64_t)) == -1)
      goto close_and_error;

   // write seq_len array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (int64_t *)txt->seq_len + e_cnt, (txt->seq_cnt - e_cnt)*sizeof(int64_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(int64_t);
   } while (e_cnt < txt->seq_cnt);


   // write seq_beg array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (int64_t *)txt->seq_beg + e_cnt, (txt->seq_cnt - e_cnt)*sizeof(int64_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(int64_t);
   } while (e_cnt < txt->seq_cnt);

   // write seq_names.
   for (int64_t i = 0; i < txt->seq_cnt; i++) {
      size_t namelen = strlen(txt->seq_name[i])+1;
      e_cnt = 0;
      do {
         b_cnt  = write(fd, (char *)txt->seq_name[i] + e_cnt, (namelen - e_cnt)*sizeof(char));
         if (b_cnt == -1)
            goto close_and_error;
         e_cnt += b_cnt / sizeof(char);
      } while (e_cnt < namelen);
   }

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

   txt_t * txt = malloc(sizeof(txt_t));
   if (txt == NULL)
      goto free_and_return;

   // Set NULL pointers.
   txt->mmap_len = 0;
   txt->mmap_ptr = NULL;
   txt->seq_len  = NULL;
   txt->seq_beg  = NULL;
   txt->seq_name = NULL;
   txt->text     = NULL;

   // Get file len and mmap file.
   struct stat sb;
   fstat(fd, &sb);
   if (sb.st_size <= 16)
      goto free_and_return;

   int64_t * data = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
   if (data == NULL)
      goto free_and_return;

   txt->mmap_len = sb.st_size;
   txt->mmap_ptr = (void *) data;

   // Read file.
   // Read magic number.
   uint64_t magic = data[0];
   if (magic != TXT_FILE_MAGICNO)
      goto free_and_return;

   // Read 'txt_len'.
   txt->txt_len = data[1];
   if (txt->txt_len < 1)
      goto free_and_return;

   // Read 'seq_cnt'.
   txt->seq_cnt = data[2];
   if (txt->seq_cnt < 1)
      goto free_and_return;

   // Read 'seq_cnt'.
   txt->wil_cnt = data[3];

   // Read 'rc_flag'.
   txt->rc_flag = data[4];
   if (txt->rc_flag < 0 || txt->rc_flag > 1)
      goto free_and_return;
   
   // Pointer to seq_beg array.
   txt->seq_len = data + 5;

   // Pointer to seq_len array.
   txt->seq_beg = txt->seq_len + txt->seq_cnt;

   // Alloc sequence names.
   txt->seq_name = malloc(txt->seq_cnt * sizeof(char *));
   if (txt->seq_name == NULL)
      return NULL;

   // Pointers to seq_names.
   char * str_ptr = (char *) (txt->seq_beg + txt->seq_cnt);
   for (int i = 0; i < txt->seq_cnt; i++) {
      txt->seq_name[i] = str_ptr;
      str_ptr += strlen(txt->seq_name[i]) + 1;
   }

   // Pointer to text.
   txt->text = (uint8_t *) str_ptr;

   // Set memory size.
   txt->mem_txt = txt->txt_len;
   txt->mem_seq = txt->seq_cnt;

   // Set symbol alphabet.
   txt->sym = sym;

   close(fd);

   return txt;
   
 free_and_return:
   close(fd);
   txt_free(txt);
   return NULL;
}


// Aux functions.

int64_t
seq_bisect
(
  int64_t   beg,
  int64_t   end,
  int64_t   target,
  int64_t * vals
)
{
   if (end - beg < 2) 
      return beg;

   int64_t mid = (beg + end) / 2;
   if (target < vals[mid])
      return seq_bisect(beg, mid, target, vals);
   else
      return seq_bisect(mid, end, target, vals);
}
