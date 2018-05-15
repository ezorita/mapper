#include "index_sar.h"

// Interface data types.
struct sar_t {
   size_t     mmap_len;
   void     * mmap_ptr;
   int64_t    sar_bits;
   int64_t    sar_len;
   int64_t  * sa;
};


// Private function headers.
int     compact_array   (sar_t * sar, uint64_t len);


// Interface functions source.
sar_t *
sar_build
(
  txt_t  * txt
)
{
   if (txt == NULL)
      return NULL;

   // Alloc memory.
   char      * data = NULL;
   saidx64_t * sa   = NULL;
   sar_t     * sar  = NULL;

   data = malloc(txt_length(txt)*sizeof(uint8_t));
   if (data == NULL)
      goto free_and_return;

   sa   = malloc(txt_length(txt)*sizeof(saidx64_t));
   if (sa == NULL)
       goto free_and_return;

   sar  = malloc(sizeof(sar_t));
   if (sar == NULL)
      goto free_and_return;

   sar->sa = NULL;
   sar->mmap_len = 0;
   sar->mmap_ptr = NULL;

   // Since wildcards are encoded as sym_count(sym), shift all symbols and force
   // wildcards to be the symbol 0.
   int n_symbols = sym_count(txt_get_symbols(txt)) + 1;
   for (uint64_t i = 0; i < txt_length(txt); i++) {
      data[i] = (txt_sym(i, txt) + 1) % n_symbols;
   }

   // Compute suffix array using divsufsort algorithm.
   divsufsort((unsigned char *) data, (saidx64_t *) sa, txt_length(txt));
   
   // Compact array.
   sar->sa = (int64_t *)sa;
   compact_array(sar, txt_length(txt));

   free(data);
   return sar;

 free_and_return:
   free(data);
   free(sa);
   sar_free(sar);
   return NULL;
}


void
sar_free
(
  sar_t   * sar
)
{
   if (sar != NULL) {
      if (sar->mmap_ptr == NULL) {
         free(sar->sa);
      } else {
         munmap(sar->mmap_ptr, sar->mmap_len);
      }
      free(sar);
   }

   return;
}


int64_t
sar_get
(
 int64_t    pos,
 sar_t    * sar
)
{
   // Check arguments.
   if (sar == NULL)
      return -1;

   int64_t sar_txt_len = (64*sar->sar_len)/sar->sar_bits;
   if (pos < 0 || pos >= sar_txt_len)
      return -1;

   // Fetch word of 'sar_bits' bits in suffix array.
   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64 - sar->sar_bits);
   uint64_t bit = pos*sar->sar_bits;
   uint64_t word = bit/64;
   bit %= 64;
   // Cast words to unsigned to avoid right shift with 1's.
   if (bit + sar->sar_bits > 64)
      return ((((uint64_t)sar->sa[word] >> bit) & mask) | ((uint64_t)sar->sa[word+1] & mask) << (64-bit)) & mask;
   else
      return ((uint64_t)sar->sa[word] >> bit) & mask;
}


int
sar_get_range
(
 int64_t    beg,
 int64_t    size,
 int64_t  * vec,
 sar_t    * sar
)
{
   // Check arguments.
   if (vec == NULL || sar == NULL)
      return -1;

   int64_t sar_txt_len = (64*sar->sar_len)/sar->sar_bits;
   if (beg < 0 || size < 1 || beg + size > sar_txt_len)
      return -1;

   // Fetch words of 'sar_bits' bits and store them in vec.
   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64 - sar->sar_bits);
   uint64_t bit = beg*sar->sar_bits;
   uint64_t word = bit/64;
   // Process bits.
   bit %= 64;
   uint64_t w = sar->sa[word];
   for (uint64_t i = 0; i < size; i++) {
      if (bit + sar->sar_bits > 64) {
         vec[i] = w >> bit;
         w = sar->sa[++word];
         vec[i] |= (w & mask) << (64-bit);
         vec[i] &= mask;
      } else {
         vec[i] = (w >> bit) & mask;
      }
      bit = (bit + sar->sar_bits)%64;
      if (bit == 0) w = sar->sa[++word];
   }
   return 0;
}


// I/O functions.
int
sar_file_write
(
  char   * filename,
  sar_t  * sar
)
{
   // Check arguments.
   if (filename == NULL || sar == NULL)
      return -1;

   // Open file.
   int fd = creat(filename, 0644);
   if (fd == -1)
      return -1;

   // Write data.
   ssize_t  e_cnt = 0;
   ssize_t  b_cnt = 0;
   uint64_t magic = SAR_FILE_MAGICNO;

   // Write magic.
   if (write(fd, &magic, sizeof(uint64_t)) == -1)
      goto close_and_error;
   
   // Write sar_bits.
   if (write(fd, (int64_t *)&(sar->sar_bits), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write sar_len.
   if (write(fd, (int64_t *)&(sar->sar_len), sizeof(int64_t)) == -1)
      goto close_and_error;

   // Write Suffix Array.
   e_cnt = 0;
   do {
      b_cnt  = write(fd, (int64_t *)sar->sa + e_cnt, (sar->sar_len - e_cnt)*sizeof(int64_t));
      if (b_cnt == -1)
         goto close_and_error;
      e_cnt += b_cnt / sizeof(int64_t);
   } while (e_cnt < sar->sar_len);

   close(fd);
   return 0;

 close_and_error:
   close(fd);
   return -1;
}


sar_t *
sar_file_read
(
  char   * filename
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
   sar_t * sar = malloc(sizeof(sar_t));
   if (sar == NULL)
      goto free_and_return;
   // Set NULL pointers.
   sar->sa = NULL;

   // Get file len and mmap file.
   struct stat sb;
   fstat(fd, &sb);
   if (sb.st_size < 24)
      goto free_and_return;

   int64_t * data = mmap(NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
   if (data == NULL)
      goto free_and_return;

   sar->mmap_len = sb.st_size;
   sar->mmap_ptr = (void *) data;

   // Read file.
   // Read magic number.
   uint64_t magic = data[0];
   if (magic != SAR_FILE_MAGICNO)
      goto free_and_return;

   sar->sar_bits = data[1];
   sar->sar_len  = data[2];
   sar->sa       = data + 3;

   close(fd);

   return sar;

 free_and_return:
   close(fd);
   sar_free(sar);
   return NULL;
}


// Private functions.
int
compact_array
(
 sar_t     * sar,
 uint64_t    len
)
{
   int bits = 0;
   while (len > ((uint64_t)1 << bits)) bits++;

   uint64_t mask = ((uint64_t)0xFFFFFFFFFFFFFFFF) >> (64-bits);
   uint64_t word = 0;
   int   lastbit = 0;

   // Clear upper bits of array[0].
   int64_t * array = sar->sa;
   array[0] &= mask;

   for (uint64_t i = 0, current; i < len; i++) {
      // Save the current value.
      current = array[i];
      // Store the compact version.
      array[word] |= (current & mask) << lastbit;
      // Update bit offset.
      lastbit += bits;
      // Word is full.
      if (lastbit >= 64) {
         lastbit = lastbit - 64;
         // Complete with remainder or set to 0 (if lastbit = 0).
         // This will clear the upper bits of array.
         array[++word] = (current & mask) >> (bits - lastbit);
      }
   }
   
   // New array size.
   sar->sar_bits = bits;
   sar->sar_len  = word + (lastbit > 0);

   // Realloc array.
   sar->sa = realloc(sar->sa, sar->sar_len*sizeof(int64_t));
   if (sar->sa == NULL)
      return -1;

   return 0;
}
