#include "io.h"

struct iostream_t {
   FILE     *  fin;
   int64_t     line;
   size_t      bytes_cnt;
   size_t      bytes_max;
   int         eof;
   seqread_t * (*parser)(char **);
   size_t      parser_lines;
};

// Private function headers.
int     io_stream_set_parser  (iostream_t * stream);


// Interface function sources.

iostream_t *
io_stream_open
(
  char    * path
)
{
   // Check arguments.
   error_test_msg(path == NULL, "argument 'path' is NULL.");

   iostream_t * stream = io_stream_open_buf(path, IO_DEFAULT_BUFSIZE);
   error_test(stream == NULL);

   return stream;

 failure_return:
   return NULL;
}


iostream_t *
io_stream_open_buf
(
 char    * path,
 size_t    bufsize
)
{
   // Declare variables.
   iostream_t  * stream = NULL;
   FILE        * fin    = NULL;

   // Check arguments.
   error_test_msg(path == NULL, "argument 'path' is NULL.");

   // Open file for reading.
   fin = fopen(path, "r");
   error_test_def(fin == NULL);

   // Alloc stream and set values.
   stream = malloc(sizeof(iostream_t));
   error_test_mem(stream);

   stream->fin       = fin;
   stream->line      = 1;
   stream->bytes_cnt = 0;
   stream->bytes_max = bufsize;
   stream->eof       = 0;
   stream->parser    = NULL;

   return stream;

 failure_return:
   if (fin != NULL)
      fclose(fin);
   free(stream);
   return NULL;
}

int
io_stream_close
(
  iostream_t  * stream
)
{
   // Check arguments.
   error_test_msg(stream == NULL, "argument 'stream' is NULL.");
   
   // Close stream.
   if (stream->fin != NULL) {
      error_test_def(fclose(stream->fin));
      stream->fin = NULL;
   }

   return 0;

 failure_return:
   stream->fin = NULL;
   return -1;
}

void
io_stream_free
(
  iostream_t * stream
)
{
   if (stream != NULL) {
      io_stream_close(stream);
      free(stream);
   }

   return;
}


gstack_t *
io_stream_read_seq
(
  iostream_t * stream
)
{
   // Declare variables.
   char      * line  = NULL;
   char     ** lines = NULL;
   gstack_t  * stack = NULL;
   
   // Check arguments.
   error_test_msg(stream == NULL, "argument 'stream' is NULL.");

   // Don't read more if EOF was found before.
   if (stream->eof == 1) {
      return NULL;
   }
   
   // Determine parse function.
   if (stream->parser == NULL) {
      error_test(io_stream_set_parser(stream) == -1);
   }

   // Alloc stack.
   stack = seqread_stack(128);
   error_test(stack == NULL);

   // Alloc memory to read lines.
   size_t linesz;
   ssize_t bytes;
   size_t lineno = 0;
   int buflineno = stream->parser_lines+1;
   
   // Alloc line buffers.
   lines = calloc(buflineno, sizeof(char *));
   error_test_mem(lines);
   
   // Fill up to max memory (stream->bytes_max).
   stream->bytes_cnt = 0;
   while ((bytes = getline(&line, &linesz, stream->fin)) > 0) {
      // Remove newline character.
      if (line[bytes-1] == '\n') {
	 line[--bytes] = 0;
      }
      // Format line
      size_t fline = lineno%stream->parser_lines;
      // Read stream->parser_lines.
      lines[fline] = strdup(line);
      error_test_mem(lines[fline]);
      // Send to parser.
      if (fline == stream->parser_lines - 1) {
	 seqread_t * read = NULL;
	 if (stream->parser_lines > 1) {
	    // Parse line.
	    read = stream->parser(lines);
	    error_test(read == NULL);
	 } else {
	    // Generate tag with line number.
	    char * tag = malloc(20);
	    error_test_mem(tag);
	    sprintf(tag, ">%ld", stream->line);
	    // Generate fake fasta format.
	    lines[1] = lines[0];
	    lines[0] = tag;
	    // Parse line.
	    read = stream->parser(lines);
	    error_test(read == NULL);
	 }
	 // Store in stack.
	 error_test(gstack_push(read, stack) == -1);
	 // Free strdups.
	 for (int i = 0; i < buflineno; i++) {
	    free(lines[i]);
	    lines[i] = NULL;
	 }
      }
      // Update counters.
      stream->bytes_cnt += bytes;
      stream->line++;
      lineno++;

      // Break conditions: reached max mem and do not split a read!
      if (stream->bytes_cnt >= stream->bytes_max && lineno%stream->parser_lines == 0) {
	 break;
      }
   }

   if (bytes == -1) {
      stream->eof = 1;
   }

   // Free mem.
   // Free strdups.
   if (lines != NULL) {
      for (int i = 0; i < buflineno; i++) {
	 free(lines[i]);
      }
   }
   free(line);
   free(lines);

   return stack;

 failure_return:
   if (stream != NULL && lines != NULL) {
      for (int i = 0; i < stream->parser_lines; i++) {
	 free(lines[i]);
      }
   }
   free(line);
   free(lines);
   gstack_free(stack);
   return NULL;
}


// Private functions.

int
io_stream_set_parser
(
  iostream_t * stream
)
{
   // Check arguments.
   error_test_msg(stream == NULL, "argument 'stream' is NULL.");
   error_test_msg(stream->fin == NULL, "stream is not open.");
   
   // Check the first character and determine the parsing function to use.
   char fchar = fgetc(stream->fin);
   switch (fchar) {
   case '@': {
      stream->parser       = seqread_parse_fastq;
      stream->parser_lines = 4;
      break;
   }
   case '>': {
      stream->parser       = seqread_parse_fasta;
      stream->parser_lines = 2;
      break;	 
   }
   default:
      stream->parser       = seqread_parse_fasta;
      stream->parser_lines = 1;
      break;
   }
   ungetc(fchar, stream->fin);

   return 0;

 failure_return:
   return -1;
}
