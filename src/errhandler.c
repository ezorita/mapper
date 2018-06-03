#include "errhandler.h"

void
warning
(
 const char  * msg,
 const char  * file,
 const char  * func,
 const int     line
)
{
   if (msg == NULL) 
      fprintf(stderr, "%s\t(%s:%d)\n", func, file, line);      
   else
      fprintf(stderr, "%s\t(%s:%d) error: %s\n", func, file, line, msg);
}
