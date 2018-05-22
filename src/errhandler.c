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
      fprintf(stderr, "%s\t(%s:%d)\n", file, func, line);      
   else
      fprintf(stderr, "%s\t(%s:%d) error: %s\n", file, func, line, msg);
}
