#include <errno.h>
#include <stdlib.h>
#include <string.h>

#ifndef
#define ERRHANDLER_H_

// Macros to throw errors.

#define error_test_mem(x) do { \
   if ((x) == NULL) { \
      warning("memory error", __FILE__, __func__, __LINE__);     \
      goto failure_return; \
   }} while (0)

#define error_test_def(x) do { \
   if (!(x)) { \
      warning(strerror(errno), __FILE__, __func__, __LINE__);    \
      goto failure_return; \
   }} while (0)

#define error_test_msg(x,m) do { \
   if (!(x)) { \
      warning(m, __FILE__, __func__, __LINE__); \
      goto failure_return; \
   }} while (0)

// Macros to propagate errors.

#define error_test(x) do { \
   if (!(x)) { \
   warning(NULL, __FILE__, __func__, __LINE__); \
   goto failure_return; \
   }} while (0)

// Helper functions.
void    warning   (const char * msg, const char * file, const char * func, const int line);

#endif
