#include <stdlib.h>
#include <string.h>

// Radix sort
#define RS_BITS 8
#define RS_MASK ((1<<RS_BITS) - 1)
#define RS_SIZE (1<<RS_BITS)
#define toDEC(a) (a/3*2+a%3-1)
#define fromDEC(a) (a*3/2+1)
#define leq2(a1,a2,b1,b2) (a1 < b1 || (a1 == b1 && a2 <= b2))
#define leq3(a1,a2,a3,b1,b2,b3) (a1 < b1 || (a1 == b1 && leq2(a2,a3,b2,b3)))
#define max(a,b) (a > b ? a : b)

long * dc3          (char* text);
void   suffixArray  (long* T, long* SA, long n, long K);
void   radixSort    (long* a, long* b, long* v, long n, long maxval, int offset);
