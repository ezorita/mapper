#include "dc3.h"
#include <stdio.h>

long *
dc3
(
 char * text
)
{
   long len = strlen(text);
   long * values = malloc((len+3)*sizeof(long));
   for (long i = 0; i < len; i++) values[i] = (long)text[i];
   values[len] = values[len+1] = values[len+2] = 0;
   long * SA = malloc(len*sizeof(long));
   suffixArray(values, SA, len, 255);
   free(values);
   return SA;
}

void
radixSort
(
 long * a,      // Indices to sort. (may be modified)
 long * b,      // Destination buffer.
 long * v,      // Values to which the indices point.
 long   n,      // Length of a.
 long   maxval, // Maximum value pointed in v.
 int    offset  // Sort using multiple values of v, starting at offset.
)
{
   int ref = 0, new = 1, it = 0;
   long cnt[RS_SIZE];
   long prf[RS_SIZE];
   long * s[2];

   // Count iterations per value.
   while ((maxval >> (RS_BITS*it)) & RS_MASK) it++;

   s[0] = a;
   s[1] = b;
   for (int o = offset; o >= 0; o--) {
   // Iterate over all bit blocks.
      for (long j = 0; j < it; j++) {
         int shift = RS_BITS * j;
         // Reset count and prefix.
         memset(cnt, 0, RS_SIZE*sizeof(long));
         prf[0] = 0;
         // Count radix RS_BITS.
         for (long i = 0; i < n; i++) cnt[(v[s[ref][i]+o] >> shift) & RS_MASK]++;
         // Prefix.
         for (int  i = 1; i < RS_SIZE; i++) prf[i] = prf[i-1] + cnt[i-1];
         // Sorted.
         for (long i = 0; i < n; i++) s[new][prf[(v[s[ref][i]+o] >> shift) & RS_MASK]++] = s[ref][i];
         // Swap buffers.
         ref = (ref+1)%2;
         new = (new+1)%2;
      }
   }
   
   // Move data to b.
   if (ref == 0) memcpy(s[1], s[0], n*sizeof(long));
}


void
suffixArray
(
 long * T,
 long * SA,
 long   n,
 long   K
)
{
   // Initialize values and alloc buffers.
   long n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2;
   long* R = malloc(sizeof(long)*(n02 + 3));
   long* SA12 = malloc(sizeof(long)*(n02 + 3)); 
   long* R0 = malloc(sizeof(long)*n0); 
   long* SA0 = malloc(sizeof(long)*n0);

   // Initialize sequence end.
   R[n02]= R[n02+1]= R[n02+2]=0;
   SA12[n02]=SA12[n02+1]=SA12[n02+2]=0;

   // B12 = i such that i mod3 = {1,2}.
   for (long i=0, j=0; i < n+(n0-n1); i++) if (i%3 != 0) R[j++] = i;

   // Sort suffixes.
   radixSort(R, SA12, T, n02, K, 3);

   // Rank sorted suffixes.
   long name = 0, c0 = -1, c1 = -1, c2 = -1;
   for (long i = 0; i < n02; i++) {
      if (T[SA12[i]] != c0 || T[SA12[i]+1] != c1 || T[SA12[i]+2] != c2) {
         name++;
         c0 = T[SA12[i]];
         c1 = T[SA12[i]+1];
         c2 = T[SA12[i]+2];
      }
      if (SA12[i] % 3 == 1) R[SA12[i]/3] = name;
      else R[SA12[i]/3 + n0] = name;
   }

   // Recur if repeated ranks exist.
   if (name < n02) {
      suffixArray(R, SA12, n02, name);
      // Recompute ranks.
      for (long i = 0; i < n02; i++) R[SA12[i]] = i + 1;
   }
   // If ranks are OK, recompute suffix array.
   else for (long i = 0; i < n02; i++) SA12[R[i] - 1] = i;

   // Pre-sort B0 using B1 ranks.
   for (long i=0, j=0; i < n02; i++) if (SA12[i] < n0) R0[j++] = 3*SA12[i]; 
   // RadixSort B0 lexicographically.
   radixSort(R0, SA0, T, n0, K, 0);

   // Merge B12 and B0.
   for (long p=0, t=n0-n1, k=0; k < n; k++) {
      // Conversion from index in SA12 to index in T.
#define GetI() (SA12[t] < n0 ? SA12[t] * 3  + 1 : (SA12[t] - n0) * 3 + 2)

      // Indices for B12 and B0.
      long i = GetI();
      long j = SA0[p];
      
      // If SA12[t] belongs to B1 compare T[i] and R[i+1], otherwise compare T[i], T[i+1] and R[i+2].
      if (SA12[t] < n0 ?  leq2(T[i], R[SA12[t] + n0], T[j], R[j/3]) : leq3(T[i],T[i+1],R[SA12[t ]-n0 +1], T[j],T[j+1],R[j/3+n0])) {
         SA[k] = i;
         t++;
         // If B12 has reached its end.
         if (t == n02) for (k++; p < n0; p++, k++) SA[k] = SA0[p];
      }
      else {
         SA[k] = j; 
         p++;
         // If B0 has reached its end.
         if (p == n0) for (k++; t < n02; t++, k++) SA[k] = GetI();
      }
   }

   // Free variables.
   free(R);
   free(SA12);
   free(SA0);
   free(R0);
}


