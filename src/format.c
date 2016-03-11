#include "format.h"

#define FEATURE_FORMAT 1
#define SAM_FORMAT 0

int
print_and_free
(
 seq_t          seq,
 matchlist_t  * matches,
 index_t      * index,
 formatopt_t    opt,
 int            seed_stage
)
{

   // Retrieve match.
   if (matches->pos > 0) {
      match_t match = matches->match[0];
      // Check mapping quality.
      // Format output.
      int dir;
      uint64_t g_start, g_end = 0;
      if (match.ref_s >= index->size/2) {
         g_start = index->size - match.ref_e - 2;
         g_end   = index->size - match.ref_s - 2;
         dir = 1;
      } else {
         g_start = match.ref_s;
         g_end   = match.ref_e;
         dir = 0;
      }
      // Search chromosome name.
      int chrnum = bisect_search(0, index->chr->nchr-1, index->chr->start, g_start+1)-1;

#if FEATURE_FORMAT == 1
      int span = match.read_e - match.read_s + 1;
      /*
      int mapq = 0;
      double perr = 0.01;
      //ident = (int)(-10*(log10(binoms[j]) + (len-n_err[j])*log10(1-err_rate) + n_err[j]*log10(err_rate)));
      double nerr_s = span-match.second;
      double nerr_b = span-match.ident;
      double p_b = log10(binom(span,round(nerr_b))) + nerr_b*log10(perr) + match.ident*log10(1-perr);
      double p_s = log10(binom(span,round(nerr_s))) + nerr_s*log10(perr) + match.second*log10(1-perr);
      //      fprintf(stdout, "span=%d, err_b=%.2f, err_s=%.2f, log10(p_b)=%.2f, p_s=%.2f, n2=%d, p_b=%f, binom_s(%d,%f)=%f,perr^n_s=%f, p_s=%f\n",span,nerr_b, nerr_s, p_b, p_s, match.flags, binom(span,round(nerr_b))*pow(perr,nerr_b), span, round(nerr_s),binom(span,round(nerr_s)),pow(perr,nerr_s),binom(span,round(nerr_s))*pow(perr,nerr_s));
      mapq = max(0,min(-10*(match.flags + p_s - p_b),60));
      */
      fprintf(stdout, "%s\t%s\t%ld\t%c\t%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n",
              seq.tag,
              index->chr->name[chrnum],
              g_start - index->chr->start[chrnum]+1,
              dir ? '-' : '+',
              span,
              seed_stage == 1,
              seed_stage == 2,
              seed_stage == 0,
              seed_stage == 3,
              match.hits*1.0/span,
              match.s_hits*1.0/span,
              match.e_exp[0],
              match.e_exp[1],
              match.e_exp[2],
              (int)match.mapq);
   }

      // Print results.
#elif SAM_FORMAT == 1
      fprintf(stdout, "%s\t%s\t%s\t%ld\t%d\t%d\t*\t%d\t%d\t%s\t%s",
              seq.tag+1, (dir == 1 ? "16" : "0"),
              index->chr->name[chrnum],
              g_start - index->chr->start[chrnum]+1,
              (int)match.mapq,
              match.score,
              (int)match.ident,
              match.flags,
              seq.seq,
              seq.q
              );
   }
   else {
      fprintf(stdout, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s", seq.tag+1, seq.seq, seq.q);
   }

#else
      fprintf(stdout, "%s\t%d\t%d\t%s:%ld-%ld:%c\t%d\t%d\t%.2f\t%d\n",
              seq.tag,
              match.read_s+1, match.read_e+1,
              index->chr->name[chrnum],
              g_start - index->chr->start[chrnum]+1,
              g_end - index->chr->start[chrnum]+1,
              dir ? '-' : '+',
              (int)match.ident, match.second,
              match.mapq,
              match.score);
   }
   else {
      fprintf(stdout, "%s\t0\t0\t*\n", seq.tag);
   }

#endif

   // Free match.
   free(seq.seq);
   free(seq.tag);
   free(seq.q);

   return 0;
}
