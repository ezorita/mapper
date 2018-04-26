#include "unittest.h"
#include "index_bwt.h"

void
test_bwt_build
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   int64_t sa_ref[] = {31,15,30,14,29,17,20,23,2,5,8,4,18,11,27,21,24,19,3,26,12,6,9,13,28,16,22,1,7,10,25,0};
   
   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   int64_t * sa_buf = malloc(50*sizeof(int64_t));
   test_assert_critical(sa_buf != NULL);

   test_assert(sar_get_range(0,32,sa_buf,sar) == 0);
   test_assert(memcmp(sa_buf, sa_ref, 32*sizeof(int64_t)) == 0);

   test_assert(bwt_build(NULL, sar) == NULL);
   test_assert(bwt_build(txt, NULL) == NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);
   
   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   // Query 'G' -> 6 hits at fp = 17, SA={19,3,26,12,6,9}.
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 6);
   test_assert(bwt_start(q) == 17);
   test_assert(bwt_depth(q) == 1);

   int64_t q00_sa[] = {19,3,26,12,6,9};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q00_sa, 6*sizeof(int64_t)) == 0);

   // Query 'GT' -> 3 hits at fp = 20, SA={12,6,9}.
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 3);
   test_assert(bwt_start(q) == 20);
   test_assert(bwt_depth(q) == 2);

   int64_t q01_sa[] = {12,6,9};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q01_sa, 3*sizeof(int64_t)) == 0);

   // Query 'GTA' -> 2 hits at fp = 20, SA={12,6}.
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2);
   test_assert(bwt_start(q) == 20);
   test_assert(bwt_depth(q) == 3);

   int64_t q02_sa[] = {12,6};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q02_sa, 2*sizeof(int64_t)) == 0);

   // Query 'GTAG' -> 1 hit at fp = 21, SA={6}.
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_start(q) == 21);
   test_assert(bwt_depth(q) == 4);

   int64_t q03_sa[] = {6};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q03_sa, 1*sizeof(int64_t)) == 0);

   // Query 'GTAGC' -> 0 hits.
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 0);

   free(q);

   // New query.
   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   // Query 'T' -> 9 hits at fp = 23, SA={13,28,16,22,1,7,10,25,0}.
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 9);
   test_assert(bwt_start(q) == 23);
   test_assert(bwt_depth(q) == 1);

   int64_t q10_sa[] = {13,28,16,22,1,7,10,25,0};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q10_sa, 9*sizeof(int64_t)) == 0);

   // Query 'TA' -> 6 hits at fp = 23, SA={13,28,16,22,1,7}.
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 6);
   test_assert(bwt_start(q) == 23);
   test_assert(bwt_depth(q) == 2);

   int64_t q11_sa[] = {13,28,16,22,1,7};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q11_sa, 6*sizeof(int64_t)) == 0);

   // Query 'TAA' -> 1 hits at fp = 24, SA={28}.
   bwtquery_t * q2 = bwt_dup_query(q);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q2, q2) == 0);
   test_assert(bwt_size(q2) == 1);
   test_assert(bwt_start(q2) == 24);
   test_assert(bwt_depth(q2) == 3);

   int64_t q12_sa[] = {28};
   test_assert(sar_get_range(bwt_start(q2), bwt_size(q2), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q12_sa, 1*sizeof(int64_t)) == 0);
   free(q2);

   // Query 'CTA' -> 2 hits at fp = 14, SA={27,21}.
   bwtquery_t * q3 = bwt_dup_query(q);
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_PREFIX, q3, q3) == 0);
   test_assert(bwt_size(q3) == 2);
   test_assert(bwt_start(q3) == 14);
   test_assert(bwt_depth(q3) == 3);

   int64_t q13_sa[] = {27,21};
   test_assert(sar_get_range(bwt_start(q3), bwt_size(q3), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q13_sa, 2*sizeof(int64_t)) == 0);
   free(q3);

   // Query 'XTA' -> 2 hits at fp = 14, SA={27,21}.
   bwtquery_t ** qv = bwt_new_vec(bwt);
   test_assert(bwt_query_all(BWT_QUERY_PREFIX, q, qv) == 0);
   // 'ATA' -> 0 hits.
   test_assert(bwt_size(qv[0]) == 0);
   // 'CTA' -> 2 hits.
   test_assert(bwt_size(qv[1]) == 2);
   test_assert(bwt_start(qv[1]) == 14);
   // 'GTA' -> 2 hits.
   test_assert(bwt_size(qv[2]) == 2);
   test_assert(bwt_start(qv[2]) == 20);
   // 'TTA' -> 2 hits.
   test_assert(bwt_size(qv[3]) == 1);
   test_assert(bwt_start(qv[3]) == 31);
   // 'NTA' -> 0 hits.
   test_assert(bwt_size(qv[4]) == 0);

   bwt_free_vec(qv);

   // Query 'TAC' -> 2 hits at fp=25, SA={16,22}.
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2);
   test_assert(bwt_start(q) == 25);
   test_assert(bwt_depth(q) == 3);

   int64_t q14_sa[] = {16,22};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q14_sa, 2*sizeof(int64_t)) == 0);

   // Query 'TACX'.
   qv = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);
   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, q, qv) == 0);

   // 'TACA' -> 0 hits.
   test_assert(bwt_size(qv[0]) == 0);
   // 'TACC' -> 0 hits.
   test_assert(bwt_size(qv[1]) == 0);
   // 'TACG' -> 1 hits.
   test_assert(bwt_size(qv[2]) == 1);
   test_assert(bwt_start(qv[2]) == 25);
   test_assert(bwt_depth(qv[2]) == 4);
   // 'TACT' -> 1 hits.
   test_assert(bwt_size(qv[3]) == 1);
   test_assert(bwt_start(qv[3]) == 26);
   test_assert(bwt_depth(qv[3]) == 4);
   // 'TACN' -> 0 hits.
   test_assert(bwt_size(qv[4]) == 0);

   free(q);
   bwt_free_vec(qv);

   // Prefix query.
   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   
   // Query 'C' -> 6 hits, fp=11, SA={4,18,11,27,21,24}.
   test_assert(bwt_prefix(sym_index('C',sym), q, q) == 0);
   test_assert(bwt_size(q) == 6);
   test_assert(bwt_start(q) == 11);
   test_assert(bwt_depth(q) == 1);
   
   int64_t q20_sa[] = {4,18,11,27,21,24};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q20_sa, 6*sizeof(int64_t)) == 0);

   // Query 'GC' -> 2 hits, fp=18, SA={3,26}.
   test_assert(bwt_prefix(sym_index('G',sym), q, q) == 0);
   test_assert(bwt_size(q) == 2);
   test_assert(bwt_start(q) == 18);
   test_assert(bwt_depth(q) == 2);
   
   int64_t q21_sa[] = {3,26};
   test_assert(sar_get_range(bwt_start(q), bwt_size(q), sa_buf, sar) == 0);
   test_assert(memcmp(sa_buf, q21_sa, 2*sizeof(int64_t)) == 0);

   // Query 'XGC'.
   qv = bwt_new_vec(bwt);
   test_assert(bwt_prefix_all(q, qv) == 0);
   // 'AGC' -> 1 hits.
   test_assert(bwt_size(qv[0]) == 1);
   test_assert(bwt_start(qv[0]) == 8);
   test_assert(bwt_depth(qv[0]) == 3);
   // 'CGC' -> 0 hits.
   test_assert(bwt_size(qv[1]) == 0);
   // 'GGC' -> 0 hits.
   test_assert(bwt_size(qv[2]) == 0);
   // 'TGC' -> 1 hits.
   test_assert(bwt_size(qv[3]) == 1);
   test_assert(bwt_start(qv[3]) == 30);
   test_assert(bwt_depth(qv[3]) == 3);
   // 'NGC' -> 0 hits.
   test_assert(bwt_size(qv[4]) == 0);

   // Duplicate vector.
   bwtquery_t ** qv2 = bwt_dup_vec(qv);

   // Query 'TAGC' -> 1 hit at fp=27, SA={1}.
   bwtquery_t * q4 = qv2[0];
   test_assert(bwt_prefix(sym_index('T', sym), q4, q4) == 0);
   test_assert(bwt_size(q4) == 1);
   test_assert(bwt_start(q4) == 27);
   test_assert(bwt_depth(q4) == 4);
   test_assert(sar_get(bwt_start(q4), sar) == 1);

   free(q);
   bwt_free_vec(qv);
   bwt_free_vec(qv2);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);

   // Test with long text.
   txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("ACGTTGCA", txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }
   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("TGCAACGT", txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }

   sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
  
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 4000);
   test_assert(bwt_start(q) == 2000);
   test_assert(bwt_depth(q) == 1);
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 2);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 3);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 4);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 5);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 6);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 0);

   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
   free(sa_buf);
}

void
test_bwt_new_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t * q = bwt_new_query(NULL);
   test_assert(q == NULL);
   
   q = bwt_new_query(bwt);   
   test_assert_critical(q != NULL);
   test_assert(bwt_start(q) == 0);
   test_assert(bwt_size(q) == 32);
   test_assert(bwt_depth(q) == 0);

   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_dup_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t * q = bwt_new_query(NULL);
   test_assert(q == NULL);
   
   q = bwt_new_query(bwt);   
   test_assert_critical(q != NULL);
   test_assert(bwt_start(q) == 0);
   test_assert(bwt_size(q) == 32);
   test_assert(bwt_depth(q) == 0);

   bwtquery_t * q2 = bwt_dup_query(NULL);
   test_assert(q2 == NULL);

   q2 = bwt_dup_query(q);
   test_assert_critical(q2 != NULL);
   test_assert(bwt_start(q2) == 0);
   test_assert(bwt_size(q2) == 32);
   test_assert(bwt_depth(q2) == 0);

   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 6);
   test_assert(bwt_start(q) == 17);
   test_assert(bwt_depth(q) == 1);

   test_assert(bwt_start(q2) == 0);
   test_assert(bwt_size(q2) == 32);
   test_assert(bwt_depth(q2) == 0);
   free(q2);
   q2 = bwt_dup_query(q);
   test_assert(bwt_size(q2) == 6);
   test_assert(bwt_start(q2) == 17);
   test_assert(bwt_depth(q2) == 1);

   free(q);
   free(q2);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_new_vec
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t ** qv = bwt_new_vec(NULL);
   test_assert(qv == NULL);
   
   qv = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);
   test_assert(bwt_start(qv[0]) == 0);
   test_assert(bwt_size(qv[0]) == 32);
   test_assert(bwt_depth(qv[0]) == 0);
   test_assert(bwt_start(qv[1]) == 0);
   test_assert(bwt_size(qv[1]) == 32);
   test_assert(bwt_depth(qv[1]) == 0);
   test_assert(bwt_start(qv[2]) == 0);
   test_assert(bwt_size(qv[2]) == 32);
   test_assert(bwt_depth(qv[2]) == 0);
   test_assert(bwt_start(qv[3]) == 0);
   test_assert(bwt_size(qv[3]) == 32);
   test_assert(bwt_depth(qv[3]) == 0);
   test_assert(bwt_start(qv[4]) == 0);
   test_assert(bwt_size(qv[4]) == 32);
   test_assert(bwt_depth(qv[4]) == 0);

   bwt_free_vec(qv);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_dup_vec
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t ** qv  = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);
   bwtquery_t ** qv2 = bwt_dup_vec(qv);
   test_assert_critical(qv2 != NULL);
   test_assert(bwt_start(qv2[0]) == 0);
   test_assert(bwt_size(qv2[0]) == 32);
   test_assert(bwt_depth(qv2[0]) == 0);
   test_assert(bwt_start(qv2[1]) == 0);
   test_assert(bwt_size(qv2[1]) == 32);
   test_assert(bwt_depth(qv2[1]) == 0);
   test_assert(bwt_start(qv2[2]) == 0);
   test_assert(bwt_size(qv2[2]) == 32);
   test_assert(bwt_depth(qv2[2]) == 0);
   test_assert(bwt_start(qv2[3]) == 0);
   test_assert(bwt_size(qv2[3]) == 32);
   test_assert(bwt_depth(qv2[3]) == 0);
   test_assert(bwt_start(qv2[4]) == 0);
   test_assert(bwt_size(qv2[4]) == 32);
   test_assert(bwt_depth(qv2[4]) == 0);

   bwt_query_all(BWT_QUERY_SUFFIX, qv[0], qv);
   test_assert(bwt_start(qv2[0]) == 0);
   test_assert(bwt_size(qv2[0]) == 32);
   test_assert(bwt_depth(qv2[0]) == 0);
   test_assert(bwt_start(qv2[1]) == 0);
   test_assert(bwt_size(qv2[1]) == 32);
   test_assert(bwt_depth(qv2[1]) == 0);
   test_assert(bwt_start(qv2[2]) == 0);
   test_assert(bwt_size(qv2[2]) == 32);
   test_assert(bwt_depth(qv2[2]) == 0);
   test_assert(bwt_start(qv2[3]) == 0);
   test_assert(bwt_size(qv2[3]) == 32);
   test_assert(bwt_depth(qv2[3]) == 0);
   test_assert(bwt_start(qv2[4]) == 0);
   test_assert(bwt_size(qv2[4]) == 32);
   test_assert(bwt_depth(qv2[4]) == 0);
   
   bwt_free_vec(qv2);
   qv2 = bwt_dup_vec(qv);
   test_assert_critical(qv2 != NULL);
   test_assert(bwt_start(qv2[0]) == 2);
   test_assert(bwt_size(qv2[0]) == 9);
   test_assert(bwt_depth(qv2[0]) == 1);
   test_assert(bwt_start(qv2[1]) == 11);
   test_assert(bwt_size(qv2[1]) == 6);
   test_assert(bwt_depth(qv2[1]) == 1);
   test_assert(bwt_start(qv2[2]) == 17);
   test_assert(bwt_size(qv2[2]) == 6);
   test_assert(bwt_depth(qv2[2]) == 1);
   test_assert(bwt_start(qv2[3]) == 23);
   test_assert(bwt_size(qv2[3]) == 9);
   test_assert(bwt_depth(qv2[3]) == 1);
   test_assert(bwt_size(qv2[4]) == 0);

   bwt_free_vec(qv);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_query
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   test_assert(bwt_query(-1, BWT_QUERY_SUFFIX, q, q) == -1);
   test_assert(bwt_query(sym_count(sym), BWT_QUERY_SUFFIX, q, q) == -1);
   test_assert(bwt_query(sym_index('A',sym), 3498, q, q) == -1);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, NULL, q) == -1);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, NULL) == -1);
   test_assert(bwt_start(q) == 0);
   test_assert(bwt_size(q) == 32);
   test_assert(bwt_depth(q) == 0);


   bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q);
   test_assert(bwt_start(q) == 2);
   test_assert(bwt_size(q) == 9);
   test_assert(bwt_depth(q) == 1);

   bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q);
   test_assert(bwt_start(q) == 4);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 2);

   bwt_query(sym_index('T',sym), BWT_QUERY_PREFIX, q, q);
   test_assert(bwt_start(q) == 24);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 3);
   
   bwt_query(sym_index('C',sym), BWT_QUERY_PREFIX, q, q);
   test_assert(bwt_start(q) == 14);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 4);

   free(q);
   
   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);

   bwt_query(sym_index('A',sym), BWT_QUERY_PREFIX, q, q);
   test_assert(bwt_start(q) == 2);
   test_assert(bwt_size(q) == 9);
   test_assert(bwt_depth(q) == 1);

   bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q);
   test_assert(bwt_start(q) == 5);
   test_assert(bwt_size(q) == 3);
   test_assert(bwt_depth(q) == 2);

   bwt_query(sym_index('T',sym), BWT_QUERY_PREFIX, q, q);
   test_assert(bwt_start(q) == 25);
   test_assert(bwt_size(q) == 2);
   test_assert(bwt_depth(q) == 3);

   bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q);
   test_assert(bwt_start(q) == 26);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 4);

   bwt_query(sym_index('C',sym), BWT_QUERY_PREFIX, q, q);
   test_assert(bwt_start(q) == 15);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 5);

   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_query_all
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t ** q0 = bwt_new_vec(bwt);
   test_assert_critical(q0 != NULL);

   test_assert(bwt_query_all(1023, q0[0], q0) == -1);
   test_assert(bwt_query_all(-1, q0[0], q0) == -1);
   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, NULL, q0) == -1);
   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, q0[0], NULL) == -1);
   free(q0[2]);
   q0[2] = NULL;
   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, q0[0], q0) == -1);
   bwt_free_vec(q0);

   q0 = bwt_new_vec(bwt);
   test_assert_critical(q0 != NULL);
   
   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, q0[0], q0) == 0);
   test_assert(bwt_start(q0[0]) == 2);
   test_assert(bwt_size(q0[0]) == 9);
   test_assert(bwt_depth(q0[0]) == 1);
   test_assert(bwt_start(q0[1]) == 11);
   test_assert(bwt_size(q0[1]) == 6);
   test_assert(bwt_depth(q0[1]) == 1);
   test_assert(bwt_start(q0[2]) == 17);
   test_assert(bwt_size(q0[2]) == 6);
   test_assert(bwt_depth(q0[2]) == 1);
   test_assert(bwt_start(q0[3]) == 23);
   test_assert(bwt_size(q0[3]) == 9);
   test_assert(bwt_depth(q0[3]) == 1);
   test_assert(bwt_size(q0[4]) == 0);

   bwtquery_t ** q1 = bwt_new_vec(bwt);
   test_assert_critical(q1 != NULL);

   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, q0[0], q1) == 0);
   test_assert(bwt_start(q1[0]) == 4);
   test_assert(bwt_size(q1[0]) == 1);
   test_assert(bwt_depth(q1[0]) == 2);
   test_assert(bwt_start(q1[1]) == 5);
   test_assert(bwt_size(q1[1]) == 3);
   test_assert(bwt_depth(q1[1]) == 2);
   test_assert(bwt_start(q1[2]) == 8);
   test_assert(bwt_size(q1[2]) == 3);
   test_assert(bwt_depth(q1[2]) == 2);
   test_assert(bwt_size(q1[3]) == 0);
   test_assert(bwt_size(q1[4]) == 0);

   bwtquery_t ** q2 = bwt_new_vec(bwt);
   test_assert_critical(q2 != NULL);

   test_assert(bwt_query_all(BWT_QUERY_SUFFIX, q1[2], q2) == 0);
   test_assert(bwt_size(q2[0]) == 0);
   test_assert(bwt_start(q2[1]) == 8);
   test_assert(bwt_size(q2[1]) == 1);
   test_assert(bwt_depth(q2[1]) == 3);
   test_assert(bwt_size(q2[2]) == 0);
   test_assert(bwt_start(q2[3]) == 9);
   test_assert(bwt_size(q2[3]) == 2);
   test_assert(bwt_depth(q2[1]) == 3);
   test_assert(bwt_size(q2[4]) == 0);
   
   bwt_free_vec(q0);
   bwt_free_vec(q1);
   bwt_free_vec(q2);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_prefix
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   test_assert(bwt_start(q) == 0);
   test_assert(bwt_size(q) == 32);
   test_assert(bwt_depth(q) == 0);
   
   test_assert(bwt_prefix(sym_index('A',sym), NULL, q) == -1);
   test_assert(bwt_prefix(sym_index('A',sym), q, NULL) == -1);
   test_assert(bwt_prefix(sym_index('A',sym), q, NULL) == -1);
   test_assert(bwt_prefix(-1, q, q) == -1);
   test_assert(bwt_prefix(sym_count(sym), q, q) == -1);
   test_assert(bwt_start(q) == 0);
   test_assert(bwt_size(q) == 32);
   test_assert(bwt_depth(q) == 0);

   test_assert(bwt_prefix(sym_index('a',sym), q, q) == 0);
   test_assert(bwt_start(q) == 2);
   test_assert(bwt_size(q) == 9);
   test_assert(bwt_depth(q) == 1);
   test_assert(bwt_prefix(sym_index('T',sym), q, q) == 0);
   test_assert(bwt_start(q) == 23);
   test_assert(bwt_size(q) == 6);
   test_assert(bwt_depth(q) == 2);
   test_assert(bwt_prefix(sym_index('C',sym), q, q) == 0);
   test_assert(bwt_start(q) == 14);
   test_assert(bwt_size(q) == 2);
   test_assert(bwt_depth(q) == 3);
   test_assert(bwt_prefix(sym_index('G',sym), q, q) == 0);
   test_assert(bwt_start(q) == 19);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 4);
   test_assert(bwt_prefix(sym_index('T',sym), q, q) == 0);
   test_assert(bwt_start(q) == 30);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 5);
   test_assert(bwt_prefix(sym_index('C',sym), q, q) == 0);
   test_assert(bwt_start(q) == 16);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 6);

   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}

void
test_bwt_prefix_all
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t ** qv = bwt_new_vec(bwt);
   test_assert_critical(qv != NULL);

   test_assert(bwt_prefix_all(NULL, NULL) == -1);
   test_assert(bwt_prefix_all(NULL, qv) == -1);
   test_assert(bwt_prefix_all(qv[0], NULL) == -1);
   free(qv[1]);
   qv[1] = NULL;
   test_assert(bwt_prefix_all(qv[0], qv) == -1);
   
   bwt_free_vec(qv);

   qv = bwt_new_vec(bwt);
   test_assert(bwt_prefix_all(qv[0], qv) == 0);
   test_assert(bwt_start(qv[0]) == 2);
   test_assert(bwt_size(qv[0]) == 9);
   test_assert(bwt_start(qv[1]) == 11);
   test_assert(bwt_size(qv[1]) == 6);
   test_assert(bwt_start(qv[2]) == 17);
   test_assert(bwt_size(qv[2]) == 6);
   test_assert(bwt_start(qv[3]) == 23);
   test_assert(bwt_size(qv[3]) == 9);
   test_assert(bwt_size(qv[4]) == 0);

   test_assert(bwt_prefix_all(qv[1], qv) == 0);
   test_assert(bwt_start(qv[0]) == 5);
   test_assert(bwt_size(qv[0]) == 3);
   test_assert(bwt_size(qv[1]) == 0);
   test_assert(bwt_start(qv[2]) == 18);
   test_assert(bwt_size(qv[2]) == 2);
   test_assert(bwt_start(qv[3]) == 29);
   test_assert(bwt_size(qv[3]) == 1);
   test_assert(bwt_size(qv[4]) == 0);

   bwt_free_vec(qv);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}


void
test_bwt_get_text
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   test_assert(bwt_get_text(NULL) == NULL);
   test_assert(bwt_get_text(bwt) == txt);

   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);
}


void
test_bwt_helpers
(void)
{
   sym_t * sym = sym_new_dna();
   test_assert_critical(sym != NULL);

   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   test_assert(txt_append("TTAGCAGTAGTCGTA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);
   test_assert(txt_append("TACGACTACTGCTAA", txt) == 0);
   test_assert(txt_append_wildcard(txt) == 0);

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   
   test_assert(bwt_start(NULL) == -1);
   test_assert(bwt_size(NULL) == -1);
   test_assert(bwt_depth(NULL) == -1);
   test_assert(bwt_start(q) == 0);
   test_assert(bwt_size(q) == 32);
   test_assert(bwt_depth(q) == 0);
   
   test_assert(bwt_query(0, BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_start(q) == 2);
   test_assert(bwt_size(q) == 9);
   test_assert(bwt_depth(q) == 1);

   test_assert(bwt_query(0, BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_start(q) == 4);
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 2);
   
   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);

}

void
test_bwt_file
(void)
{
   sym_t * sym = sym_new_dna();
   // Test with long text.
   txt_t * txt = txt_new(sym);
   test_assert_critical(txt != NULL);
   // Bidirectional BWT serach requires FW and RC text.
   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("ACGTTGCA", txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }
   for (int i = 0; i < 1000; i++) {
      test_assert(txt_append("TGCAACGT", txt) == 0);
      test_assert(txt_append_wildcard(txt) == 0);
   }

   sar_t * sar = sar_build(txt);
   test_assert_critical(sar != NULL);

   bwt_t * bwt = bwt_build(txt, sar);
   test_assert_critical(bwt != NULL);

   bwtquery_t * q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 4000);
   test_assert(bwt_start(q) == 2000);
   test_assert(bwt_depth(q) == 1);
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 2);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 3);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 4);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 5);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 6);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 0);

   test_assert(bwt_file_write("test00.bwt", bwt) == 0);
   free(q);
   bwt_free(bwt);
   
   bwt = bwt_file_read("test00.bwt", txt);
   test_assert_critical(bwt != NULL);

   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 4000);
   test_assert(bwt_start(q) == 2000);
   test_assert(bwt_depth(q) == 1);
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 2);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 3);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 4);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 5);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 6);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 0);

   free(q);
   bwt_free(bwt);

   // Repeat with different mark interval.

   bwt = bwt_build_opt(txt, sar, 32);
   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 4000);
   test_assert(bwt_start(q) == 2000);
   test_assert(bwt_depth(q) == 1);
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 2);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 3);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 4);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 5);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 6);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 0);

   test_assert(bwt_file_write("test01.bwt", bwt) == 0);
   free(q);
   bwt_free(bwt);
   
   bwt = bwt_file_read("test01.bwt", txt);
   test_assert_critical(bwt != NULL);

   q = bwt_new_query(bwt);
   test_assert_critical(q != NULL);
   
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 4000);
   test_assert(bwt_start(q) == 2000);
   test_assert(bwt_depth(q) == 1);
   test_assert(bwt_query(sym_index('C',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 2);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 3);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 2000);
   test_assert(bwt_start(q) == 4000);
   test_assert(bwt_depth(q) == 4);
   test_assert(bwt_query(sym_index('T',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 5);
   test_assert(bwt_query(sym_index('G',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 1000);
   test_assert(bwt_start(q) == 5000);
   test_assert(bwt_depth(q) == 6);
   test_assert(bwt_query(sym_index('A',sym), BWT_QUERY_SUFFIX, q, q) == 0);
   test_assert(bwt_size(q) == 0);

   free(q);
   bwt_free(bwt);
   sar_free(sar);
   txt_free(txt);
   sym_free(sym);

}

// Define test cases to be run (for export).
const test_case_t test_cases_index_bwt[] = {
   {"index_bwt/bwt_build",            test_bwt_build},
   {"index_bwt/bwt_new_query",        test_bwt_new_query},
   {"index_bwt/bwt_dup_query",        test_bwt_dup_query},
   {"index_bwt/bwt_new_vec",          test_bwt_new_vec},
   {"index_bwt/bwt_dup_vec",          test_bwt_dup_vec},
   {"index_bwt/bwt_query",            test_bwt_query},
   {"index_bwt/bwt_query_all",        test_bwt_query_all},
   {"index_bwt/bwt_prefix",           test_bwt_prefix},
   {"index_bwt/bwt_prefix_all",       test_bwt_prefix_all},
   {"index_bwt/bwt_get_text",         test_bwt_get_text},
   {"index_bwt/bwt_helpers",          test_bwt_helpers},
   {"index_bwt/bwt_file",             test_bwt_file},
   {NULL, NULL}, // Sentinel. //
};
