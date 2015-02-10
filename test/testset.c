#include <glib.h>
#include <stdio.h>
#include <string.h>
#include "faultymalloc.h"
#include "bwmapper.h"

typedef struct {
} fixture;


// --  ERROR MESSAGE HANDLING FUNCTIONS  -- //

char ERROR_BUFFER[1024];
char OUTPUT_BUFFER[1024];
int BACKUP_FILE_DESCRIPTOR;
int BACKUP_STDOUT;

void
mute_stderr
(void)
{
   fflush(stderr);
   int fd = open("/dev/null", O_WRONLY);
   dup2(fd, STDERR_FILENO);
}

void
mute_stdout
(void)
{
   fflush(stderr);
   int fd = open("/dev/null", O_WRONLY);
   dup2(fd, STDOUT_FILENO);
}

void
redirect_stdout_to
(
   char buffer[]
)
{
   // Flush stderr, redirect to /dev/null and set buffer.
   fflush(stdout);
   int temp = open("/dev/null", O_WRONLY);
   dup2(temp, STDOUT_FILENO);
   memset(buffer, '\0', 1024 * sizeof(char));
   setvbuf(stdout, buffer, _IOFBF, 1024);
   close(temp);
   fflush(stdout);
}


void
redirect_stderr_to
(
   char buffer[]
)
{
   // Flush stderr, redirect to /dev/null and set buffer.
   fflush(stderr);
   int temp = open("/dev/null", O_WRONLY);
   dup2(temp, STDERR_FILENO);
   memset(buffer, '\0', 1024 * sizeof(char));
   setvbuf(stderr, buffer, _IOFBF, 1024);
   close(temp);
   // Fill the buffer (needed for reset).
   fprintf(stderr, "fill the buffer");
   fflush(stderr);
}

void
unredirect_sderr
(void)
{
   fflush(stderr);
   dup2(BACKUP_FILE_DESCRIPTOR, STDERR_FILENO);
   setvbuf(stderr, NULL, _IONBF, 0);
}

// Declare test functions.
// test-hitmap:
void test_combine_matches(void);
void test_feedback_gaps(void);
void test_fill_gaps(void);
void test_matchlist_new(void);
void test_matchlist_add(void);
void test_compar_seqsort(void);
void test_compar_matchid(void);
void test_compar_readstart(void);
void test_compar_readend(void);
void test_compar_matchspan(void);
void test_compar_refstart(void);

void
test_fuse_matches
(void)
{
   // Read of 500 nt.
   int slen = 500;
   hmargs_t args = {.read_ref_ratio = 2, .dist_accept = 10};
   
   matchlist_t * list = matchlist_new(5);
   match_t * match1 = malloc(sizeof(match_t));
   match_t * match2 = malloc(sizeof(match_t));
   matchlist_add(&list, match1);
   matchlist_add(&list, match2);

   // TEST 1. Forward strand.
   match1->read_s = 0;
   match1->read_e = 100;
   match1->ref_s  = 1000;
   match1->ref_e  = 1115;
   match1->dir    = 0;
   match1->ident  = 0.8;

   match2->read_s = 150;
   match2->read_e = 300;
   match2->ref_s  = 1170;
   match2->ref_e  = 1350;
   match2->dir    = 0;
   match2->ident  = 0.8;

   fuse_matches(&list, slen, args);

   g_assert(list->pos == 1);
   g_assert(list->match[0]->read_s == 0);
   g_assert(list->match[0]->read_e == 300);
   g_assert(list->match[0]->ref_s  == 1000);
   g_assert(list->match[0]->ref_e  == 1350);
   g_assert(list->match[0]->dir    == 0);
   
   //reset
   free(list->match[0]);
   list->pos = 0;
   matchlist_add(&list, match1);
   matchlist_add(&list, match2);

   // TEST 2. Reverse strand.
   match1->read_s = 100;
   match1->read_e = 170;
   match1->ref_s  = 1000;
   match1->ref_e  = 1080;
   match1->dir    = 1;

   match2->read_s = 250;
   match2->read_e = 350;
   match2->ref_s  = 760;
   match2->ref_e  = 920;
   match2->dir    = 1;

   fuse_matches(&list, slen, args);

   g_assert(list->pos == 1);
   g_assert(list->match[0]->read_s == 100);
   g_assert(list->match[0]->read_e == 350);
   g_assert(list->match[0]->ref_s  == 760);
   g_assert(list->match[0]->ref_e  == 1080);
   g_assert(list->match[0]->dir    == 1);

   //reset
   free(list->match[0]);
   list->pos = 0;
   matchlist_add(&list, match1);
   matchlist_add(&list, match2);


   // TEST 3. Different strands.
   match1->dir = 0;

   fuse_matches(&list, slen, args);

   g_assert(list->pos == 2);

   free(list);
   free(match1);
   free(match2);
}

void
test_find_repeats
(void)
{
   match_t * matches = malloc(5*sizeof(match_t));
   for (int i = 0; i < 10; i++) matches[i].repeats = matchlist_new(5);
   matchlist_t * list = matchlist_new(5);
   
   double overlap = 0.9;

   // TEST 1.
   // [0]-[1] overlap=95, span=105, 95/105 = 0.907
   // [0]-[2] overlap=92, span=110, 92/110 = 0.83
   // [1]-[2] overlap=97, span=105, 97/105 = 0.92
   matches[0].read_s = 0;
   matches[0].read_e = 100;
   matches[1].read_s = 5;
   matches[1].read_e = 105;
   matches[2].read_s = 8;
   matches[2].read_e = 110;
   matchlist_add(&list, matches);
   matchlist_add(&list, matches + 1);
   matchlist_add(&list, matches + 2);

   g_assert(find_repeats(list, overlap) == 0);
   
   g_assert(matches[0].repeats->pos == 2);
   g_assert(matches[1].repeats->pos == 3);
   g_assert(matches[2].repeats->pos == 2);

   g_assert(matches[0].repeats->match[1] == matches + 1);
   g_assert(matches[1].repeats->match[0] == matches);
   g_assert(matches[1].repeats->match[2] == matches + 2);
   g_assert(matches[2].repeats->match[0] == matches + 1);

   // TEST 2.
   // [3]-[4] 165/180 = 0.917
   // [0]-[2] 100/110 = 0.909
   // [1]-[2] 100/110 = 0.909
   matches[3].read_s = 70;
   matches[3].read_e = 250;
   matches[4].read_s = 75;
   matches[4].read_e = 240;
   matches[2].read_s = 1;
   matchlist_add(&list, matches + 3);
   matchlist_add(&list, matches + 4);

   g_assert(find_repeats(list, overlap) == 0);
   
   g_assert(matches[0].repeats->pos == 3);
   g_assert(matches[1].repeats->pos == 3);
   g_assert(matches[2].repeats->pos == 3);
   g_assert(matches[3].repeats->pos == 2);
   g_assert(matches[4].repeats->pos == 2);

   g_assert(matches[0].repeats->match[1] == matches + 2);
   g_assert(matches[0].repeats->match[2] == matches + 1);
   g_assert(matches[1].repeats->match[0] == matches);
   g_assert(matches[1].repeats->match[1] == matches + 2);
   g_assert(matches[2].repeats->match[0] == matches);
   g_assert(matches[2].repeats->match[2] == matches + 1); 
   g_assert(matches[3].repeats->match[1] == matches + 4); 
   g_assert(matches[4].repeats->match[0] == matches + 3); 

   free(matches);
   free(list);
}

void
test_combine_matches
(void)
{
   match_t * matches  = malloc(5*sizeof(match_t));
   matchlist_t * list = matchlist_new(5);
   matchlist_t * intervals;

   double overlap_tolerance = 0.1;
   
   // TEST 1. Two intervals, overlapped but less than 10%.   
   matches[0].read_s = 91; 
   matches[0].read_e = 200;
   matches[1].read_s = 0;
   matches[1].read_e = 100; // match[1]-[2] overlap 9nt. 9% of [1], < 9% of [2]. 
   matchlist_add(&list, matches);
   matchlist_add(&list, matches + 1);

   intervals = combine_matches(list, overlap_tolerance);
   g_assert(intervals->pos == 2);
   g_assert(intervals->match[0] == matches + 1);
   g_assert(intervals->match[1] == matches);
   free(intervals);

   // TEST 2. Two intervals, overlapped more than 10% wrt the smallest.
   matches[1].read_e = 102;

   intervals = combine_matches(list, overlap_tolerance);
   g_assert(intervals->pos == 1);
   g_assert(intervals->match[0] == matches);
   free(intervals);
   
   // TEST 3. Three intervals, 2 small overlapped with one big.
   matches[0].read_s = 0;
   matches[0].read_e = 100;
   matches[1].read_s = 75;
   matches[1].read_e = 125;
   matches[2].read_s = 15;
   matches[2].read_e = 70;
   matchlist_add(&list, matches + 2);

   intervals = combine_matches(list, overlap_tolerance);
   g_assert(intervals->pos == 1);
   g_assert(intervals->match[0] == matches);
   free(intervals);

   // TEST 4. Five intervals, some overlap.
   matches[1].read_s = 97;
   matches[1].read_e = 134;
   matches[4].read_s = 136;
   matches[4].read_e = 194;
   matches[2].read_s = 190;
   matches[2].read_e = 290;
   matches[3].read_s = 282;
   matches[3].read_e = 389;
   matchlist_add(&list, matches + 3);
   matchlist_add(&list, matches + 4);
   
   intervals = combine_matches(list, overlap_tolerance);
   g_assert(intervals->pos == 5);
   g_assert(intervals->match[0] == matches);
   g_assert(intervals->match[1] == matches + 1);
   g_assert(intervals->match[2] == matches + 4);
   g_assert(intervals->match[3] == matches + 2);
   g_assert(intervals->match[4] == matches + 3);
   free(intervals);
   
   free(matches);
   free(list);
   
}

void
test_feedback_gaps
(void)
// This function feeds the non-matched gaps back to the sub_t stack. The non-matched
// nucleotides will be searched again with a higher tau.
// Feedback conditions:
//  - A non-matched region longer than 'match_min_len'.
//  - A region matched with identity less than 'feedback_id_thr'.
{
   seq_t seq;
   seq.seq  = "ATGACTCGACGATCGACTAGCACGATCGAC"; // len 30.
   seq.rseq = "TACGATCGTACGACTGCTACGACGACTAGG"; // It's not true but it doesn't need to be.
   int nsubs = 2*strlen(seq.seq);
   
   // Allocate structs.
   hmargs_t args = {.kmer_size = 5, .feedback_id_thr = 0.7, .match_min_len = 5};
   matchlist_t * intervals = matchlist_new(10);
   sublist_t   * slist     = malloc(sizeof(sublist_t) + nsubs*sizeof(sub_t));
   slist->size = 0;

   // Add intervals.
   match_t * match1 = malloc(sizeof(match_t));
   match1->repeats = matchlist_new(2);
   matchlist_add(&(match1->repeats), match1);
   matchlist_add(&intervals, match1);

   // TEST 1. Nucleotides 10 to 25 covered.
   match1->read_s = 10;
   match1->read_e = 25;
   match1->ident  = 0.75;
   
   feedback_gaps(0, seq, intervals, slist, NULL, args);

   // there must be (10-5+1 + 5-5+1) * 2 = 14 subsequences.
   g_assert(slist->size == 14);
   // Check subsequence positions.
   for (int i = 0; i <= 5; i++) {
      // Check forward.
      g_assert(slist->sub[2*i].seqid == 2*i);
      // Check reverse.
      g_assert(slist->sub[2*i+1].seqid == (30 - args.kmer_size - i)*2 + 1);
   }
   g_assert(slist->sub[12].seqid == 50);
   g_assert(slist->sub[13].seqid == 1);

   // TEST 2. Both matches are long enough but match1 is below the id threshold.
   slist->size = 0;
   match_t * match2 = malloc(sizeof(match_t));
   match2->repeats = matchlist_new(2);
   matchlist_add(&(match2->repeats), match2);
   matchlist_add(&intervals, match2);

   match1->read_s = 5;
   match1->read_e = 12;
   match1->ident  = 0.68;
   match2->read_s = 15;
   match2->read_e = 26;
   match2->ident  = 0.80;
   
   feedback_gaps(0, seq, intervals, slist, NULL, args);

   // there must be (15-5+1) * 2 = 22 subsequences.
   g_assert(slist->size == 22);
   // Check subsequence positions.
   for (int i = 0; i <= 10 ; i++) {
      // Check forward.
      g_assert(slist->sub[2*i].seqid == 2*i);
      // Check reverse.
      g_assert(slist->sub[2*i+1].seqid == (30 - args.kmer_size - i)*2 + 1);
   }

   // TEST 3. Same as before, but now match1 has a repeat with identity 90%.
   slist->size = 0;
   match_t * match3 = malloc(sizeof(match_t));
   match3->repeats = matchlist_new(2);
   matchlist_add(&(match3->repeats), match3);
   matchlist_add(&(match3->repeats), match1);
   matchlist_add(&(match1->repeats), match3);

   match3->read_s = 5;
   match3->read_e = 13;
   match3->ident  = 0.9;

   feedback_gaps(0, seq, intervals, slist, NULL, args);

   // there must be (5-5+1) * 2 = 2 subsequences.
   g_assert(slist->size == 2);
   g_assert(slist->sub[0].seqid == 0);
   g_assert(slist->sub[1].seqid == 51);

   // TEST 4. Now the repeat is not long enough and triggers a subseq feedback between match3 and match2.
   slist->size = 0;
   match3->read_s = 4;
   match3->read_e = 10;
   
   feedback_gaps(0, seq, intervals, slist, NULL, args);

   // there must be ((15-10)-5+1) * 2 = 2 subsequences.
   g_assert(slist->size == 2);
   g_assert(slist->sub[0].seqid == 20);
   g_assert(slist->sub[1].seqid == 31);

   // TEST 5. Match 1 and 2 overlap, and match 1 is below id_thr, the feedback should stop where match2 starts.
   slist->size = 0;
   match1->repeats->pos = 1;
   match1->read_e = 19;
   match1->ident  = 0.69;

   feedback_gaps(0, seq, intervals, slist, NULL, args);

   // there must be (15-5+1) * 2 = 22 subsequences.
   g_assert(slist->size == 22);
   // Check subsequence positions.
   for (int i = 0; i <= 10 ; i++) {
      // Check forward.
      g_assert(slist->sub[2*i].seqid == 2*i);
      // Check reverse.
      g_assert(slist->sub[2*i+1].seqid == (30 - args.kmer_size - i)*2 + 1);
   }

   // Free memory.
   free(match1);
   free(match2);
   free(match3);
   free(intervals);
   free(slist);

}

void
test_fill_gaps
(void)
// Fills the gaps of the estimated read with greater overlap tolerance.
// For instance...
// [--*****************..........***************----] current estimate (intervals)
//               *****************                    match
//               ^^^^^^          ^                    overlap (< overlap_max_tolerance)
//                     ^^^^^^^^^^                     gap coverage (> gap_min_coverage)
// The match will be used to fill the middle gap, even though it overlaps with the other
// reads. This is done once the search is finished.
{
   // Small example with a sequence of 100 nt, allowed gaps of 10 nt.
   double overlap_max_tolerance = 0.5;
   double gap_min_coverage      = 0.9;
   int    match_minlen          = 10;
   int    seqlen                = 100;

   match_t * interv = malloc(sizeof(match_t));
   match_t * umatch = malloc(sizeof(match_t));
   
   matchlist_t * intervals = matchlist_new(5);
   matchlist_t * matches   = matchlist_new(5);

   // umatch 50 nt long. 24 samples of overlap. Gap is covered 100% - should insert.
   interv->read_s = 0;
   interv->read_e = 74;
   interv->ident  = 0.8;
   umatch->read_s = 50;
   umatch->read_e = 100;
   umatch->ident = 0.75;
   
   matchlist_add(&intervals, interv);
   matchlist_add(&matches, interv);
   matchlist_add(&matches, umatch);
   
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 2);
   g_assert(intervals->match[1] == umatch);

   // umatch 50 nt long. 26 samples of overlap (>50%). - should not insert.
   intervals->pos = 1;
   interv->read_e = 76;
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 1);
   g_assert(intervals->match[0] == interv);

   // umatch 30 nt long. 10 samples of overlap. Gap is covered 20/25 samples (<90%). - should not insert.
   umatch->read_s = 65;
   umatch->read_e = 94;
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 1);
   g_assert(intervals->match[0] == interv);

   // fill at the beginning. check final sort as well.
   interv->read_s = 25;
   interv->read_e = 100;
   umatch->read_s = 0;
   umatch->read_e = 25;
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 2);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == interv);

   // fill in between two matches.
   interv->read_s = 50;
   match_t * umatch2 = malloc(sizeof(match_t));
   umatch2->read_s = 25;
   umatch2->read_e = 50;
   umatch2->ident = 0.77;
   matchlist_add(&matches, umatch2);
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 3);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == umatch2);
   g_assert(intervals->match[2] == interv);

   // reset intervals.
   intervals->pos = 0;
   matchlist_add(&intervals, umatch);
   matchlist_add(&intervals, interv);
   
   // fill in between with overlap in both sides.
   // overlaps (exceed):
   umatch2->read_s = 12;
   umatch2->read_e = 70;
   // umatch1-umatch2 = 13 nt. (13/25 > 0.5 of umatch1) -- should not insert.
   // umatch2-interv  = 20 nt. (20/50 < 0.5 of interv)
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 2);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == interv);

   // overlaps (ok):
   umatch2->read_s = 22;
   umatch2->read_e = 53;
   // umatch1-umatch2 = 3 nt. (3/25 < 0.5 of umatch1, 3/31 < 0.5 of umatch2)
   // umatch2-interv  = 3 nt. (3/31 < 0.5 of umatch2, 3/50 < 0.5 of interv)  -- should insert.
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 3);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == umatch2);
   g_assert(intervals->match[2] == interv);

   // reset intervals.
   intervals->pos = 0;
   matchlist_add(&intervals, umatch);
   matchlist_add(&intervals, interv);

   // fill in between with insufficient gap coverage.
   umatch2->read_s = 30;
   umatch2->read_e = 60;
   // umatch2-interv = 10 nt (10/30 > 0.5)
   // gap coverage   = 20 nt (20/25 < 0.9) -- should not insert.
   g_assert(intervals->pos == 2);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == interv);

   // fill with two candidates with different identity.
   match_t * umatch2b = malloc(sizeof(match_t));
   matchlist_add(&matches, umatch2b);

   umatch2b->read_s = umatch2->read_s = 25;
   umatch2b->read_e = umatch2->read_e = 50;
   umatch2b->ident = 0.8;

   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 3);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == umatch2b);
   g_assert(intervals->match[2] == interv);

   // reset intervals.
   intervals->pos = 0;
   matchlist_add(&intervals, umatch);
   matchlist_add(&intervals, interv);

   // now umatch2.id > umatch2b.id
   umatch2->ident = 0.85;
   g_assert(fill_gaps(&intervals, matches, seqlen, match_minlen, gap_min_coverage, overlap_max_tolerance) == 0);
   g_assert(intervals->pos == 3);
   g_assert(intervals->match[0] == umatch);
   g_assert(intervals->match[1] == umatch2);
   g_assert(intervals->match[2] == interv);

   // free structures.
   free(interv);
   free(umatch);
   free(umatch2);
   free(umatch2b);
   free(matches);
   free(intervals);

}

void
test_matchlist_new
(void)
{
   matchlist_t * ml;
   
   // List with 100 elements.
   ml = matchlist_new(100);
   g_assert(ml != NULL);
   g_assert(ml->pos == 0);
   g_assert(ml->size == 100);
   free(ml);

   // List with '-1' element.
   ml = matchlist_new(-1);
   g_assert(ml != NULL);
   g_assert(ml->pos == 0);
   g_assert(ml->size == 1);
   free(ml);
}

void
test_matchlist_add
(void)
{
   matchlist_t * ml = matchlist_new(2);
   match_t     * m0 = malloc(sizeof(match_t));
   match_t     * m1 = malloc(sizeof(match_t));   
   // Insert matches.
   g_assert_cmpint(matchlist_add(&ml, m0), ==, 0);
   g_assert(ml->match[0] == m0);
   g_assert(ml->size == 2);
   g_assert(ml->pos  == 1);

   g_assert_cmpint(matchlist_add(&ml, m1), ==, 0);
   g_assert(ml->match[1] == m1);
   g_assert(ml->size == 2);
   g_assert(ml->pos  == 2);

   // Force realloc.
   g_assert_cmpint(matchlist_add(&ml, m1), ==, 0);
   g_assert(ml->match[2] == m1);
   g_assert(ml->size == 4);
   g_assert(ml->pos  == 3);

   free(m0);
   free(m1);
   free(ml);
   
}

// mergesort_mt compar functions.
//int           compar_seqsort   (const void * a, const void * b, const int val);
//int           compar_matchid   (const void * a, const void * b, const int param);
//int           compar_readstart(const void * a, const void * b, const int param);
//int           compar_readend  (const void * a, const void * b, const int param);
//int           compar_matchsize (const void * a, const void * b, const int param);
//int           compar_refstart  (const void * a, const void * b, const int param);

void
test_compar_seqsort
(void)
{
   // Alphabetical comparison of the first K nucleotides of two subsequences.
   // If the K nucleotides are equal, sorting is based on seqid, lower first.
   sub_t * sub1 = malloc(sizeof(sub_t));
   sub_t * sub2 = malloc(sizeof(sub_t));
   
   // sub1.seq < sub2.seq. (22 nucleotides compared)
   sub1->seq = "ACACATCATCCAGACGTACGCGTATATTACGATC";
   sub2->seq = "ATCATCCGACTACGCAGTCGATCACACGCAGTCA";
   g_assert_cmpint(compar_seqsort((void *)sub1, (void *)sub2, 22), ==, -1);

   // sub1.seq > sub2.seq.
   sub2->seq = "ACACATCATCCACACAGCATCACGACGCTACGAC";
   g_assert_cmpint(compar_seqsort((void *)sub1, (void *)sub2, 22), ==, 1);
   
   // sub1.seq == sub2.seq. sub1.seqid > sub2.seqid.
   sub2->seq = "ACACATCATCCAGACGTACGCGAAAAAAAAAAAA";
   sub1->seqid = 321414;
   sub2->seqid = 129133;
   g_assert_cmpint(compar_seqsort((void *)sub1, (void *)sub2, 22), ==, 1);

   // sub1.seq == sub2.seq. sub1.seqid < sub2.seqid.
   sub1->seqid = 122333;
   g_assert_cmpint(compar_seqsort((void *)sub1, (void *)sub2, 22), ==, -1);

   // sub1.seq == sub2.seq. sub1.seqid < sub2.seqid.
   sub1->seqid = 122333;
   g_assert_cmpint(compar_seqsort((void *)sub1, (void *)sub2, 22), ==, -1);

   // sub1.seq == sub2.seq. sub1.seqid == sub2.seqid.
   sub1->seqid = sub2->seqid;
   g_assert_cmpint(compar_seqsort((void *)sub1, (void *)sub2, 22), ==, -1);
   
   free(sub1);
   free(sub2);
}

void
test_compar_matchid
(void)
{
   // This funcrion compares a pair of matches based on their identity:
   // if identity(arg1) < identity(arg2) return +1.
   // else return -1.

   match_t * m1 = malloc(sizeof(match_t));
   match_t * m2 = malloc(sizeof(match_t));
   
   // m1.ident > m2.ident.
   m1->ident = 0.854362;
   m2->ident = 0.845413;
   g_assert_cmpint(compar_matchid((void *)(&m1), (void *)(&m2), 0), ==, -1);
   // m1.ident < m2.ident.
   m1->ident = 0.845412;
   g_assert_cmpint(compar_matchid((void *)(&m1), (void *)(&m2), 0), ==, 1);
   
   // m1.ident == m2.ident. [Warning, using == with double]
   m2->ident = m1->ident;
   g_assert_cmpint(compar_matchid((void *)(&m1), (void *)(&m2),0), ==, -1);
   
   free(m1);
   free(m2);
}


void
test_compar_readstart
(void)
{
   // This funcrion compares a pair of matches based on their read start position:
   // if read_start(arg1) < read_start(arg2) return -1.
   // if read_start(arg1) > read_start(arg2) return +1.
   // if read_start(arg1) = read_start(arg2) then:
   //   if read_end(arg1) > read_end(arg2) return +1.
   //   else return +1.

   match_t * m1 = malloc(sizeof(match_t));
   match_t * m2 = malloc(sizeof(match_t));
   
   // m1.start < m2.start.
   m1->read_s = 10324;
   m2->read_s = 33255;
   g_assert_cmpint(compar_readstart((void *)(&m1), (void *)(&m2), 0), ==, -1);
   // m1.start > m2.start.
   m1->read_s = 44813;
   g_assert_cmpint(compar_readstart((void *)(&m1), (void *)(&m2), 0), ==, 1);
   
   // m1.start == m2.start, m1.end < m2.end.
   m1->read_s = 33255;
   m1->read_e = 33343;
   m2->read_e = 34345;
   g_assert_cmpint(compar_readstart((void *)(&m1), (void *)(&m2),0), ==, -1);
   
   // m1.start == m2.start, m1.end > m2.end.
   m1->read_e = 44235;
   g_assert_cmpint(compar_readstart((void *)(&m1), (void *)(&m2),0), ==, 1);
   
   // m1.start == m2.start, m1.end == m2.end.
   m1->read_e = 34345;
   g_assert_cmpint(compar_readstart((void *)(&m1), (void *)(&m2),0), ==, -1);

   free(m1);
   free(m2);
}

void
test_compar_readend
(void)
{
   // This funcrion compares a pair of matches based on their read end position:
   // if read_end(arg1) < read_end(arg2) return -1.
   // if read_end(arg1) > read_end(arg2) return +1.
   // if read_end(arg1) = read_end(arg2) then:
   //   if read_start(arg1) > read_start(arg2) return +1.
   //   else return +1.

   match_t * m1 = malloc(sizeof(match_t));
   match_t * m2 = malloc(sizeof(match_t));
   
   // m1.end < m2.end.
   m1->read_s = 1000;
   m1->read_e = 1500;
   m2->read_s = 1505;
   m2->read_e = 1650;
   g_assert_cmpint(compar_readend((void *)(&m1), (void *)(&m2), 0), ==, -1);

   // m1.end > m2.end.
   m1->read_e = 1800;
   g_assert_cmpint(compar_readend((void *)(&m1), (void *)(&m2), 0), ==, 1);
   
   // m1.end == m2.end, m1.start < m2.start.
   m1->read_e = 1650;
   g_assert_cmpint(compar_readend((void *)(&m1), (void *)(&m2),0), ==, -1);
   
   // m1.end == m2.end, m1.stat > m2.start.
   m1->read_s = 1550;
   g_assert_cmpint(compar_readend((void *)(&m1), (void *)(&m2),0), ==, 1);
   
   // m1.end == m2.end, m1.start == m2.start.
   m1->read_s = 1505;
   g_assert_cmpint(compar_readend((void *)(&m1), (void *)(&m2),0), ==, -1);

   free(m1);
   free(m2);
   
}

void
test_compar_matchspan
(void)
{
   // This function compares a pair of matches based on their read span.
   // if span(arg1) > span(arg2) return -1.
   // if span(arg1) < span(arg2) return +1.
   // if span(arg1) = span(arg2) then:
   //   if start(arg1) > start(arg2) return +1.
   //   else return -1.

   match_t * m1 = malloc(sizeof(match_t));
   match_t * m2 = malloc(sizeof(match_t));
   
   // m1.span > m2.span.
   m1->read_s = 100;
   m1->read_e = 350; // span(arg1) = 250.
   m2->read_s = 400;
   m2->read_e = 500; // span(arg2) = 100.
   g_assert_cmpint(compar_matchspan((void *)(&m1), (void *)(&m2), 0), ==, -1);
   // m1.span < m2.span.
   m2->read_e = 800; // span(arg2) = 400;
   g_assert_cmpint(compar_matchspan((void *)(&m1), (void *)(&m2), 0), ==, 1);
   
   // m1.span == m2.span, m1.start < m2.start.
   m1->read_e = 500; // span(arg1) = span(arg2) = 400.
   g_assert_cmpint(compar_matchspan((void *)(&m1), (void *)(&m2),0), ==, -1);
   
   // m1.span == m2.span, m1.start > m2.start.
   m1->read_s += 1000;
   m1->read_e += 1000;
   g_assert_cmpint(compar_matchspan((void *)(&m1), (void *)(&m2),0), ==, 1);
   
   // m1.span == m2.span, m1.start == m2.start.
   m2->read_s += 700;
   m2->read_e += 700;
   g_assert_cmpint(compar_matchspan((void *)(&m1), (void *)(&m2),0), ==, -1);

   free(m1);
   free(m2);
   
}

void
test_compar_refstart
(void)
{
   // This funcrion compares a pair of matches based on their reference (genome) start position:
   // if ref_start(arg1) < ref_start(arg2) return -1.
   // if ref_start(arg1) > ref_start(arg2) return +1.
   // if ref_start(arg1) = ref_start(arg2) then:
   //   if ref_end(arg1) > ref_end(arg2) return +1.
   //   else return +1.

   match_t * m1 = malloc(sizeof(match_t));
   match_t * m2 = malloc(sizeof(match_t));
   
   // m1.start < m2.start.
   m1->ref_s = 10324;
   m2->ref_s = 33255234;
   g_assert_cmpint(compar_refstart((void *)(&m1), (void *)(&m2), 0), ==, -1);
   // m1.start > m2.start.
   m1->ref_s = 248134323;
   g_assert_cmpint(compar_refstart((void *)(&m1), (void *)(&m2), 0), ==, 1);
   
   // m1.start == m2.start, m1.end < m2.end.
   m1->ref_s = 33255234;
   m1->ref_e = 33343341;
   m2->ref_e = 34345453;
   g_assert_cmpint(compar_refstart((void *)(&m1), (void *)(&m2),0), ==, -1);
   
   // m1.start == m2.start, m1.end > m2.end.
   m1->ref_e = 4423432434;
   g_assert_cmpint(compar_refstart((void *)(&m1), (void *)(&m2),0), ==, 1);
   
   // m1.start == m2.start, m1.end == m2.end.
   m1->ref_e = 34345453;
   g_assert_cmpint(compar_refstart((void *)(&m1), (void *)(&m2),0), ==, -1);

   free(m1);
   free(m2);
}


//

int
main
(
   int argc,
   char **argv
)
{
   // Save the stderr file descriptor upon start.
   BACKUP_FILE_DESCRIPTOR = dup(STDERR_FILENO);
   BACKUP_STDOUT = dup(STDOUT_FILENO);

   g_test_init(&argc, &argv, NULL);
   g_test_add_func("/hitmap/fuse_matches", test_fuse_matches);
   g_test_add_func("/hitmap/find_repeats", test_find_repeats);
   g_test_add_func("/hitmap/combine_matches", test_combine_matches);
   g_test_add_func("/hitmap/fill_gaps", test_fill_gaps);
   g_test_add_func("/hitmap/test_feedback_gaps", test_feedback_gaps);
   g_test_add_func("/hitmap/matchlist_new", test_matchlist_new);
   g_test_add_func("/hitmap/matchlist_add", test_matchlist_add);
   g_test_add_func("/hitmap/test_compar_seqsort", test_compar_seqsort);
   g_test_add_func("/hitmap/test_compar_matchid", test_compar_matchid);
   g_test_add_func("/hitmap/test_compar_readstart", test_compar_readstart);
   g_test_add_func("/hitmap/test_compar_readend", test_compar_readend);
   g_test_add_func("/hitmap/test_compar_matchspan", test_compar_matchspan);
   g_test_add_func("/hitmap/test_compar_refstart", test_compar_refstart);
   return g_test_run();
}
