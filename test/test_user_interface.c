#include <unistd.h>
#include <getopt.h>
#include "unittest.h"
#include "user_interface.h"
#include "version.h"


void
test_ui_parse
(void)
{
   char * buffer = malloc(10000);
   test_assert_critical(buffer != NULL);
   
   char * argv0[] = {"./mapper", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(1, argv0) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_MAP);

   char * argv1[] = {"./mapper", "-h", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(2, argv1) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_MAP);

   char * argv2[] = {"./mapper", "--help", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(2, argv2) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_MAP);

   char * argv3[] = {"./mapper", "-v", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(2, argv3) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr("mapper version: "MAPPER_VERSION);

   char * argv4[] = {"./mapper", "--version", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(2, argv4) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr("mapper version: "MAPPER_VERSION);

   char * argv5[] = {"./mapper", "index", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(2, argv5) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_INDEX);

   char * argv6[] = {"./mapper", "index", "-h", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv6) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_INDEX);

   char * argv7[] = {"./mapper", "index", "--help", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv7) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_INDEX);

   char * argv8[] = {"./mapper", "index", "add", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv8) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_ADD);

   char * argv9[] = {"./mapper", "index", "add", "-h", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv9) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_ADD);

   char * argv10[] = {"./mapper", "index", "add", "--help", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv10) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_ADD);

   char * argv11[] = {"./mapper", "index", "build", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv11) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_BUILD);

   char * argv12[] = {"./mapper", "index", "build", "-h", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv12) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_BUILD);

   char * argv13[] = {"./mapper", "index", "build", "--help", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv13) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_BUILD);

   char * argv14[] = {"./mapper", "index", "view", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv14) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_VIEW);

   char * argv15[] = {"./mapper", "index", "view", "-h", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv15) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_VIEW);

   char * argv16[] = {"./mapper", "index", "view", "--help", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv16) == EXIT_SUCCESS);
   unredirect_stderr();
   test_assert_stderr(USAGE_VIEW);


   char * argv17[] = {"./mapper", "fake-index", "fake-input", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv17) == EXIT_FAILURE);
   unredirect_stderr();

   char * argv18[] = {"./mapper", "index", "add", "-k10", "-d-2", "fake-index", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(6, argv18) == EXIT_FAILURE);
   unredirect_stderr();

   char * argv19[] = {"./mapper", "index", "build", "--output", "ui_test00", "examples/repeats.fa", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(6, argv19) == EXIT_SUCCESS);
   unredirect_stderr();

   // Test index.
   index_t * index = index_read("ui_test00");
   test_assert(strcmp(index->fname_base, "ui_test00") == 0);
   test_assert_critical(index != NULL);
   test_assert(index->sym != NULL);
   test_assert(index->txt != NULL);
   test_assert(index->sar != NULL);
   test_assert(index->bwt != NULL);
   test_assert(index->ann == NULL);
   test_assert(index->ann_cnt == 0);

   index_free(index);

   char * argv20[] = {"./mapper", "index", "add", "-k25", "-d1", "-t1", "ui_test00", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(7, argv20) == EXIT_SUCCESS);
   unredirect_stderr();

   // Test index.
   index = index_read("ui_test00");
   test_assert(strcmp(index->fname_base, "ui_test00") == 0);
   test_assert_critical(index != NULL);
   test_assert(index->sym != NULL);
   test_assert(index->txt != NULL);
   test_assert(index->sar != NULL);
   test_assert(index->bwt != NULL);
   test_assert(index->ann != NULL);
   test_assert(index->ann_cnt == 1);

   // Test index content.
   uint8_t seq_one_suf[25] = {0,3,1,2,0,3,0,3,1,0,2,1,1,0,1,3,0,1,2,0,2,0,1,0,0};
   uint8_t seq_two_pre[25] = {0,0,1,0,2,0,2,1,0,3,1,0,1,2,2,0,1,3,0,3,0,2,1,3,0};

   bwtquery_t * q = bwt_new_query(index->bwt);
   for (int i = 0; i < 25; i ++) {
      test_assert(bwt_query(seq_one_suf[i], BWT_QUERY_SUFFIX, q, q) == 0);
   }
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 25);
   test_assert(strcmp(txt_pos_to_str(sar_get(bwt_start(q), index->sar), index->txt), "one:1:+") == 0);
   free(q);

   q = bwt_new_query(index->bwt);
   for (int i = 0; i < 25; i ++) {
      test_assert(bwt_query(seq_two_pre[i], BWT_QUERY_PREFIX, q, q) == 0);
   }
   test_assert(bwt_size(q) == 1);
   test_assert(bwt_depth(q) == 25);
   test_assert(strcmp(txt_pos_to_str(sar_get(bwt_start(q), index->sar), index->txt), "two:1:+") == 0);
   free(q);

   // Test annotation.
   locinfo_t * lci = ann_query(txt_str_to_pos("six:1:+", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 1);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 10);
   free(lci);

   lci = ann_query(txt_str_to_pos("eight:1:+", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 2);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 8);
   free(lci);

   lci = ann_query(txt_str_to_pos("three:1:-", index->txt), index->ann[0]);
   test_assert_critical(lci != NULL);
   test_assert(lci->dist == 1);
   test_assert(lci->neigh_cnt == 1);
   test_assert(lci->align_cnt == 1);
   test_assert(lci->align_pos[0] == 4);
   free(lci);

   index_free(index);


   // Force too many/too few arguments.
   char * argv21[] = {"./mapper", "a", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(2, argv21) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_INSUF_ARG, USAGE_MAP);
   test_assert_stderr(buffer);

   char * argv22[] = {"./mapper", "a", "b", "c", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv22) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_MAX_ARG, USAGE_MAP);
   test_assert_stderr(buffer);

   char * argv23[] = {"./mapper", "index", "build", "c", "d", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv23) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_MAX_ARG, USAGE_BUILD);
   test_assert_stderr(buffer);
   
   char * argv24[] = {"./mapper", "index", "build", "-o", "d", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv24) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_INSUF_ARG, USAGE_BUILD);
   test_assert_stderr(buffer);

   char * argv25[] = {"./mapper", "index", "add", "c", "d", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv25) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_INSUF_ARG, USAGE_ADD);
   test_assert_stderr(buffer);

   char * argv26[] = {"./mapper", "index", "add", "-k13", "d", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv26) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_INSUF_ARG, USAGE_ADD);
   test_assert_stderr(buffer);

   char * argv27[] = {"./mapper", "index", "add", "-k13", "-d1", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv27) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_INSUF_ARG, USAGE_ADD);
   test_assert_stderr(buffer);

   char * argv28[] = {"./mapper", "index", "add", "-k13", "-d1", "c", "d", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(7, argv28) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_MAX_ARG, USAGE_ADD);
   test_assert_stderr(buffer);

   char * argv29[] = {"./mapper", "index", "view", "c", "d", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv29) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_MAX_ARG, USAGE_VIEW);
   test_assert_stderr(buffer);


   // Incorrect arguments.
   char * argv30[] = {"./mapper", "index", "not-a-command", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv30) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%s", ERROR_COMMAND, USAGE_INDEX);
   test_assert_stderr(buffer);

   char * argv31[] = {"./mapper", "-g", "a", "b", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv31) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%c\n%s", ERROR_INCORRECT_OPT, 'g', USAGE_MAP);
   test_assert_stderr(buffer);

   char * argv32[] = {"./mapper", "index", "add", "-j", "a", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv32) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%c\n%s", ERROR_INCORRECT_OPT, 'j', USAGE_ADD);
   test_assert_stderr(buffer);

   char * argv33[] = {"./mapper", "index", "build", "-j", "a", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(5, argv33) == EXIT_FAILURE);
   unredirect_stderr();
   sprintf(buffer, "%s%c\n%s", ERROR_INCORRECT_OPT, 'j', USAGE_BUILD);
   test_assert_stderr(buffer);

   // Other tests.
   char * argv34[] = {"./mapper", "index", "view", "ui_test00", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv34) == EXIT_SUCCESS);
   unredirect_stderr();

   char * argv35[] = {"./mapper", "index", "view", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(4, argv35) == EXIT_FAILURE);
   unredirect_stderr();


   char * argv36[] = {"./mapper", "ui_test00", "fake-input-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(3, argv36) == EXIT_FAILURE);
   unredirect_stderr();

   // Repeated or incorrect option values.
   char * argv37[] = {"./mapper", "index", "add", "-k10", "-d1", "-t-12", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(7, argv37) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_THREAD_POSITIVE);

   char * argv38[] = {"./mapper", "index", "add", "-k10", "-d1", "-t12", "-t8", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(8, argv38) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_THREAD_REPEAT);

   char * argv39[] = {"./mapper", "index", "add", "-k-10", "-d1", "-t12", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(7, argv39) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_KMER_POSITIVE);

   char * argv40[] = {"./mapper", "index", "add", "-k10", "-d1", "-t12", "-k32", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(8, argv40) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_KMER_REPEAT);

   char * argv41[] = {"./mapper", "index", "add", "-k10", "-d-1", "-t12", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(7, argv41) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_DIST_POSITIVE);

   char * argv42[] = {"./mapper", "index", "add", "-k10", "-d1", "-t12", "-d4", "fake-index-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(8, argv42) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_DIST_REPEAT);

   char * argv43[] = {"./mapper", "index", "build", "-o", "output-f1", "-o", "output-f2", "fake-fasta-file", NULL};
   redirect_stderr();
   optind = 1;
   test_assert(ui_parse(8, argv43) == EXIT_FAILURE);
   unredirect_stderr();
   test_assert_stderr(OPT_OUTPUT_REPEAT);


}

// Define test cases to be run (for export).
const test_case_t test_cases_ui[] = {
   {"ui/ui_parse",         test_ui_parse},
   {NULL, NULL}, // Sentinel. //
};

