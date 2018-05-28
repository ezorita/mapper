#include "unittest.h"
#include "user_interface.h"
#include "version.h"

void
test_mem_ui_parse
(void)
{
   char ** argvec[45];
   argvec[0] = (char *[]){"./mapper", NULL};
   argvec[1] = (char *[]){"./mapper", NULL};
   argvec[2] = (char *[]){"./mapper", "-h", NULL};
   argvec[3] = (char *[]){"./mapper", "--help", NULL};
   argvec[4] = (char *[]){"./mapper", "-v", NULL};
   argvec[5] = (char *[]){"./mapper", "--version", NULL};
   argvec[6] = (char *[]){"./mapper", "index", NULL};
   argvec[7] = (char *[]){"./mapper", "index", "-h", NULL};
   argvec[8] = (char *[]){"./mapper", "index", "--help", NULL};
   argvec[9] = (char *[]){"./mapper", "index", "add", NULL};
   argvec[10] = (char *[]){"./mapper", "index", "add", "-h", NULL};
   argvec[11] = (char *[]){"./mapper", "index", "add", "--help", NULL};
   argvec[12] = (char *[]){"./mapper", "index", "build", NULL};
   argvec[13] = (char *[]){"./mapper", "index", "build", "-h", NULL};
   argvec[14] = (char *[]){"./mapper", "index", "build", "--help", NULL};
   argvec[15] = (char *[]){"./mapper", "index", "view", NULL};
   argvec[16] = (char *[]){"./mapper", "index", "view", "-h", NULL};
   argvec[17] = (char *[]){"./mapper", "index", "view", "--help", NULL};
   argvec[18] = (char *[]){"./mapper", "fake-index", "fake-input", NULL};
   argvec[19] = (char *[]){"./mapper", "index", "add", "-k10", "-d-2", "fake-index", NULL};
   argvec[20] = (char *[]){"./mapper", "index", "build", "--output", "ui_test10", "examples/repeats.fa", NULL};
   argvec[21] = (char *[]){"./mapper", "index", "add", "-k25", "-d1", "-t1", "ui_test10", NULL};
   argvec[22] = (char *[]){"./mapper", "a", NULL};
   argvec[23] = (char *[]){"./mapper", "a", "b", "c", NULL};
   argvec[24] = (char *[]){"./mapper", "index", "build", "c", "d", NULL};
   argvec[25] = (char *[]){"./mapper", "index", "build", "-o", "d", NULL};
   argvec[26] = (char *[]){"./mapper", "index", "add", "c", "d", NULL};
   argvec[27] = (char *[]){"./mapper", "index", "add", "-k13", "d", NULL};
   argvec[28] = (char *[]){"./mapper", "index", "add", "-k13", "-d1", NULL};
   argvec[29] = (char *[]){"./mapper", "index", "add", "-k13", "-d1", "c", "d", NULL};
   argvec[30] = (char *[]){"./mapper", "index", "view", "c", "d", NULL};
   argvec[31] = (char *[]){"./mapper", "index", "not-a-command", NULL};
   argvec[32] = (char *[]){"./mapper", "-g", "a", "b", NULL};
   argvec[33] = (char *[]){"./mapper", "index", "add", "-j", "a", NULL};
   argvec[34] = (char *[]){"./mapper", "index", "build", "-j", "a", NULL};
   argvec[35] = (char *[]){"./mapper", "index", "view", "ui_test00", NULL};
   argvec[36] = (char *[]){"./mapper", "index", "view", "fake-index-file", NULL};
   argvec[37] = (char *[]){"./mapper", "ui_test00", "fake-input-file", NULL};
   argvec[38] = (char *[]){"./mapper", "index", "add", "-k10", "-d1", "-t-12", "fake-index-file", NULL};
   argvec[39] = (char *[]){"./mapper", "index", "add", "-k10", "-d1", "-t12", "-t8", "fake-index-file", NULL};
   argvec[40] = (char *[]){"./mapper", "index", "add", "-k-10", "-d1", "-t12", "fake-index-file", NULL};
   argvec[41] = (char *[]){"./mapper", "index", "add", "-k10", "-d1", "-t12", "-k32", "fake-index-file", NULL};
   argvec[42] = (char *[]){"./mapper", "index", "add", "-k10", "-d-1", "-t12", "fake-index-file", NULL};
   argvec[43] = (char *[]){"./mapper", "index", "add", "-k10", "-d1", "-t12", "-d4", "fake-index-file", NULL};
   argvec[44] = (char *[]){"./mapper", "index", "build", "-o", "output-f1", "-o", "output-f2", "fake-fasta-file", NULL};
   
   redirect_stderr();

   // Set alloc failure rate to 0.1.
   set_alloc_failure_rate_to(0.01);
   for (int i = 0; i < 1000; i++) {
      for (int j = 0; j < 45; j++) {
         int argc = 0;
         while(argvec[j][argc] != NULL) argc++;
         optind = 1;
         ui_parse(argc, argvec[j]);
      }
   }
   reset_alloc();

   // Set alloc countdown 0->10.
   for (int i = 0; i <= 1000; i++) {
      set_alloc_failure_countdown_to(i);
      for (int j = 0; j < 45; j++) {
         int argc = 0;
         while(argvec[j][argc] != NULL) argc++;
         optind = 1;
         ui_parse(argc, argvec[j]);
      }
   }
   reset_alloc();
   unredirect_stderr();

}

// Define test cases to be run (for export).
const test_case_t test_mem_ui[] = {
   {"mem/ui/ui_parse",         test_mem_ui_parse},
   {NULL, NULL}, // Sentinel. //
};
