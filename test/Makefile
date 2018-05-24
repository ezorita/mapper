vpath %.c lib
vpath %.h lib

P= runtests
SOURCE_DIR= ../src
HEADER_DIR= $(SOURCE_DIR)

OBJECTS= test_index_sym.o test_index_txt.o test_index_sar.o test_index_bwt.o test_index_ann.o test_index.o test_user_interface.o test_mem_index_sym.o test_mem_index_txt.o libunittest.so
SOURCE= index_sym.c index_txt.c index_sar.c index_bwt.c index_ann.c index.c blocksearch.c errhandler.c user_interface.c mapper.c divsufsort.c
HEADER= index_sym.h index_txt.h index_sar.h index_bwt.h index_ann.h index.h blocksearch.h errhandler.h user_interface.h mapper.h divsufsort.h divsufsort_private.h
INCLUDE= lib $(SOURCE_DIR) $(HEADER_DIR)
COVFILE= blocksearch.gcno index_sym.gcno index_txt.gcno index_bwt.gcno index_sar.gcno index_ann.gcno index.gcno

CC= gcc
INCLUDES= $(addprefix -I, $(INCLUDE))
COVERAGE= -fprofile-arcs -ftest-coverage
SOURCES= $(addprefix $(SOURCE_DIR)/,$(SOURCE))
HEADERS= $(addprefix $(HEADER_DIR)/,$(HEADER))
CFLAGS= -std=gnu99 -g -Wall -O0 $(INCLUDES) $(COVERAGE)
LDLIBS= -L`pwd` -Wl,-rpath=`pwd` -lunittest -lpthread

$(P): $(OBJECTS) $(SOURCES) $(HEADERS) runtests.c
	$(CC) $(CFLAGS) runtests.c $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

clean: clean-examples
	rm -f $(P) $(OBJECTS) *.gcda *.gcno *.gcov .inspect.gdb

clean-examples:
	rm -f *.sym *.txt *.bwt *.sym *.txt *.sar *.ann *.ann.*

libunittest.so: unittest.c
	$(CC) -fPIC -shared $(CFLAGS) -o libunittest.so lib/unittest.c

test: clean-examples $(P)
	./$(P)

$(COVFILE): test


coverage: $(COVFILE)
	gcov $(COVFILE)

inspect: $(P)
	gdb --command=.inspect.gdb --args $(P)

valgrind: $(P)
	#valgrind --vgdb=yes --vgdb-error=0 ./$(P)
	valgrind --leak-check=full ./$(P)

vgdb: $(P)
	valgrind --vgdb=yes --vgdb-error=0 ./$(P)
