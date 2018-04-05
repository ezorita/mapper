vpath %.c lib
vpath %.h lib

P= runtests
SOURCE_DIR= ../src
HEADER_DIR= $(SOURCE_DIR)

OBJECTS= test_indexbuild.o libunittest.so
SOURCE= indexbuild.c index.c divsufsort.c bwtquery.c annotate.c blocksearch.c
HEADER= indexbuild.h index.h divsufsort.h bwtquery.h annotate.h blocksearch.h definitions.h
INCLUDE= lib $(SOURCE_DIR) $(HEADER_DIR)

CC= gcc
INCLUDES= $(addprefix -I, $(INCLUDE))
COVERAGE= -fprofile-arcs -ftest-coverage
SOURCES= $(addprefix $(SOURCE_DIR)/,$(SOURCE))
HEADERS= $(addprefix $(HEADER_DIR)/,$(HEADER))
CFLAGS= -std=gnu99 -g -Wall -O0 $(INCLUDES) $(COVERAGE)
LDLIBS= -L`pwd` -Wl,-rpath=`pwd` -lunittest -lpthread

$(P): $(OBJECTS) $(SOURCES) $(HEADERS) runtests.c
	$(CC) $(CFLAGS) runtests.c $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

clean:
	rm -f $(P) $(OBJECTS) *.gcda *.gcno *.gcov .inspect.gdb

libunittest.so: unittest.c
	$(CC) -fPIC -shared $(CFLAGS) -o libunittest.so lib/unittest.c

test: $(P)
	./$(P)

inspect: $(P)
	gdb --command=.inspect.gdb --args $(P)

valgrind: $(P)
	#valgrind --vgdb=yes --vgdb-error=0 ./$(P)
	valgrind --leak-check=full ./$(P)

vgdb: $(P)
	valgrind --vgdb=yes --vgdb-error=0 ./$(P)
