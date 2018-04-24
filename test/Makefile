vpath %.c lib
vpath %.h lib

P= runtests
SOURCE_DIR= ../rsrc
HEADER_DIR= $(SOURCE_DIR)

OBJECTS= test_index_sym.o test_index_txt.o test_index_sar.o libunittest.so
SOURCE= index_sym.c index_txt.c index_sar.c divsufsort.c
HEADER= index_sym.h index_txt.h index_sar.h divsufsort.h divsufsort_private.h
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
	rm -f $(P) $(OBJECTS) *.gcda *.gcno *.gcov .inspect.gdb *.sym *.txt

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
