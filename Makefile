SRC_DIR= src
INC_DIR= src
OBJ_DIR= obj
OBJECT_FILES= align.o
SOURCE_FILES= mapper.c bwtquery.c algs.c seed.c filter.c format.c annotate.c xxhash.c blocksearch.c indexbuild.c divsufsort.c index.c interface.c
HEADER_FILES= mapper.h bwtquery.h algs.h seed.h filter.h format.h index.h annotate.h xxhash.h blocksearch.h indexbuild.h divsufsort.h divsufsort_private.h interface.h definitions.h
GSOURCE_FILES= gbrowser.c algs.c
GHEADER_FILES= bwmapper.h algs.h
ASOURCE_FILES= aligner.c align.c
AHEADER_FILES= align.h
DSOURCE_FILES = ann_dump.c index.c annotate.c bwtquery.c blocksearch.c
DHEADER_FILES = $(HEADER_FILES)
OTHER_FILES = Makefile

OBJECTS= $(addprefix $(OBJ_DIR)/,$(OBJECT_FILES))
SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
HEADERS= $(addprefix $(INC_DIR)/,$(HEADER_FILES))
GSOURCES= $(addprefix $(SRC_DIR)/,$(GSOURCE_FILES))
GHEADERS= $(addprefix $(INC_DIR)/,$(GHEADER_FILES))
ASOURCES= $(addprefix $(SRC_DIR)/,$(ASOURCE_FILES))
AHEADERS= $(addprefix $(INC_DIR)/,$(AHEADER_FILES))
DSOURCES= $(addprefix $(SRC_DIR)/,$(DSOURCE_FILES))
DHEADERS= $(addprefix $(INC_DIR)/,$(DHEADER_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))

#CFLAGS= -std=c99 -Wall -O3 -mpopcnt
#CFLAGS= -std=c99 -Wall -g -pg -O0 -mpopcnt
CFLAGS= -std=c99 -Wall -g -O0 -mpopcnt
LDLIBS= -lm -lpthread
CC= gcc

all: mapper

mapper: $(OBJECTS) $(SOURCES) $(HEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

anndump: $(DSOURCES) $(DHEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(DSOURCES) $(LDLIBS) -o $@

gbrowser: $(GSOURCES) $(GHEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(GSOURCES) $(LDLIBS) -o $@

aligner: $(ASOURCES) $(AHEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(ASOURCES) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h $(OTHER_FILES)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
clean:
	rm -f $(OBJECTS) mapper buildindex
