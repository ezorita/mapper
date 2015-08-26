SRC_DIR= src
INC_DIR= src
OBJ_DIR= obj
OBJECT_FILES= dc3.o align.o
SOURCE_FILES= mapper-main.c bwmapper.c hitmap.c algs.c poucet.c
HEADER_FILES= bwmapper.h hitmap.h poucet.h definitions.h
ISOURCE_FILES= index.c indexquery.c divsufsort.c
IHEADER_FILES= index.h indexquery.h divsufsort.h divsufsort_private.h algs.h
GSOURCE_FILES= gbrowser.c algs.c
GHEADER_FILES= bwmapper.h algs.h
ASOURCE_FILES= aligner.c align.c
AHEADER_FILES= align.h

OBJECTS= $(addprefix $(OBJ_DIR)/,$(OBJECT_FILES))
SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
HEADERS= $(addprefix $(INC_DIR)/,$(HEADER_FILES))
GSOURCES= $(addprefix $(SRC_DIR)/,$(GSOURCE_FILES))
GHEADERS= $(addprefix $(INC_DIR)/,$(GHEADER_FILES))
ISOURCES= $(addprefix $(SRC_DIR)/,$(ISOURCE_FILES))
IHEADERS= $(addprefix $(INC_DIR)/,$(IHEADER_FILES))
ASOURCES= $(addprefix $(SRC_DIR)/,$(ASOURCE_FILES))
AHEADERS= $(addprefix $(INC_DIR)/,$(AHEADER_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))

CFLAGS= -std=c99 -Wall -g #-O3 -march=core-avx-i -mpopcnt
LDLIBS=
ILDLIBS= -lpthread
CC= gcc

all: bwmapper

bwmapper: $(OBJECTS) $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

buildindex: $(ISOURCES) $(IHEADERS)
	$(CC) $(CFLAGS) $(ISOURCES) $(ILDLIBS) -o $@

gbrowser: $(GSOURCES) $(GHEADERS)
	$(CC) $(CFLAGS) $(GSOURCES) $(LDLIBS) -o $@

aligner: $(ASOURCES) $(AHEADERS)
	$(CC) $(CFLAGS) $(ASOURCES) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
clean:
	rm -f $(OBJECTS) bwmapper
