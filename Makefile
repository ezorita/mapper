SRC_DIR= src
INC_DIR= src
OBJ_DIR= obj
OBJECT_FILES= align.o
SOURCE_FILES= mapper.c indexquery.c algs.c seed.c filter.c format.c
HEADER_FILES= mapper.h indexquery.h algs.h seed.h filter.h format.h index.h
ISOURCE_FILES= indexbuild.c indexquery.c divsufsort.c
IHEADER_FILES= indexbuild.h indexquery.h divsufsort.h index.h divsufsort_private.h algs.h
GSOURCE_FILES= gbrowser.c algs.c
GHEADER_FILES= bwmapper.h algs.h
ASOURCE_FILES= aligner.c align.c
AHEADER_FILES= align.h
OTHER_FILES = Makefile

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

CFLAGS= -std=c99 -Wall -g -O3 -march=core-avx-i -mpopcnt
#CFLAGS= -std=c99 -Wall -g -pg -O0 -march=core-avx-i -mpopcnt
#CFLAGS= -std=c99 -Wall -g -O0 -march=core-avx-i -mpopcnt
LDLIBS= -lm -lpthread
ILDLIBS= 
CC= gcc

all: mapper buildindex

mapper: $(OBJECTS) $(SOURCES) $(HEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

buildindex: $(ISOURCES) $(IHEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(ISOURCES) $(ILDLIBS) -o $@

gbrowser: $(GSOURCES) $(GHEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(GSOURCES) $(LDLIBS) -o $@

aligner: $(ASOURCES) $(AHEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(ASOURCES) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h $(OTHER_FILES)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
clean:
	rm -f $(OBJECTS) mapper buildindex
