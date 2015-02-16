SRC_DIR= src
INC_DIR= src
OBJ_DIR= obj
OBJECT_FILES= dc3.o align.o
SOURCE_FILES= mapper-main.c bwmapper.c hitmap.c algs.c poucet.c
HEADER_FILES= bwmapper.h hitmap.h algs.h poucet.h definitions.h
GSOURCE_FILES= gbrowser.c algs.c
GHEADER_FILES= bwmapper.h algs.h

OBJECTS= $(addprefix $(OBJ_DIR)/,$(OBJECT_FILES))
SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
HEADERS= $(addprefix $(INC_DIR)/,$(HEADER_FILES))
GSOURCES= $(addprefix $(SRC_DIR)/,$(GSOURCE_FILES))
GHEADERS= $(addprefix $(INC_DIR)/,$(GHEADER_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))

CFLAGS= -std=c99 -g -Wall -O3
LDLIBS= -lpthread -lm
CC= gcc

all: bwmapper

bwmapper: $(OBJECTS) $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

gbrowser: $(GSOURCES) $(GHEADERS)
	$(CC) $(CFLAGS) $(GSOURCES) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
clean:
	rm -f $(OBJECTS) bwmapper
