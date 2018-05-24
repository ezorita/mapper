SRC_DIR= src
INC_DIR= src
OBJ_DIR= obj
OBJECT_FILES= index_bwt.o index_sar.o index_txt.o index_sym.o index_ann.o index.o blocksearch.o errhandler.o user_interface.o mapper.o divsufsort.o
SOURCE_FILES= main.c
HEADER_FILES= version.h
OTHER_FILES = Makefile

OBJECTS= $(addprefix $(OBJ_DIR)/,$(OBJECT_FILES))
SOURCES= $(addprefix $(SRC_DIR)/,$(SOURCE_FILES))
HEADERS= $(addprefix $(INC_DIR)/,$(HEADER_FILES))
INCLUDES= $(addprefix -I, $(INC_DIR))

#CFLAGS= -std=c99 -Wall -O3 -mpopcnt
#CFLAGS= -std=c99 -Wall -g -pg -O0 -mpopcnt
CFLAGS= -std=c99 -Wall -g -O0 -mpopcnt
LDLIBS= -lm -lpthread
CC= gcc

all: dirs mapper

dirs: $(OBJ_DIR)

mapper: $(OBJECTS) $(SOURCES) $(HEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

clean:
	rm -f $(OBJECTS) mapper
