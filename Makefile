SRC_DIR= rsrc
INC_DIR= rsrc
OBJ_DIR= obj
OBJECT_FILES= index_bwt.o
SOURCE_FILES= 
HEADER_FILES= 
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

all: mapper

mapper: $(OBJECTS) $(SOURCES) $(HEADERS) $(OTHER_FILES)
	$(CC) $(CFLAGS) $(SOURCES) $(OBJECTS) $(LDLIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(SRC_DIR)/%.h
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@
clean:
	rm -f $(OBJECTS) mapper
