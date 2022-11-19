CC=g++
#CC=clang
TARGET = netAlign
CFLAGS=-fopenmp -O3 -Wall

LDFLAGS=$(CFLAGS)
INCLUDES_MATCH=$(MATCHBOX_DIR)/include
LIBSMATCH=$(MATCHBOX_DIR)/lib/optimized/libmatchbox.a
OBJECTS = \
  matchingHalfApproxDominatingNew.o \
  parseInputFiles.o  \
  sparseImpl.o \
  netAlignExactMatch.o \
  netAlignImplMP.o \
  netAlignImplMR.o \
  netAlignImpl.o $(TARGET).o

all: $(TARGET) message

# fix for a quick hack to avoid lots of 
# code duplication
netAlignImplMP.o : netAlignImplMPTasked.cpp

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) -o $@ $(OBJECTS) 
.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(TARGET) $(OBJECTS)

message:
	echo "Executable: $(TARGET) has been created"
