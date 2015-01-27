CC=g++
CFLAGS=-std=c++11 -O3 -c -Wall
LDFLAGS=-std=c++11 -O3
SOURCES=src/data/Data.cc \
	src/mesh/Finite_Element.cc \
	src/manatee/Manatee.cc \
	src/mesh/Mesh.cc \
	src/neutronics/Neutronics.cc \
	src/parser/Parser.cc 
OBJECTS=$(SOURCES:.cc=.o)
EXECUTABLE=manatee

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cc.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm src/*/*.o

