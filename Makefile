CC=g++
#CFLAGS=-c -ffast-math -O3 -std=c++0x -Wall -pg -lm
#CFLAGS=-c -std=c++0x -DDOUBLE -O3 -Wall -pg -lm
#LDFLAGS= -pg
CFLAGS=-c -std=c++0x -DDOUBLE -O3 -Wall -lm
LDFLAGS= 
SOURCES=main.cpp QuadNode.cpp QuadTree.cpp Solver.cpp Fields.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=gala2d

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE)

