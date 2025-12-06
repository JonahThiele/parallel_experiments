CC      = mpicc
CFLAGS  = -Wall -std=c11

SRC_CPP = P1_MPI.c
SRC_C   = ../graph.c ../xxhash.c

OBJ_CPP = $(SRC_CPP:.cpp=.o)
OBJ_C   = $(SRC_C:.c=.o)

OBJ     = $(OBJ_CPP) $(OBJ_C)

TARGET  = color_refinement_MPI

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

# Rule for C++ files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.o: %.c
	mpicc -Wall -std=c11 -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
