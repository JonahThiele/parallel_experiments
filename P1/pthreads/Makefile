CC      = gcc
CFLAGS  =  -g -pthread
SRC     = ./P1_pthread.c ../graph.c ../xxhash.c
TARGET  = color_refinement_pthreads

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $@ $(SRC)

clean:
	rm -f $(TARGET)