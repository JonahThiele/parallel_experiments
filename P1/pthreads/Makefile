CC      = gcc
# we always add debuggers in this household
CFLAGS  = -std=c11 -g -pthread

SRC     = P1_pthreads.c ../graph.c ../xxhash.c
OBJ     = $(SRC:.c=.o)

TARGET  = color_refinement_pthreads

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
