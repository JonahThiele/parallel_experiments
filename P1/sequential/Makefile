CC      = gcc
# we always add debugging flags in this household
CFLAGS  = -Wall -std=c11 -g

SRC     = P1.c ../graph.c ../xxhash.c
OBJ     = $(SRC:.c=.o)

TARGET  = color_refinement

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET)
