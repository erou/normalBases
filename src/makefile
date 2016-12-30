CC=gcc
CFLAGS=-I/usr/include/flint -lflint -lmpfr -lgmp -lpthread
DEPS = main.h
OBJ = main.o random.o is_normal.o
TEST = random.o is_normal.o test.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)
test: $(TEST)
	gcc -o $@ $^ $(CFLAGS)
