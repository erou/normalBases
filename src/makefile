CC=gcc
CFLAGS=-I/usr/include/flint -lflint -lmpfr -lgmp -lpthread
DEPS = main.h
OBJ = main.o random.o is_normal.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)
test: test.o
	gcc -o $@ $^ $(CFLAGS)
