CC=gcc
CFLAGS=-I/usr/include/flint -lflint -lmpfr -lgmp -lpthread
DEPS = main.h
OBJ = main.o random.o is_normal.o sigma_order.o factor_refinement.o luneburg.o frobenius_composition.o
TEST = is_normal.o test.o random.o sigma_order.o factor_refinement.o luneburg.o frobenius_composition.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	gcc -o $@ $^ $(CFLAGS)
test: $(TEST)
	gcc -o $@ $^ $(CFLAGS)
