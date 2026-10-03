CC = gcc
CFLAGS = -std=c99 -Wall -Wextra -Wpedantic
LDLIBS = -lm

matrix_calculator: sparse_matrix_calculator.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

matrix_regressions: tests/test_regressions.c
	$(CC) $(CFLAGS) -I. -o $@ $< $(LDLIBS)

asan: tests/test_regressions.c
	$(CC) $(CFLAGS) -g -fsanitize=address,undefined -I. -o matrix_regressions_asan $< $(LDLIBS)

clean:
	rm -f matrix_calculator matrix_regressions matrix_regressions_asan *.o
