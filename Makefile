CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors
LDFLAGS = -lm

symnmf: symnmf.o
	gcc $(CFLAGS) symnmf.o -o symnmf $(LDFLAGS)

symnmf.o: symnmf.c symnmf.h
	gcc $(CFLAGS) -c symnmf.c

clean:
	rm -f symnmf.o symnmf wrap_alloc.o symnmf_wrap