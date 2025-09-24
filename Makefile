CFLAGS = -ansi -Wall -Wextra -Werror -pedantic-errors -g -O0 -fno-omit-frame-pointer -ggdb3
LDFLAGS = -lm
GFLAGS = 

symnmf: symnmf.o
	gcc $(CFLAGS) symnmf.o -o symnmf $(LDFLAGS)

symnmf.o: symnmf.c symnmf.h
	gcc $(CFLAGS) -c symnmf.c

clean:
	rm -f symnmf.o symnmf wrap_alloc.o symnmf_wrap

symnmf_wrap: symnmf.o wrap_alloc.o
	gcc $(CFLAGS) symnmf.o wrap_alloc.o -o symnmf_wrap $(LDFLAGS) \
	    -Wl,--wrap=malloc -Wl,--wrap=calloc -Wl,--wrap=realloc

wrap_alloc.o: wrap_alloc.c
	gcc $(CFLAGS) -c wrap_alloc.c

.PHONY: clean