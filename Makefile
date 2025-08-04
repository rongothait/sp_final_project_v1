symnmf: symnmf.o symnmf.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.o -o symnmf -lm

symnmf.o: symnmf.c
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c symnmf.c