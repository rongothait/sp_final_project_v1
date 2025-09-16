symnmf: symnmf.o
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors symnmf.o -o symnmf -lm

symnmf.o: symnmf.c symnmf.h
	gcc -ansi -Wall -Wextra -Werror -pedantic-errors -c symnmf.c