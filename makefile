# --- macros -------------------------------------------------------------

SCC	= scc
SCCOPT	= -vv -w -g -d -ww

CC=gcc
CFLAGS=-lm

# --- SpecC rules --------------------------------------------------------

.SUFFIXES:
.SUFFIXES:	.sc .c


make: canny.c canny.sc
	$(CC) canny.c $(CFLAGS) -o canny_gcc_o
	$(SCC) canny -sc2out $(SCCOPT)
   
clean:
	rm canny canny_gcc_o

test:
	./canny beachbus.pgm; 	
	
