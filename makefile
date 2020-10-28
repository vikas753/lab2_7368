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
	rm canny canny_gcc_o beachbus.pgm_* sldl_output*

test:
	./canny beachbus.pgm;
	./canny_gcc_o beachbus.pgm 0.6 0.3 0.8;
	cmp -l beachbus.pgm_s_0.60_l_0.30_h_0.80.pgm sldl_output_0_sc.pgm;
	@true;
 
     	
	
