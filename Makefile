    CC     = gcc -g -std=c11
    CFLAGS =
    LFLAGS = -lm

      PROG = matrixInv
      OBJS = utils.o  matriz.o


$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

clean:
	@rm -f *~ *.bak
	@rm -f *.o
	@rm -f $(PROG)

purge:   clean
	@rm -f core *.out
	


