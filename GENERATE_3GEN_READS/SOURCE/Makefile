CSRC= generate_3gen_reads.c arguments.c util.c string_util.c

COBJ= $(CSRC:.c=.o)

OPTIM= -g -DKLUDGE
DEBUG= -Wall

CC = gcc

CFLAGS= $(OPTIM) $(DEBUG) 
LFLAGS = -lm

generate_3gen_reads: $(COBJ)
	$(CC) -o $@ $(COBJ) $(LFLAGS) 
	@echo
	@echo    "******************************************************"
	@echo    "*         program compilation successfull            *"
	@echo    "******************************************************"
	@echo

clean:
	rm -f *.o  generate_3gen_reads

generate_3gen_reads.o: alloc.h
arguments.o: alloc.h
util.o: alloc.h
util.o: alloc.h
string_util.o: alloc.h
