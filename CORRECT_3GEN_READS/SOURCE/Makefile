CSRC= correct_3gen_reads.c arguments.c hash_table.c read_sequences.c sequences_data.c string_util.c util.c

COBJ= $(CSRC:.c=.o)

OPTIM= -g
DEBUG= -Wall

CC = gcc

CFLAGS= $(OPTIM) $(DEBUG) 
LFLAGS = -lm

correct_3gen_reads: $(COBJ)
	$(CC) -o $@ $(COBJ) $(LFLAGS) 
	@echo
	@echo    "******************************************************"
	@echo    "*         program compilation successfull            *"
	@echo    "******************************************************"
	@echo

clean:
	rm -f *.o  correct_3gen_reads

correct_3gen_reads.o: alloc.h
argument.o: alloc.h
hash_table.o: alloc.h struct.h
read_sequences.o: alloc.h
sequences_data.o: alloc.h
string_util.o: alloc.h
util.o: alloc.h
