CSRC= correct_3gen_reads.c arguments.c hash_table.c read_sequences.c sequences_data.c string_util.c util.c \
	align_short_reads.c band_dyn_prog.c

COBJ= $(CSRC:.c=.o)

OPTIM= -g -DDEBUG
DEBUG= -Wall
#DEBUG= -Wall -std=c99  # -std= c99 allows declaration of variable i such as in for(int i = 0; i < 100; i++)

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
align_short_reads.o: alloc.h struct.h
band_dyn_prog.o: alloc.h struct.h
