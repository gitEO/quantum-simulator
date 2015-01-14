# use "gcc" to compile source files.
CC = g++
# the linker is also "gcc". It might be something else with other compilers.
LD = g++
# Compiler flags go here.
CFLAGS = -Wno-deprecated -O3 -I/site/Blitz++-0.8/include
# Linker flags go here. Currently there aren't any, but if we'll switch to
# code optimization, we might add "-s" here to strip debug info and symbols.
LDFLAGS = -L/site/Blitz++-0.8/lib -lblitz -lm
# use this command to erase files.
RM = /bin/rm -f
# list of generated object files.
OBJS = main.o ny.o
# program executable file name.
PROG = main

# top-level rule, to compile everything.
all: $(PROG)

# rule to link the program
$(PROG): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(PROG)

# rule for file "main.o".
main.o: main.cpp main.h ny.h 
	$(CC) $(CFLAGS) -c main.cpp

# rule for file "file1.o".
ny.o: ny.cpp main.h ny.h
	$(CC) $(CFLAGS) -c ny.cpp


# rule for cleaning re-compilable files.
clean:
	$(RM) $(PROG) $(OBJS)
