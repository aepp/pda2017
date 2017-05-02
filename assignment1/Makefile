MAIN_PROGRAM=assignment

IDIR=include
CC=mpicc
CFLAGS=-I $(IDIR)

SDIR=src
ODIR=obj
TEST=tests
#LDIR =lib
#BDIR=bin

#LIBS=-lm

_DEPS = help.h util.h task1.h task2a.h task2b.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = help.o util.o task1.o task2a.o task2b.o $(MAIN_PROGRAM).o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(MAIN_PROGRAM): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	find . -type f | xargs touch
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~
	rm -f $(MAIN_PROGRAM)