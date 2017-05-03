# name of the main program (file with int main() function
MAIN_PROGRAM=assignment

# header files for include
IDIR=include
# compiler
CC=mpicc
# compiler flags
CFLAGS=-I $(IDIR) # include header files for compilation

# directory with source code
SDIR=src
# directory to store object files
ODIR=obj

# create directory to store compiled object files
dir_guard=@mkdir -p $(@D)

# generate paths of included files
_DEPS = help.h util.h task1.h task2a.h task2b.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# generate paths of object files
_OBJ = help.o util.o task1.o task2a.o task2b.o $(MAIN_PROGRAM).o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# compile dependencies
$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(dir_guard)
	$(CC) -c -o $@ $< $(CFLAGS)

# compile main files with dependencies
$(MAIN_PROGRAM): $(OBJ)
	$(dir_guard)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

# clean task to clean compiled files
clean:
	find . -type f | xargs touch # to get rid of timestamp warnings due to remote usage
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~
	rm -f $(MAIN_PROGRAM)