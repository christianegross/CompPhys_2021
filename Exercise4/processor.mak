# General settings:
PROFILING :=

# Compiler settings:
CC := gcc
CFLAGS = -c
DEBUG :=
WARNING := -Wall -Wpedantic
OPTIMIZATION := -Og
CFLAGS += $(PROFILING) $(DEBUG) $(WARNING) $(OPTIMIZATION)

# Directories:
OBJDIR := obj
BINDIR := bin
SRCDIR := src
DISTDIR := dist
DATADIR := data
GPDIR := gplot
LIBDIR := libs

# Linker settings:
LDFLAGS = -L libs
LIBS := -lm -lgsl -lgslcblas -fopenmp
LDFLAGS += $(PROFILING) $(LIBS)



EXECUTIONFILE := ising1d

objects := $(subst $(SRCDIR)/,$(OBJDIR)/,$(cfiles:.c=.o))
deps := $(objects:.o=.d)

.PHONY:	all clean dist cleanall

all: $(BINDIR)/$(EXECUTIONFILE)
-include $(deps)

# List all make variables
list-variables:
	@echo cfiles=$(cfiles)
	@echo objects=$(objects)
	@echo deps=$(deps)
	@echo CFLAGS=$(CFLAGS)
	@echo LDFLAGS=$(LDFLAGS)

# generate dependencies
$(OBJDIR)/%.d : $(SRCDIR)/%.c
	mkdir -p $(@D)
	@echo [make] finding headers of $< ...
	$(CC) -MM -MT "$@ $(patsubst %.d,%.o,$@)" -MF $@ $<

# compiling
$(OBJDIR)/%.o : $(SRCDIR)/%.c
	@echo [make] compiling $< ...
	$(CC) $(CFLAGS) $< -o $@

# linking
$(BINDIR)/$(EXECUTIONFILE) : $(objects) processor.mak Makefile
	@echo [make] linking $@ ...
	mkdir -p $(BINDIR)
	mkdir -p $(DATADIR)
	$(CC) -o $@ $(objects) $(LDFLAGS)

#distribute
$(DISTDIR)/$(EXECUTIONFILE).zip:
	mkdir -p $(DISTDIR)
	@echo [make] zipping $< ...
	zip -r $(DISTDIR)/$(EXECUTIONFILE).zip $(SRCDIR) $(GPDIR) $(LIBDIR) processor.mak Makefile CPHomework6.pdf

dist: $(DISTDIR)/$(EXECUTIONFILE).zip


# Builder will call this to install the application before running.
install:
	echo "Installing is not supported"

#plot
plot:
	mkdir -p $(DATADIR)
	@echo [make] plotting ...
	gnuplot $(GPDIR)/script.plt

#main calculation
calc:$(BINDIR)/$(EXECUTIONFILE)
	@echo [make] calculating ...
	$(BINDIR)/$(EXECUTIONFILE)

# Builder uses this target to run the application.
run:
	make calc

#clean directrories
clean:
	$(RM) -r -f $(OBJDIR)
	$(RM) -r -f $(BINDIR)

cleanall:
	$(RM) -r -f $(OBJDIR)
	$(RM) -r -f $(BINDIR)
	$(RM) -r -f $(DATADIR)



