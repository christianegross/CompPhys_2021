# Compiler settings:
CC := gcc
CFLAGS = -c
DEBUG :=
WARNING := 
OPTIMIZATION := -O2
CFLAGS += $(PROFILING) $(DEBUG) $(WARNING) $(OPTIMIZATION)

# Archiver settings:
AR := ar
ARFLAGS = -rcs

# Directories:
OBJDIR := obj
SRCDIR := src
DISTDIR := dist

LIBFILE := libgmm.a
HEADERFILE := generalmemorymanager.h

cfiles := $(wildcard $(SRCDIR)/*.c)

objects := $(subst $(SRCDIR)/,$(OBJDIR)/,$(cfiles:.c=.o))
deps := $(objects:.o=.d)

.PHONY:	all clean dist cleanall

all: $(DISTDIR)/$(LIBFILE)
-include $(deps)

# generate dependencies
$(OBJDIR)/%.d : $(SRCDIR)/%.c
	mkdir -p $(@D)
	@echo [make] finding headers of $< ...
	$(CC) -MM -MT "$@ $(patsubst %.d,%.o,$@)" -MF $@ $<

# compiling
$(OBJDIR)/%.o : $(SRCDIR)/%.c
	@echo [make] compiling $< ...
	$(CC) $(CFLAGS) $< -o $@

# packaging
$(DISTDIR)/$(LIBFILE) : $(objects) Makefile
	@echo [make] packaging $@ ...
	mkdir -p $(DISTDIR)
	$(AR) $(ARFLAGS) $@ $(objects)
	cp $(SRCDIR)/$(HEADERFILE) $(DISTDIR)

#clean directrories
clean:
	$(RM) -r -f $(OBJDIR)
	$(RM) -r -f $(DISTDIR)
