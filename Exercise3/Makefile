SRCDIR := src
# generate a complete directrory list
dirs := $(wildcard $(SRCDIR)/*)
#dirs := $(foreach dir,$(SRCDIR),$(wildcard $(dir)/*))
cfiles := $(foreach dir,$(dirs),$(wildcard $(dir)/*.c))
-include processor.mak