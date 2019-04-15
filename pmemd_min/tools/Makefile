dirs := . $(shell cat Filepath)
cpp_dirs := $(dirs)

VPATH := $(foreach dir,$(cpp_dirs),$(wildcard $(dir)))
VPATH := $(subst $(space),:,$(VPATH))

SOURCES := $(shell cat Srcfiles)
BASENAMES := $(basename $(basename $(SOURCES)))
OBJS    := $(addsuffix .o, $(BASENAMES))
INCLUDES    := $(foreach dir,$(cpp_dirs),-I$(dir))

include $(CASEROOT)/Macros
all: $(BINDIR)/pmemd.MPI.$(BLDID)
$(BINDIR)/pmemd.MPI.$(BLDID): $(BINDIR)/pmemd.MPI
	cp $< $@
$(BINDIR)/pmemd.MPI: pmemd.MPI
	cp $< $@
pmemd.MPI: $(OBJS)
	$(FC) $(LDFLAGS) $(OBJS) -o $@ $(LIBS)
.PHONY: clean
# Filepath:
# 	echo $(CASEROOT)/SourceMods > Filepath
# 	echo $(PMEMDSRC) >> Filepath
Srcfiles: Filepath
	$(CASETOOL)/mkSrcfiles > Srcfiles
Depends: Filepath Srcfiles
	$(CASETOOL)/mkDepends Filepath Srcfiles > Depends
%.o: %.F90
	$(FC) $(FFLAGS) $< -c -o $@
%_sw64.o: %_sw64.c
	$(SWCC) -host  -DMPE $(CFLAGS) $< -c -o host.$@
	$(SWCC) -slave -DCPE $(CFLAGS) $< -c -o slave.$@
	sw5ld -r host.$@ slave.$@ -o $@
	rm host.$@ slave.$@
%_cpe.o: %_cpe.S
	$(SWCC) -slave -c $< -o $@
%_mpe.o: %_mpe.S
	$(SWCC) -host -c $< -o $@
%.o: %.c
	$(CC) $(CFLAGS) $< -c -o $@
-include Depends
