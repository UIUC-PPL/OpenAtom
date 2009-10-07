# Setp the basic compiler variables, flags, implicit rules etc
########### This stuff should be able take care of itself ############

# Required directory structure within the charm installation
CHARMBIN  = $(CHARMBASE)/bin
CHARMINC  = $(CHARMBASE)/include

# Compilers
CC        = $(CHARMBIN)/charmc
CXX       = $(CHARMBIN)/charmc
FC        = $(CHARMBIN)/charmc

# Other tools
LN        = ln -sf
DOXYGEN   = doxygen

# Basic compiler/linker flags
CFLAGS   += -language charm++
CXXFLAGS += -language charm++
FFLAGS   += -f77

####### Canned verbosity settings
v ?= 0

ifeq ($(strip $v),0)
  info-dep = @echo Generating dependencies for $(<F)
  info-ci  = @echo Parsing interface definitions in $(<F)
  info-cpp = @echo Compiling $(<F)
  info-c   = @echo Compiling $(<F)
  info-f   = @echo Compiling $(<F)
  info-ar  = @echo =========== Producing archive $@
  info-ld  = @echo =========== Linking to produce $@
  q = @
else 
ifeq ($(strip $v),1)
  info-dep = @echo Generating dependencies for $(<F)
  info-ci  = @echo Parsing interface definitions in $(<F)
  info-cpp = @echo Compiling $(<F) with options $(CXXFLAGS)
  info-c   = @echo Compiling $(<F) with options $(CFLAGS)
  info-f   = @echo Compiling $(<F) with options $(FFLAGS)
  info-ar  = @echo =========== Producing archive $@ containing $^
  info-ld  = @echo =========== Linking to produce $@ with extra libs: $(LDLIBS)
  q = @
else 
ifeq ($(strip $v),2)
  info-dep = @echo Generating dependencies for $(<F)
  info-ci  = @echo Parsing interface definitions in $(<F)
  info-cpp =
  info-c   =
  info-f   =
  info-ar  = @echo ===================================================
  info-ld  = @echo ===================================================
  q =
else
  $(error Wrong level of verboseness input. Use a value between 0-2)
endif
endif
endif

####### Pattern rules
# Rule to generate dependency information for C++ source files
%.d: %.C
	$(info-dep)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

%.d: %.cpp
	$(info-dep)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

%.d: %.cxx
	$(info-dep)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

# Rule to generate dependency information for C source files
%.d: %.c
	$(info-dep)
	@g++ -MM -MG $(CPPFLAGS) $(INCDIRS:%=-I%) -I$(CHARMINC) $< | perl $(CHARMBIN)/dep.pl $(CHARMINC) $(DEPSTRIPDIRS) \
	| sed 's|$*.o[ :]*|$*.o $@ : |g' > $@

# Rule to generate dependency info for charm++ interface (ci) definition files
# @note: Need to handle pathological cases like multi-line module declarations 
%.di: %.ci
	$(info-dep)
	@grep -oE "(extern[ ]+)?module[ ]+\w+" $< | \
	awk ' function printExternDeps(nExt,module,externs) { if (nExt>0) { printf "%s.decl.h: ",module; for (i=0;i<nExt;i++) printf "%s.decl.h ",externs[i]; printf "\n\n" } }   { if ($$1 ~ /extern/) { externs[nExt++] = $$3 } else { printExternDeps(nExt,modules[cnt-1],externs); nExt=0; modules[cnt++] = $$2 } }    END { printExternDeps(nExt,modules[cnt-1],externs); for (i=0;i<cnt;i++) printf "%s.decl.h %s.def.h ",modules[i],modules[i]; printf ": $<\n\t$$(info-ci)\n\t$q$(CXX) -c $$<\n\n" }' > $@


# Pattern rules copied from the built-in make rules
%.o: %.C
	$(info-cpp)
	$q$(COMPILE.C) $(OUTPUT_OPTION) $<

%.o: %.cpp
	$(info-cpp)
	$q$(COMPILE.cpp) $(OUTPUT_OPTION) $<

%.o: %.cc
	$(info-cpp)
	$q$(COMPILE.cc) $(OUTPUT_OPTION) $<

%.o: %.c
	$(info-c)
	$q$(COMPILE.c) $(OUTPUT_OPTION) $<

%.o: %.f
	$(info-f)
	$q$(COMPILE.f) $(OUTPUT_OPTION) $<

%.o: %.F
	$(info-f)
	$q$(COMPILE.F) $(OUTPUT_OPTION) $<

%.f: %.F
	@echo Preprocessing $<
	$q$(PREPROCESS.F) $(OUTPUT_OPTION) $<

