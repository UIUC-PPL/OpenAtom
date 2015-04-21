include config.mk

C_FILES=$(wildcard *.C)
CI_FILES=$(wildcard *.ci)
MODULES=$(CI_FILES:.ci=.decl.h)
OBJS=$(addprefix $(BUILDDIR)/, $(C_FILES:.C=.o))

default: all
all: $(BUILDDIR) $(BUILDDIR)/gw_bse

.PHONY : $(DIRS)

$(BUILDDIR):
	mkdir $(BUILDDIR)

$(BUILDDIR)/gw_bse: make_ci $(OBJS) $(DIRS)
	$(CHARMC) $(LFLAGS) $(OBJS) $(BUILDDIR)/*.o -o $@
	mv charmrun $(BUILDDIR)

make_ci: $(MODULES)
	for dir in $(DIRS); do make -C $${dir} make_ci; done

$(DIRS):
	make -C $@

%.decl.h: %.ci
	$(CHARMC) $<

$(BUILDDIR)/%.o: %.C %.h $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

$(BUILDDIR)/%.o: %.C $(MODULES)
	$(CHARMC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) *.decl.h *.def.h
	for dir in $(DIRS); do make -C $${dir} clean; done
