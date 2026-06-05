CXX ?= g++
CXXFLAGS ?= -O3 -std=c++17
PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

HTSLIB_DIR ?=
USE_RPATH ?= 0

INCLUDES :=
LIBDIRS :=
RPATH :=
LIBS := -lhts -lz -lbz2 -llzma -lcurl -lpthread

ifneq ($(strip $(HTSLIB_DIR)),)
INCLUDES += -I$(HTSLIB_DIR)
LIBDIRS += -L$(HTSLIB_DIR)
ifeq ($(USE_RPATH),1)
RPATH += -Wl,-rpath,$(HTSLIB_DIR)
endif
endif

all: allele_counter

allele_counter: allele_counter.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(INCLUDES) $(LIBDIRS) $(RPATH) $(LIBS)

install: allele_counter
	install -d $(DESTDIR)$(BINDIR)
	install -m 0755 allele_counter $(DESTDIR)$(BINDIR)/allele_counter

clean:
	rm -f allele_counter

.PHONY: all install clean
