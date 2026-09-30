# ---------------------------------------------------------------------------
# Makefile for crashe (colrad / CRM code)
#
#   make            release build  -> bin/crashe
#   make debug      debug build    -> bin/crashe_debug
#   make test       release build, then python3 test.py
#   make clean      remove obj/ and bin/
#   make -j         parallel builds are safe
#
# Release and debug objects live in separate directories (obj/release,
# obj/debug), so switching between them never reuses stale objects.
#
# Module dependencies are generated automatically from the `use` statements
# in src/*.f90 (see "Dependencies" below), so adding a module or a `use` line
# needs no Makefile edit. The one convention this relies on is that a module
# lives in a file of the same name: module foo  ->  src/foo.f90
# ---------------------------------------------------------------------------

# --- Compiler and tools ------------------------------------------------------
FC       := gfortran
MKDIR_P  := mkdir -p
RM       := rm -rf

# --- Build type --------------------------------------------------------------
BUILD    ?= release

# --- Paths -------------------------------------------------------------------
SRCDIR   := src
OBJDIR   := obj/$(BUILD)
BINDIR   := bin

# --- Flags -------------------------------------------------------------------
# -ffree-line-length-none: several source lines exceed gfortran's default
# 132-character limit, which is a hard error by default.
BASE_FFLAGS := -fopenmp -ffree-line-length-none -J$(OBJDIR) -I$(OBJDIR)
RELEASE_FLG := -O3 -fbacktrace -fcheck=all -g -Warray-temporaries
DEBUG_FLG   := -Og -g -Wall -Wextra -fcheck=all -fbacktrace \
               -ffpe-trap=invalid,zero,overflow -Warray-temporaries

ifeq ($(BUILD),release)
  MODE_FLAGS := $(RELEASE_FLG)
  SUFFIX     :=
else ifeq ($(BUILD),debug)
  MODE_FLAGS := $(DEBUG_FLG)
  SUFFIX     := _debug
else
  $(error BUILD must be 'release' or 'debug', got '$(BUILD)')
endif

FFLAGS   := $(BASE_FFLAGS) $(MODE_FLAGS)
LDFLAGS  := -llapack -lblas

# --- Files -------------------------------------------------------------------
TARGET   := $(BINDIR)/crashe$(SUFFIX)
SRCS     := $(wildcard $(SRCDIR)/*.f90)
OBJS     := $(SRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)
DEPFILE  := $(OBJDIR)/deps.mk

# --- Rules -------------------------------------------------------------------
.PHONY: all debug test clean

all: $(TARGET)

debug:
	$(MAKE) BUILD=debug

$(TARGET): $(OBJS) | $(BINDIR)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

# Order-only prerequisites (after the |): the directory must exist, but its
# timestamp never triggers a rebuild. Safe under make -j.
$(OBJDIR) $(BINDIR):
	@$(MKDIR_P) $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# --- Dependencies ------------------------------------------------------------
# For every src/foo.f90, emit   obj/<build>/foo.o: obj/<build>/bar.o ...
# for each `use bar` that has a matching src/bar.f90. Intrinsic modules
# (iso_fortran_env, omp_lib, ...) have no source file and are skipped, as is a
# file "using" itself. Module names are case-insensitive, hence the tr.
# Compiling a file writes its .mod alongside its .o, so depending on the
# .o of each used module is enough to get the order (and rebuilds) right.
$(DEPFILE): $(SRCS) | $(OBJDIR)
	@echo "Generating $@"
	@for f in $(SRCS); do \
	  b=$$(basename $$f .f90); \
	  mods=$$(grep -i -E '^[[:space:]]*use[[:space:]]+[a-z0-9_]+' $$f \
	          | tr 'A-Z' 'a-z' \
	          | sed -E 's/^[[:space:]]*use[[:space:]]+([a-z0-9_]+).*/\1/' \
	          | sort -u); \
	  deps=""; \
	  for m in $$mods; do \
	    if [ "$$m" != "$$b" ] && [ -f $(SRCDIR)/$$m.f90 ]; then \
	      deps="$$deps $(OBJDIR)/$$m.o"; \
	    fi; \
	  done; \
	  echo "$(OBJDIR)/$$b.o:$$deps"; \
	done > $@

ifneq ($(MAKECMDGOALS),clean)
-include $(DEPFILE)
endif

# --- Housekeeping ------------------------------------------------------------
clean:
	$(RM) obj $(BINDIR)

test: all
	python3 test.py