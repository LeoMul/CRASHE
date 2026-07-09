# --- Compiler and Tools ---
FC       := gfortran
MKDIR_P  := mkdir -p
RM       := rm -rf

# --- Paths ---
SRCDIR   := src
OBJDIR   := obj
BINDIR   := bin

# --- Flags ---
# Using immediate assignment (:=) for performance
BASE_FFLAGS := -fopenmp -J$(OBJDIR) -I$(OBJDIR)
RELEASE_FLG := -O3 -fbacktrace -fcheck=all -g -Warray-temporaries 
DEBUG_FLG   := -Og -g -Wall -Wextra -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -Warray-temporaries

# Set default FFLAGS to Release
FFLAGS      := $(BASE_FFLAGS) $(RELEASE_FLG)
LDFLAGS     := -llapack -lblas

# --- Files ---
TARGET      := $(BINDIR)/crashe
SRCS        := $(wildcard $(SRCDIR)/*.f90)
OBJS        := $(SRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)

# --- Rules ---
.PHONY: all clean debug dirs

all: dirs $(TARGET)

debug: FFLAGS := $(BASE_FFLAGS) $(DEBUG_FLG)
debug: dirs $(TARGET)

dirs:
	@$(MKDIR_P) $(OBJDIR) $(BINDIR)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

# --- Dependency Tree ---
# Simplified: only list what actually "USEs" what.
$(OBJDIR)/periodic_table.o:   $(OBJDIR)/types.o
$(OBJDIR)/input.o:   $(OBJDIR)/types.o

$(OBJDIR)/plasma_module.o:    $(OBJDIR)/types.o \
                              $(OBJDIR)/readadf04_module.o \
                              $(OBJDIR)/periodic_table.o

$(OBJDIR)/readadf04_module.o: $(OBJDIR)/types.o \
                              $(OBJDIR)/interpolation_module.o \
                              $(OBJDIR)/input.o 

$(OBJDIR)/crm_module.o:       $(OBJDIR)/types.o \
                              $(OBJDIR)/interpolation_module.o \
                              $(OBJDIR)/readadf04_module.o    \
                              $(OBJDIR)/plasma_module.o \
                              $(OBJDIR)/sorting.o

$(OBJDIR)/colradfort.o:       $(OBJDIR)/crm_module.o \
                              $(OBJDIR)/input.o 

$(OBJDIR)/onion_module.o:     $(OBJDIR)/colradfort.o


$(OBJDIR)/main.o:             $(OBJDIR)/colradfort.o\
                              $(OBJDIR)/onion_module.o

clean:
	$(RM) $(OBJDIR) $(BINDIR)

.PHONY: test

test: all
	python3 test.py