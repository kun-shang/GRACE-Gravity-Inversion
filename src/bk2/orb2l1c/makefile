#	create executable

SRC  = main.c l1b2l1c.c
LIBS = libnovas.a libcoord.a libgrvts.a libnumrs.a
EXE  = orb2l1c.e
XDIR    = ../../exe
DDIR    = .

#	create executable

DIR = $(HOME)/my_packages
LDIR = $(DIR)/lib
IDIR = $(DIR)/include
VPATH = .:$(LDIR):

FFLAGS   = -Wall -I$(IDIR) -O3
FC       = gfortran
CFLAGS   = -Wall -I$(IDIR) -O3 -fopenmp
CC       = gcc
LDFLAGS  = -L. -L$(LDIR)
LDLIBS   = -lm -lgfortran
AR       = ar
ARFLAGS  = rcs

#EXE  = $(SRC:%.c=%.x)
#DBG  = $(SRC:%.c=%.d)

all: $(XDIR)/$(EXE)
$(XDIR)/$(EXE) : $(SRC) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

debug: $(DBG)
$(DDIR)/$(DBG) : $(SRC) $(LIBS)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm  $(XDIR)/$(EXE)
.PHONY : all debug clean
