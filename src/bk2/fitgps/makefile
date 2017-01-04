#	create executable

SRC  = main.c smod.c
LIBS = libnovas.a libcoord.a libgrvts.a libnumrs.a
EXE  = fitgps.e
XDIR    = ../../exe
DDIR    = .

#	create executable

DIR = $(HOME)/my_packages
LDIR = $(DIR)/lib
IDIR = $(DIR)/include
VPATH = .:$(LDIR):

FFLAGS   = -Wall -I$(IDIR) -O3
FC       = gfortran
CFLAGS   = -Wall -I$(IDIR) -O3
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
$(DBG) : $(SRC) $(LIBS)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm  $(XDIR)/$(EXE) 
.PHONY : all debug clean
