#	create executable

SRC  = main.c smod.c
LIBS = libnovas.a libcoord.a libgrvts.a
EXE  = fitgps.e

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

all: $(EXE)
$(EXE) : $(SRC) $(LIBS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

debug: $(DBG)
$(DBG) : $(SRC) $(LIBS)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LDFLAGS) $(LDLIBS)

clean:
	rm  $(EXE) $(DBG)
.PHONY : all debug clean
