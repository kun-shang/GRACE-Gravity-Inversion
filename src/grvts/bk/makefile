
MODS = alltide ataod

LIBS = libgrvts.a

###############################################################################
LDIR = ../../lib
LIB  = $(LDIR)/$(LIBS)

IDIR =../../include
FFLAGS   = -Wall -I$(IDIR) -O3
FC       = gfortran
CFLAGS   = -Wall -I$(IDIR) -O0 -g
CC       = gcc
LDFLAGS  = -lm
AR       = ar
ARFLAGS  = rcs

#OBJS = $(patsubst %,%.o,$(MODS))
OBJS = $(MODS:%=%.o)


all : $(LIB)

$(LIB) : $(LIB)($(OBJS))

clean:
	rm  $(LIB)
.PHONY : all clean


#%.o : %.c
#	$(CC) -c $(CFLAGS) $< -o $@
#	ar -rv $(LIB) *.o

#%.o : %.f
#	$(FC) -c $(FFLAGS) $< -o $@
#	ar -rv $(LIB) *.o

#%.o : %.f90
#	$(FC) -c $(FFLAGS) $< -o $@
#	ar -rv $(LIB) *.o
