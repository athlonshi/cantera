cc_sources = gri_pairs.cpp

INC                  = -I. -I$(top_builddir)/include/
AM_CPPFLAGS          = $(INC)
AM_CXXFLAGS          = $(AM_CPPFLAGS)

LINK                 = -loneD -lzeroD -lequil -lkinetics -ltransport -lthermo 
LINK                += -lctnumerics -lctmath -ltpx -lctspectra -lctbase -lcvode 
LINK                += -lctlapack -lctblas -lctf2c -lm  -lctf2c -lstdc++ -lexecstream -pthread
AM_LDFLAGS           = -L$(top_builddir)/build/lib/ 
LIBS                 =  $(LINK)

bin_PROGRAMS         = gri_pairs
library_includedir   = $(INC)

#-----------------------
# Cantera DH graph test
#-----------------------

gri_pairs_SOURCES = $(cc_sources)

TESTS_ENVIRONMENT = 
TESTS             = runtest
