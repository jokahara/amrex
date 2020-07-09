
CMP = mpic++
LNK = mpic++

INC_DCCRG = -I$(HOME)/AMR/dccrg

INC_ZOLTAN = -I$(HOME)/vlasiator/libraries/zoltan/src/include
LIB_ZOLTAN = -L$(HOME)/vlasiator/libraries/zoltan/build/src/ -LIB_ZOLTAN

INCS += ${INC_BOOST}
INCS += ${INC_ZOLTAN}

LIBS += ${LIB_BOOST}
LIBS += ${LIB_ZOLTAN}

CXXFLAGS = -I $(HOME)/include  -L $(HOME)/lib -O3  -funroll-loops -std=c++17 -fopenmp -W -Wall -Wno-unused -fabi-version=0 -mavx -DPAPI_MEM
MATHFLAGS = -ffast-math

clean:
	rm -rf *.o dccrg amrex

dccrg.o: dccrg.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c dccrg.cpp ${INCS} ${INC_DCCRG}
dccrg: dccrg.o
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -o dccrg dccrg.o ${LIBS}

AMREX_LIBRARY_HOME ?= ${HOME}/AMR/amrex/tmp_install_dir/

LIBDIR := $(AMREX_LIBRARY_HOME)/lib
INCDIR := $(AMREX_LIBRARY_HOME)/include

COMPILE_CPP_FLAGS ?= $(shell awk '/Cflags:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)
COMPILE_LIB_FLAGS ?= $(shell awk '/Libs:/ {$$1=$$2=""; print $$0}' $(LIBDIR)/pkgconfig/amrex.pc)

CFLAGS := -I$(INCDIR) $(COMPILE_CPP_FLAGS)
LFLAGS := -L$(LIBDIR) $(COMPILE_LIB_FLAGS)

amrex.o: amrex.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c amrex.cpp ${INCS} ${INC_AMREX} ${INC_DCCRG} $(CFLAGS)
amrex: amrex.o
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -o amrex amrex.o ${LIBS} $(LFLAGS)
