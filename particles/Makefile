
CMP = mpic++
LNK = mpic++

INC_DCCRG = -I$(HOME)/AMR/dccrg

INC_BOOST = -I/usr/include/boost/
LIB_BOOST = -L/usr/lib/x86_64-linux-gnu/ -lboost_program_options

INC_ZOLTAN = -I$(HOME)/vlasiator/libraries/zoltan/src/include
LIB_ZOLTAN = -L$(HOME)/vlasiator/libraries/zoltan/build/src/ -lzoltan

INC_SAMRAI = -I$(HOME)/AMR/SAMRAI/build/include
LIB_SAMRAI = -L$(HOME)/AMR/SAMRAI/build/lib/ -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi/ -lSAMRAI_geom -lSAMRAI_mesh -lSAMRAI_pdat -lSAMRAI_xfer -lSAMRAI_hier -lSAMRAI_tbox -lhdf5

INC_CORGI = -I$(HOME)/AMR/corgi 
INC_CORGI += -I$(HOME)/AMR/corgi/mpi4cpp/include/
INC_CORGI += -I$(HOME)/AMR/corgi/mpi4cpp/tools/optional-lite/include/

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

corgi.o: corgi.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c corgi.cpp ${INCS} ${INC_CORGI} ${INC_DCCRG} $(CFLAGS)
corgi: corgi.o
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -o corgi corgi.o ${LIBS} $(LFLAGS)
	
samrai.o: samrai.cpp 
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -c samrai.cpp ${INCS} ${INC_SAMRAI} ${INC_DCCRG} $(CFLAGS)
samrai: samrai.o
	${CMP} ${CXXFLAGS} ${FLAGS} ${MATHFLAGS} -o samrai samrai.o ${LIBS} ${LIB_SAMRAI}