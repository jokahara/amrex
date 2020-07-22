
#set default architecture, can be overridden from the compile line
ARCH = ${VLASIATOR_ARCH}

include ../MAKE/Makefile.${ARCH}

INCS += ${INC_BOOST}
INCS += ${INC_ZOLTAN}

LIBS += ${LIB_BOOST}
LIBS += ${LIB_ZOLTAN}

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