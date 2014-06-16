#
# This is the makefile for compiling PETScExt. See 
# http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html for directions on installing PETSc.
# See also conf for additional commands.
#
ALL: all
LOCDIR	 = ./
DIRS	 = src include
CFLAGS	 = 
FFLAGS	 = 

INSTALL_DIR=${PETSC_DIR}/${PETSC_ARCH}

version:
	version_release=`grep '^#define PETSC_VERSION_RELEASE ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
	version_major=`grep '^#define PETSC_VERSION_MAJOR ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
	version_minor=`grep '^#define PETSC_VERSION_MINOR ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
	version_subminor=`grep '^#define PETSC_VERSION_SUBMINOR ' include/petscversion.h |tr -s ' ' | cut -d ' ' -f 3`; \
	echo ${version_release}

#include ${PETSC_DIR}/conf/base
#include ${PETSC_DIR}/conf/test
include ${PETSCEXT_DIR}/conf/petscext_common
#include ${PETSC_DIR}/conf/rules

#
# Basic targets to build PETSc libraries.
# all: builds the c, fortran, and f90 libraries
all: 
	-@echo "BUILDSHAREDLIB: ${BUILDSHAREDLIB}"
	-@${MAKE} clean_check
	@${OMAKE}  PETSC_ARCH=${PETSC_ARCH}  chkpetscdirpext
	-@${OMAKE} all_build 2>&1 | tee ${PETSC_ARCH}/conf/make.log
	-@echo "BUILDSHAREDLIB: ${BUILDSHAREDLIB}"
	-@egrep -i "( error | error:)" ${PETSC_ARCH}/conf/make.log > /dev/null; if [ "$$?" = "0" ]; then \
           echo "********************************************************************"; \
           echo "  Error during compile, check ${PETSC_ARCH}/conf/make.log"; \
           echo "  Send it and ${PETSC_ARCH}/conf/configure.log to dave.mayhem23@gmail.com";\
           echo "********************************************************************"; \
           exit 1; fi

all_build: chk_petsc_dir chk_petscext_dir chklib_dir info info_h deletelibs  build shared_pext

chkpetscdirpext:
	@mypwd=`pwd`; cd ${PETSC_DIR} 2>&1 > /dev/null; true_PETSC_DIR=`pwd`; cd $${mypwd} 2>&1 >/dev/null; \
        newpwd=`echo $${mypwd} | sed "s+$${true_PETSC_DIR}+DUMMY+g"`;\
        haspetsc=`echo $${mypwd} | sed "s+petsc-+DUMMY+g"`;\
        if [ $${mypwd} = $${newpwd} -a $${haspetsc} != $${mypwd} ]; then \
          printf ${PETSC_TEXT_HILIGHT}"*********************W-a-r-n-i-n-g*************************\n" ; \
          echo "Your PETSC_DIR may not match the directory you are in";\
          echo "PETSC_DIR " $${true_PETSC_DIR} "Current directory" $${mypwd};\
          echo "Ignore this if you are running make test             ";\
          printf "******************************************************"${PETSC_TEXT_NORMAL}"\n" ; \
        fi

#
# Prints information about the system and version of PETSc being compiled
#
info:
	-@echo "=========================================="
	-@echo On `date` on `hostname`
	-@echo Machine characteristics: `uname -a`
	-@echo "-----------------------------------------"
	-@echo "Using PETScExt directory ${PETSCEXT_DIR}"
	-@echo "Using PETSc directory: ${PETSC_DIR}"
	-@echo "Using PETSc arch: ${PETSC_ARCH}"
	-@echo "-----------------------------------------"
	-@grep "define PETSCEXT_VERSION" ${PETSCEXT_DIR}/include/petscextversion.h | ${SED} "s/........//"
	-@echo "-----------------------------------------"
	-@echo "Using configure Options: ${CONFIGURE_OPTIONS}"
	-@echo "Using configuration flags:"
	-@grep "\#define " ${INSTALL_DIR}/include/petscconf.h
	-@echo "-----------------------------------------"
	-@echo "Using include paths: ${PETSC_INCLUDE} ${PETSCEXT_INCLUDE}"
	-@echo "------------------------------------------"
	-@echo "Using C/C++ compiler: ${PCC} ${PCC_FLAGS} ${COPTFLAGS} ${CFLAGS}"
	-@echo "C/C++ Compiler version: " `${CCV}`
	-@if [ "${FC}" != "" ]; then \
	   echo "Using Fortran compiler: ${FC} ${FC_FLAGS} ${FFLAGS} ${FPP_FLAGS}";\
	   echo "Fortran Compiler version: " `${FCV}`;\
         fi
	-@echo "-----------------------------------------"
	-@echo "Using C/C++ linker: ${CC_LINKER}"
	-@if [ "${FC}" != "" ]; then \
	   echo "Using Fortran linker: ${FC_LINKER}";\
         fi
	-@echo "-----------------------------------------"
	-@echo "Using libraries: ${PETSC_LIB} ${PETSCEXT_LIB}"
	-@echo "------------------------------------------"
	-@echo "Using mpiexec: ${MPIEXEC}"
        -@echo "BUILDSHAREDLIB: ${BUILDSHAREDLIB}"
	-@echo "=========================================="
#
#
MINFO = ${INSTALL_DIR}/include/petscmachineinfo.h
info_h:
	-@cp ${MINFO} MINFO
	-@echo  "static const char *petscmachineinfo = \"\__n__\"" >> MINFO
	-@echo  "\"-----------------------------------------\__n__\"" >> MINFO
	-@if [ -f /usr/bin/cygcheck.exe ]; then \
	  echo  "\"Libraries compiled on `date` on `hostname|/usr/bin/dos2unix` \__n__\"" >> MINFO; \
          else \
	  echo  "\"Libraries compiled on `date` on `hostname` \__n__\"" >> MINFO; \
          fi
	-@echo  "\"Machine characteristics: `uname -a` \__n__\"" >> MINFO
	-@echo  "\"Using PETSc directory: ${PETSC_DIR}\__n__\"" >> MINFO
	-@echo  "\"Using PETSc arch: ${PETSC_ARCH}\__n__\"" >> MINFO
	-@echo  "\"-----------------------------------------\"; " >> MINFO
	-@echo  "static const char *petsccompilerinfo = \"\__n__\"" >> MINFO
	-@echo  "\"Using C compiler: ${PCC} ${PCC_FLAGS} ${COPTFLAGS} ${CFLAGS}\__n__\"" >> MINFO
	-@echo  "\"Using Fortran compiler: ${FC} ${FC_FLAGS} ${FFLAGS} ${FPP_FLAGS}\__n__\"" >> MINFO
	-@echo  "\"-----------------------------------------\"; " >> MINFO
	-@echo  "static const char *petsccompilerflagsinfo = \"\__n__\"" >> MINFO
	-@echo  "\"Using include paths: ${PETSC_INCLUDE}\__n__\"" >> MINFO
	-@echo  "\"------------------------------------------\"; " >> MINFO
	-@echo  "static const char *petsclinkerinfo = \"\__n__\"" >> MINFO
	-@echo  "\"Using C linker: ${CLINKER}\__n__\"" >> MINFO
	-@echo  "\"Using Fortran linker: ${FLINKER}\__n__\"" >> MINFO
	-@echo  "\"Using libraries: ${PETSC_LIB} \__n__\"" >> MINFO
	-@echo  "\"------------------------------------------\"; " >> MINFO
#	-@cat MINFO | ${SED} -e 's/\\ /\\\\ /g' | ${SED} -e 's/__n__/n/g' > ${MINFO}
	-@ if [ -f /usr/bin/cygcheck.exe ]; then /usr/bin/dos2unix ${MINFO} 2> /dev/null; fi
	-@rm -f MINFO

#
# Builds the PETSc libraries
# This target also builds fortran77 and f90 interface
# files and compiles .F files
#
build:
	-@echo "BEGINNING TO COMPILE PETScExt LIBRARIES IN ALL DIRECTORIES"
	-@echo "========================================="
	-@${OMAKE}  PETSC_ARCH=${PETSC_ARCH} PETSCEXT_DIR=${PETSCEXT_DIR} ACTION=libfast tree
	-@${RANLIB} ${PETSCEXT_LIB_DIR}/*.${AR_LIB_SUFFIX}
	-@echo "Completed building libraries"
	-@echo "========================================="

build_check:
	-@${MAKE} all -C src/vec/tests
	-@${MAKE} all -C src/tests

#	-@${MAKE} all -C src/mat/tests
#	-@${MAKE} all -C src/pc/tests
#	-@${MAKE} all -C src/ksp/tests
#	-@${MAKE} all -C src/snes/tests

check: build_check
	-@echo "======================================================"
	-@echo "BEGINNING TO COMPILE AND RUN PETScExt TEST EXAMPLES"
	-@echo "Due to different numerical round-off on certain"
	-@echo "machines some of the numbers may not match exactly."
	-@echo "======================================================"
	-@echo ""
	-@${MAKE} check -C src/vec/tests
	-@${MAKE} check -C src/tests
	-@echo ""
	-@echo "Completed compiling and running test examples"
	-@echo "======================================================"

clean_check:
	-@${MAKE} clean_test -C src/vec/tests
	-@${MAKE} clean_test -C src/mat/tests
	-@${MAKE} clean_test -C src/pc/tests
	-@${MAKE} clean_test -C src/ksp/tests
	-@${MAKE} clean_test -C src/snes/tests
	-@${MAKE} clean_test -C src/tests

#	-@${MAKE} check -C src/mat/tests
#	-@${MAKE} check -C src/pc/tests
#	-@${MAKE} check -C src/ksp/tests
#	-@${MAKE} check -C src/snes/tests






# Ranlib on the libraries
ranlib:
	${RANLIB} ${PETSCEXT_LIB_DIR}/*.${AR_LIB_SUFFIX}

# Deletes PETSc libraries
deletelibs: 
	-${RM} -rf ${PETSCEXT_LIB_DIR}/libpetscext*

# Cleans up build
allclean: deletelibs
	-@${OMAKE} PETSC_ARCH=${PETSC_ARCH} PETSCEXT_DIR=${PETSCEXT_DIR} ACTION=clean tree


#
# Check if PETSC_DIR variable specified is valid
#
chk_petsc_dir:
	@if [ ! -f ${PETSC_DIR}/include/petscversion.h ]; then \
	  echo "Incorrect PETSC_DIR specified: ${PETSC_DIR}!"; \
	  echo "You need to use / to separate directories, not \\!"; \
	  echo "Aborting build"; \
	  false; fi
#
#

#
# Check if PETSCEXT_DIR variable specified is valid
#
chk_petscext_dir:
	@if [ ! -f ${PETSCEXT_DIR}/include/petscextversion.h ]; then \
	  echo "Incorrect PETSCEXT_DIR specified: ${PETSCEXT_DIR}!"; \
	  echo "You need to use / to separate directories, not \\!"; \
	  echo "Aborting build"; \
	  false; fi
#
#

install:
	-@if [ "${PETSCEXT_INSTALL_DIR}" = "${PETSCEXT_DIR}" ]; then \
		echo "Install directory is current directory; nothing needs to be done";\
	else \
		echo Installing PETScExt at ${PETSCEXT_INSTALL_DIR};\
		echo "  1) Creating directories for installation";\
		if [ ! -d `dirname ${PETSCEXT_INSTALL_DIR}` ]; then \
			${MKDIR} `dirname ${PETSCEXT_INSTALL_DIR}` ; \
		fi;\
		if [ ! -d ${PETSCEXT_INSTALL_DIR} ]; then \
			${MKDIR} ${PETSCEXT_INSTALL_DIR} ; \
		fi;\
		echo "  2) Creating include/";\
		if [ ! -d ${PETSCEXT_INSTALL_DIR}/include ]; then \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/vec ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/mat ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/ksp ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/pc ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/snes ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/utils ; \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/include/private/helpers ; \
		fi;\
		find include/ -name "*.h" | xargs -I {} cp {} ${PETSCEXT_INSTALL_DIR}/{};\
		echo "  3) Creating conf/";\
		if [ ! -d ${PETSCEXT_INSTALL_DIR}/conf ]; then \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/conf ; \
		fi;\
		cp -f conf/petscext_common* ${PETSCEXT_INSTALL_DIR}/conf;\
		echo "  4) Creating PETSC_ARCH/";\
		if [ ! -d ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH} ]; then \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH} ; \
		fi;\
		echo "  5) Copying PETSC_ARCH/conf";\
		cp -f ${PETSC_ARCH}/petscextconf ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH}/petscextconf;\
		cp -rf ${PETSC_ARCH}/conf ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH}/;\
		if [ ! -d ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH}/lib ]; then \
			${MKDIR} ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH}/lib ; \
		fi;\
		if [ -d ${PETSC_ARCH}/lib ]; then \
			cp -f ${PETSC_ARCH}/lib/* ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH}/lib;\
			${RANLIB} ${PETSCEXT_INSTALL_DIR}/${PETSC_ARCH}/lib/*.a ;\
                        ${OMAKE} PETSC_ARCH=${PETSC_ARCH} PETSCEXT_DIR=${PETSCEXT_INSTALL_DIR} shared_pext; \
		fi;\
		echo "sh/bash: PETSCEXT_DIR="${PETSCEXT_INSTALL_DIR}"; export PETSCEXT_DIR";\
		echo "csh/tcsh: setenv PETSCEXT_DIR "${PETSCEXT_INSTALL_DIR} ;\
	fi;


install_src:
	-@if [ "${INSTALL_DIR}" = "${PETSC_DIR}" ]; then \
	  echo "You did not set a directory to install to";\
        else \
	  echo Installing PETSc source at ${INSTALL_DIR};\
          if [ ! -d `dirname ${INSTALL_DIR}` ]; then \
	    ${MKDIR} `dirname ${INSTALL_DIR}` ; \
          fi;\
          if [ ! -d ${INSTALL_DIR} ]; then \
	    ${MKDIR} ${INSTALL_DIR} ; \
          fi;\
          ${CP} -fr src ${INSTALL_DIR};\
        fi;


# ----------------------------------------------------

distclean:
	-@${MAKE} clean_check
	-@rm -rf ${PETSC_ARCH}

purge:
	-@${MAKE} clean_check
	-@rm -f config/*.pyc
	-@rm -f config/checklink.c
	-@rm -rf ${PETSC_ARCH}




petscao : petscmat petscao.f90.h
petscda : petscksp petscda.f90.h
petscdraw : petsc petscdraw.f90.h
petscis : petsc petscis.f90.h
petscksp : petscpc  petscksp.f90.h
petsclog : petsc petsclog.f90.h
petscmat : petscvec petscmat.f90.h
petscmg : petscksp petscmg.f90.h
petscpc : petscmat petscpc.f90.h
petscsnes : petscksp petscsnes.f90.h
petscsys : petsc petscsys.f90.h
petscts : petscsnes petscts.f90.h
petsc : petsc.f90.h
petscvec : petscis petscvec.f90.h
petscviewer : petsc petscviewer.f90.h
modules : petscao petscda petscdraw petscis petscksp petsclog petscmat petscmg petscpc petscsnes petscsys petscts petsc petscvec petscviewer
