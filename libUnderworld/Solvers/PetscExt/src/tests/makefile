#################################################################################
#
#	Build the tests for PETScExt/mat
#
#################################################################################

include ${PETSCEXT2_DIR}/conf/petscext_common_variables

SRC		=	stokes-operator.c \
			magma-operator.c \
			test-MatrixMarketIO.c \
			convert-petsc-MatVec-MatrixMarket.c \
			test-BlockSNESMFFD.c test-BlockMFFD.c \
			test-MatFreeExplicit.c

OBJ		=	$(SRC:.c=.o)
APP		=	$(SRC:.c=.app)

INC		=	${PETSCEXT_INCLUDE}
LIB		=	${PETSCEXT_LIB}
CFLAGS		=	-g -Wall

all:
	@for a in ${APP} ; do \
		${MAKE} $$a; \
	done

%.app: %.o petsccompat.o
	-@echo ""
	-${CLINKER} ${CFLAGS} -o $*.app $*.o petsccompat.o ${INC} ${LIB}

%.o: %.c ../compat/petsccompat.c
	-@echo ""
	${PCC} ${CFLAGS} -c $*.c ${INC}

check:
	-@${PETSCEXT2_DIR}/testing/ExecuteTests.pl


clean_test:
	-@rm -f ${APP}
	-@rm -f ${OBJ}
	-@rm -f petsc.mtx petsc_rhs1.mtx
	-@rm -f tmp-*
	-@${PETSCEXT2_DIR}/testing/ExecuteTests.pl clean

