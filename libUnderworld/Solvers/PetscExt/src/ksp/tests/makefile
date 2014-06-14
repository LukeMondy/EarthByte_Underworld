#################################################################################
#
#	Build the tests for PETScExt/mat
#
#################################################################################

include ${PETSCEXT_DIR}/conf/petscext_common_variables

SRC		=	test-ksp-cgr.c

OBJ		=	$(SRC:.c=.o)
APP		=	$(SRC:.c=.app)

INC		=	${PETSCEXT_INCLUDE}
LIB		=	${PETSCEXT_LIB}
CFLAGS		=	-g -Wall

all:
	@for a in ${APP} ; do \
		${MAKE} $$a; \
	done

%.app: %.o
	-@echo ""
	-${CLINKER} ${CFLAGS} -o $*.app $*.o ${INC} ${LIB}

%.o: %.c
	-@echo ""
	${PCC} ${CFLAGS} -c $*.c ${INC}

check:
	-@${PETSCEXT_DIR}/testing/ExecuteTests.pl


clean_test:
	-@rm -f ${APP}
	-@rm -f ${OBJ}
	-@${PETSCEXT_DIR}/testing/ExecuteTests.pl clean
