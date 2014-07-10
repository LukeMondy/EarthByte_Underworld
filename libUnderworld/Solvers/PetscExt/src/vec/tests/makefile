#################################################################################
#
#	Build the tests for PETScExt/vec
#
#################################################################################

include ${PETSCEXT2_DIR}/conf/petscext_common_variables

SRC		=	test-VecBlock.c \
			test-VecBlock-merge.c \
			test-VecBlock-mem.c

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
	-@${PETSCEXT2_DIR}/testing/ExecuteTests.pl


clean_test:
	-@rm -f ${APP}
	-@rm -f ${OBJ}
	-@${PETSCEXT2_DIR}/testing/ExecuteTests.pl clean
