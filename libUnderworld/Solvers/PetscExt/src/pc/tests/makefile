#################################################################################
#
#	Build the tests for PETScExt/pc
#
#################################################################################

include ${PETSCEXT2_DIR}/conf/petscext_common_variables

SRC		=	test-PCBlock.c \
			test-PCBlock-Apply.c \
			test-PCBlock-block-ksp-solve.c \
			test-PCBlock-Apply-parallel.c \
			test-PCBlock-ApplyQ_QT.c \
			test-PCBlock-opts.c \
			test-PCBlock-Apply2.c \
			test-PCBlock-Apply3.c \
			test-PCBSchur-Apply.c

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

