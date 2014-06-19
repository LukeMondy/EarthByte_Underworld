include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules

checklink: checklink.o chkopts
	${CLINKER} -o checklink checklink.o ${TESTFLAGS} ${PETSC_KSP_LIB}
	@${RM} -f checklink checklink.o
