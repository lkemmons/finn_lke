.SUFFIXES: .o .f90

F90    = $(FC)
FFLAGS = -g

EXEC = finn_fire_emis

OBJS = finn_fire_emis.o

##dependencies
#$(OBJECTS) : makefile

${EXEC} :       ${OBJS}
		${F90} -o $@ ${OBJS} -L/usr/lib64 -lnetcdf -lnetcdff

.f90.o:
		${F90} ${FFLAGS} -I/usr/include -c $<
cleanup:
		rm -f ${OBJS}

clean:
		rm -f core ${EXEC} ${OBJS}
