PROGRAM = domain
#PROGRAM = domain_dos

PROGRAMc = $(PROGRAM).c
PROGRAMo = $(PROGRAM).o

#SUBR = nrutil.c ricparam.c timer.c bulkgap_sbr.c OP_update.c op_guess.c step_integrators.c domainDOS.c
SUBR = nrutil.c ricparam_noPI.c mpi_timer.c bulkgap_sbr.c op_update.c op_guess.c step_integrators_triplet.c get_OP_FE.c FSsphere_lebedev.c

SUBRo = $(SUBR:.c=.o)

LIBS = -lm

FLAGS = -Wall -O
#FLAGS = -C -E -Wall -O 

#C = /usr/local/mpich-1.2.7/C/bin
#C = /usr/bin
#C = /usr/local/mpich-1.2.7p1/bin/
#Compile = /opt/local/bin/mpicc-mpich-mp
Compile = mpicc

.SUFFIXES : .c 
.c.o :
	$(Compile) $(FLAGS) -c $<

$(PROGRAM) : $(PROGRAMo) $(SUBRo)
	$(Compile) $(FLAGS) -o $(PROGRAM) $(PROGRAMo) $(SUBRo) $(LIBS)

#copy:
#	rm ../copy/*; cp * ../copy; rm *~
#
clean:
	rm *.o;
#
#data:
#	rm outdata.dat; rm sceq_evolution.dat; rm gap_evolution.dat
