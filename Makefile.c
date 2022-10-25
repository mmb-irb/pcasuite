include config.mk

SCANNER_OBJS=scanner/parser.tab.o scanner/lex.yy.o \
	scanner/ast.o scanner/aux_parser.o

all:    pcazip pczdump pcaunzip genpcz

# Linking targets
pcazip: pcazip.o traj_io.o binpos_io.o x_io.o matfit.o trajutils.o utils.o gaussianrms.o vector.o pcz_io.o netcdf_io.o $(SCANNER_OBJS)
	$(CC) $(LDFLAGS) -o pcazip pcazip.o traj_io.o binpos_io.o \
	x_io.o trajutils.o matfit.o netcdf_io.o \
	utils.o gaussianrms.o vector.o pcz_io.o $(SCANNER_OBJS) $(LIBLAPACK) $(LDLIBS)

pcazip_mpi: pcazip_mpi.o traj_io.o binpos_io.o x_io.o matfit.o trajutils.o utils.o gaussianrms.o vector.o pcz_io.o mpi_routines.o netcdf_io.o $(SCANNER_OBJS)
	$(MPICC) $(LDFLAGS) -o pcazip_mpi pcazip_mpi.o traj_io.o binpos_io.o \
	x_io.o trajutils.o matfit.o mpi_routines.o netcdf_io.o \
	utils.o gaussianrms.o vector.o pcz_io.o $(SCANNER_OBJS) $(LIBSCALAPACK) $(LDLIBS)
	
pczdump: pczdump.o utils.o pcz_io.o pczcompute.o traj_io.o netcdf_io.o binpos_io.o x_io.o matfit.o trajutils.o mahalanobis.o clusterlist.o $(SCANNER_OBJS)
pcaunzip: pcaunzip.o pcz_io.o utils.o pczcompute.o $(SCANNER_OBJS)
#pcaunzip: pcaunzip.o pcz_io.o utils.o pczcompute.o traj_io.o netcdf_io.o binpos_io.o x_io.o matfit.o trajutils.o $(SCANNER_OBJS)
genpcz: genpcz.o utils.o pcz_io.o

# Compilation targets
pcazip.o: pcazip.c pcazip.h traj_io.h binpos_io.h x_io.h trajutils.h constants.h pcz_io.h netcdf_io.h
pcazip_mpi.o: pcazip.c pcazip.h traj_io.h binpos_io.h x_io.h trajutils.h constants.h pcz_io.h netcdf_io.h
	$(MPICC) -c -o pcazip_mpi.o -DMPI_VERSION_ pcazip.c
mpi_routines.o: mpi_routines.c mpi_routines.c scalapack.h utils.h
	$(MPICC) -c -o mpi_routines.o $(FUNDERSCORES) -DMPI_VERSION_ mpi_routines.c
trajutils.o: trajutils.c trajutils.h pcazip.h
x_io.o: x_io.c x_io.h pcazip.h traj_io.h
utils.o: utils.c utils.h
#matfit.o: matfit.f
matfit.o: matfit.c
mahalanobis.o: mahalanobis.h pcz_io.h traj_io.h x_io.h trajutils.h
pczdump.o: pczdump.c pczdump.h constants.h traj_io.h binpos_io.h pcz_io.h scanner/y.tab.h clusterlist.h
pczcompute.o: pczcompute.c pczcompute.h utils.h
pcz_io.o: pcz_io.c pcz_io.h pcazip.h traj_io.h
netcdf_io.o: netcdf_io.c netcdf_io.h pcazip.h traj_io.h
gaussianrms.o: gaussianrms.c gaussianrms.h
vector.o: vector.c vector.h
clusterlist.o: clusterlist.c clusterlist.h utils.h
pcaunzip.o: pcaunzip.c pcaunzip.h utils.h constants.h pcz_io.h traj_io.h binpos_io.h pcz_io.h scanner/y.tab.h
genpcz.o: genpcz.c genpcz.h utils.h constants.h pcz_io.h
traj_io.o: traj_io.c traj_io.h x_io.h binpos_io.h netcdf_io.h pcazip.h utils.h

# Scanner targets
scanner/y.tab.h:
	(cd scanner; make y.tab.h; cd ..)

scanner/parser.tab.o:
	(cd scanner; make parser.tab.o; cd ..)

scanner/lex.yy.o:
	(cd scanner; make lex.yy.o; cd ..)

scanner/ast.o:
	(cd scanner; make ast.o; cd ..)

scanner/aux_parser.o: scanner/aux_parser.c scanner/aux_parser.h pcz_io.h
	(cd scanner; make aux_parser.o; cd ..)

# Cleaning target
.PHONY: clean
clean:
	rm -f pcazip$(SFX) pcaunzip$(SFX) pczdump$(SFX) core* *.o mpcc_*
	(cd scanner; make clean; cd ..)
