FC	= /hosts/pluto/soft/intel_fc_80/bin/ifort
#FC	= g77
FFLAGS	= -O
SFX	= 
LIBLAPACK	= /usr/lib/liblapack.so.3
LIBBLAS =	/usr/lib/libblas.so.3

all:    pcazip.fort pcaunzip.fort pczdump.fort

pcazip.fort:
	$(FC) $(FFLAGS) -o pcazip.fort mygetargs.f90 trajutils.f90 x_io.f90 matfit.f  pcazip.f90 $(LIBLAPACK) $(LIBBLAS)

pczdump.fort:
	$(FC) $(FFLAGS) -o pczdump.fort mygetargs.f90 pcz_io.f90 pczdump.f90

pcaunzip.fort:
	$(FC) $(FFLAGS) -o pcaunzip.fort mygetargs.f90 pcz_io.f90 pcaunzip.f90

clean:
	rm -f pcazip.fort$(SFX) pcaunzip.fort$(SFX) pcadump.fort$(SFX)


