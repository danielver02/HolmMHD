AM_FCFLAGS = -fcheck=all -fbacktrace -Wconversion -Wsurprising  -fbounds-check -ffixed-line-length-none -ffree-line-length-none -lm -I$(HDF5_INC)
AM_LDFLAGS = -ldl -L$(HDF5_LIB) $(HDF5_LDFLAGS)

bin_PROGRAMS = HolmMHD
HolmMHD_SOURCES = HolmMHD.f90 BField.f90 divF.f90 initialise.f90 mpi_routines.f90 output.f90 parameters.f90 wvel.f90
