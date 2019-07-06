! This file is part of HolmMHD
! Copyright (C) 2019 Daniel Verscharen (d.verscharen@ucl.ac.uk)
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
!ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
!WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
!ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
!(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
!LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
!ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!The views and conclusions contained in the software and documentation are those
!of the authors and should not be interpreted as representing official policies,
!either expressed or implied, of the HolmMHD project.

subroutine output(Nrun,rho,U,B,P,divB)
use globals
use HDF5
implicit none
include "mpif.h"      !Include MPI library variables

double precision :: rho(-1:Nx+2,-1:Ny+2,-1:Nz+2),P(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2), B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: divB(Nx,Ny,Nz)

double precision, allocatable, dimension(:,:,:) :: rhoout,Pout
double precision, allocatable, dimension(:,:,:,:) :: Uout,Bout

integer :: i,j,k,Nrun,error,Nout,comp,RBnumrun,RBi,RBj,RBk,num_RB,source
integer  :: status(MPI_STATUS_SIZE)

integer(HID_T) :: file_id		! File identifier
integer(HID_T) :: dset4_id,dset5_id  ! Dataset identifier
integer(HID_T) :: dset6_id,dset7_id,dset8_id       ! Dataset identifier
integer(HID_T) :: dspace4_id,dspace5_id ! Dataspace identifier
integer(HID_T) :: dspace6_id,dspace7_id,dspace8_id       ! Dataspace identifier
integer(HID_T) :: memspace4, memspace5, memspace6, memspace7, memspace8 ! Memory space identifier
integer(HSIZE_T) :: dims(3),dimsfield(4) ! Dataset dimensions
integer(HSIZE_T) :: count(3),offset(3) ! For hyperslab writing
integer(HSIZE_T) :: countfield(4),offsetfield(4) ! For hyperslab writing of fields

integer ::   rank=3, rankfield=4 ! Dataset rank
character(LEN=27) :: filename ! File name for the HDF5 file
character(LEN=4) :: timestring
character(LEN=3), PARAMETER :: dset4name = "rho"     ! Dataset name
character(LEN=1), PARAMETER :: dset5name = "U"     ! Dataset name
character(LEN=1), PARAMETER :: dset6name = "B"     ! Dataset name
character(LEN=1), PARAMETER :: dset7name = "P"     ! Dataset name
character(LEN=4), PARAMETER :: dset8name = "divB"     ! Dataset name


! Gather information from all root blocks:
num_RB=RBx*RBy*RBz


Nout=Nrun/outsteps
filename="./data/output_fields0000.h5"

if (Nout.LT.10) then
	 write ( timestring, '(i1)' ) Nout
	 filename(24:25)=timestring
	 filename(25:27)='.h5'
else if ((Nout.GE.10).AND.(Nout.LT.100)) then
	 write ( timestring, '(i2)' ) Nout
	 filename(23:24)=timestring
else if ((Nout.GE.100).AND.(Nout.LT.1000)) then
	 write ( timestring, '(i3)' ) Nout
	filename(22:24)=timestring
else if ((Nout.GE.1000).AND.(Nout.LT.10000)) then
	 write ( timestring, '(i4)' ) Nout
	filename(21:24)=timestring
else if (Nout.GE.10000) then
	write (*,*) "ERROR: More than 9999 timesteps."
endif

dims=(/ Nx*RBx,Ny*RBy,Nz*RBz /)
dimsfield=(/ 3,Nx*RBx,Ny*RBy,Nz*RBz /)

call h5open_f(error)

! The following part only if you are the master process:
if (master) then ! Create the HDF5 file:

	! Initialize HDF5 interface and create file:
	call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

	! Create the dataspace and dataset with default properties:
	call h5screate_simple_f(rank, dims, dspace4_id, error) ! for rho
	call h5dcreate_f(file_id, dset4name, H5T_NATIVE_DOUBLE, dspace4_id,dset4_id, error) ! for rho
	call h5screate_simple_f(rankfield, dimsfield, dspace5_id, error) ! for U
	call h5dcreate_f(file_id, dset5name, H5T_NATIVE_DOUBLE, dspace5_id,dset5_id, error) ! for U
	call h5screate_simple_f(rankfield, dimsfield, dspace6_id, error) ! for B
	call h5dcreate_f(file_id, dset6name, H5T_NATIVE_DOUBLE, dspace6_id,dset6_id, error) ! for B
	call h5screate_simple_f(rank, dims, dspace7_id, error) ! for P
	call h5dcreate_f(file_id, dset7name, H5T_NATIVE_DOUBLE, dspace7_id,dset7_id, error) ! for P
	call h5screate_simple_f(rank, dims, dspace8_id, error) ! for P
	call h5dcreate_f(file_id, dset8name, H5T_NATIVE_DOUBLE, dspace8_id,dset8_id, error) ! for divB

	call h5sclose_f(dspace4_id, error)
	call h5sclose_f(dspace5_id, error)
	call h5sclose_f(dspace6_id, error)
	call h5sclose_f(dspace7_id, error)
	call h5sclose_f(dspace8_id, error)
	call h5dclose_f(dset4_id, error)
	call h5dclose_f(dset5_id, error)
	call h5dclose_f(dset6_id, error)
	call h5dclose_f(dset7_id, error)
	call h5dclose_f(dset8_id, error)


	call h5fclose_f(file_id, error)
endif

call mpi_barrier(MPI_COMM_WORLD,ierror)

! This is done by all processes:

count=(/ Nx, Ny, Nz /)
countfield=(/ 3, Nx, Ny, Nz /)

! The following needs to happen sequentially:

do RBnumrun=1,num_RB
	call RBcoord(RBnumrun,RBi,RBj,RBk)

	offset=(/ (RBi-1)*Nx, (RBj-1)*Ny, (RBk-1)*Nz /)
	offsetfield=(/ 0, (RBi-1)*Nx, (RBj-1)*Ny, (RBk-1)*Nz /)

if (iproc.EQ.proc_for_block(RBnumrun)) then
! Open the file:
call h5fopen_f(filename,H5F_ACC_RDWR_F, file_id, error)

! Open the data sets:
call h5dopen_f(file_id, dset4name, dset4_id, error)
call h5dopen_f(file_id, dset5name, dset5_id, error)
call h5dopen_f(file_id, dset6name, dset6_id, error)
call h5dopen_f(file_id, dset7name, dset7_id, error)
call h5dopen_f(file_id, dset8name, dset8_id, error)

! Get data set identifiers and select subset:
call h5dget_space_f(dset4_id, dspace4_id, error)
call h5dget_space_f(dset5_id, dspace5_id, error)
call h5dget_space_f(dset6_id, dspace6_id, error)
call h5dget_space_f(dset7_id, dspace7_id, error)
call h5dget_space_f(dset8_id, dspace8_id, error)


call h5sselect_hyperslab_f(dspace4_id, H5S_SELECT_SET_F, offset, count, error)
call h5sselect_hyperslab_f(dspace5_id, H5S_SELECT_SET_F, offsetfield, countfield, error)
call h5sselect_hyperslab_f(dspace6_id, H5S_SELECT_SET_F, offsetfield, countfield, error)
call h5sselect_hyperslab_f(dspace7_id, H5S_SELECT_SET_F, offset, count, error)
call h5sselect_hyperslab_f(dspace8_id, H5S_SELECT_SET_F, offset, count, error)

! Create memory data space:
call h5screate_simple_f(rank, count, memspace4, error)
call h5screate_simple_f(rankfield, countfield, memspace5, error)
call h5screate_simple_f(rankfield, countfield, memspace6, error)
call h5screate_simple_f(rank, count, memspace7, error)
call h5screate_simple_f(rank, count, memspace8, error)


! Allocate temporary contiguous arrays and write subset to dataset:
allocate(rhoout(Nx,Ny,Nz))
rhoout(:,:,:)=rho(1:Nx,1:Ny,1:Nz)
call h5dwrite_f(dset4_id, H5T_NATIVE_DOUBLE, rhoout, count, error, memspace4, dspace4_id)
deallocate(rhoout)

allocate(Uout(3,Nx,Ny,Nz))
Uout(:,:,:,:)=U(1:3,1:Nx,1:Ny,1:Nz)
call h5dwrite_f(dset5_id, H5T_NATIVE_DOUBLE, Uout, countfield, error, memspace5, dspace5_id)
deallocate(Uout)

allocate(Bout(3,Nx,Ny,Nz))
Bout(:,:,:,:)=B(1:3,1:Nx,1:Ny,1:Nz)
call h5dwrite_f(dset6_id, H5T_NATIVE_DOUBLE, Bout, countfield, error, memspace6, dspace6_id)
deallocate(Bout)

allocate(Pout(Nx,Ny,Nz))
Pout(:,:,:)=P(1:Nx,1:Ny,1:Nz)
call h5dwrite_f(dset7_id, H5T_NATIVE_DOUBLE, Pout, count, error, memspace7, dspace7_id)
deallocate(Pout)

call h5dwrite_f(dset8_id, H5T_NATIVE_DOUBLE, divB, count, error, memspace8, dspace8_id)

! Close data spaces and data sets
call h5sclose_f(dspace4_id, error)
call h5sclose_f(dspace5_id, error)
call h5sclose_f(dspace6_id, error)
call h5sclose_f(dspace7_id, error)
call h5sclose_f(dspace8_id, error)
call h5sclose_f(memspace4, error)
call h5sclose_f(memspace5, error)
call h5sclose_f(memspace6, error)
call h5sclose_f(memspace7, error)
call h5sclose_f(memspace8, error)
call h5dclose_f(dset4_id, error)
call h5dclose_f(dset5_id, error)
call h5dclose_f(dset6_id, error)
call h5dclose_f(dset7_id, error)
call h5dclose_f(dset8_id, error)


! Close HDF5 interface:
call h5fclose_f(file_id, error)
call h5close_f(error)
endif

call mpi_barrier(MPI_COMM_WORLD,ierror)

enddo

end subroutine



subroutine generate_xdmf(Nrun,time_array)
use globals
implicit none
integer :: i,Nrun
character(LEN=14) :: form
character(LEN=1) :: xstring,ystring,zstring
double precision :: time_array(0:Nt+1)

open(unit=100,file='visualize_output.xmf',status='replace',action='write')

!define the filenames for the output
  write ( xstring, '(i1)' ) int(log10(1.*RBx*Nx))+2
  write ( ystring, '(i1)' ) int(log10(1.*RBy*Ny))+2
  write ( zstring, '(i1)' ) int(log10(1.*RBz*Nz))+1

  form(1:4)="(A,I"
  form(5:6)=zstring
  form(6:7)=",I"
  form(8:9)=ystring
  form(9:10)=",I"
  form(11:12)=xstring
  form(12:14)=",A)"


write (100,'(A)') '<?xml version="1.0" ?>'
write (100,'(A)') '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
write (100,'(A)') '<Xdmf Version="2.0">'
write (100,'(A)') '<Domain>'
write (100,'(A)') '<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'

do i=0,Nrun/outsteps
	write (100,'(A)') '  <Grid Name="XYZ" GridType="Uniform">'
	write (100,'(A,F15.5,A)') '  <Time Type="Single" Value="',time_array(i*outsteps),'" />'
	write (100,form) '    <Topology TopologyType="3DCoRectMesh" NumberOfElements="',RBz*Nz+1,RBy*Ny+1,RBx*Nx+1,'"/>'
	write (100,'(A)') '     <Geometry GeometryType="Origin_DxDyDz">'
	write (100,'(A)') '       <DataItem Name="Origin" DataType="Float" Dimensions="3" Format="XML">'
	write (100,'(A)')  '  0 0 0 '
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '       <DataItem Name="DxDyDz" DataType="Float" Dimensions="3" Format="XML">'
	write (100,*) '       ',dz,dy,dx
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Geometry>'
	write (100,'(A)') '     <Attribute Name="rho" AttributeType="Scalar" Center="Cell">'
	write (100,form) '       <DataItem Dimensions="',RBz*Nz,RBy*Ny,RBx*Nx,'" NumberType="Float" Precision="8" Format="HDF">'
	if (i.LT.10) then
		write (100,'(A,I1,A)') '        ./data/output_fields000',i,'.h5:/rho'
	else if ((i.GE.10).AND.(i.LT.100)) then
		write (100,'(A,I2,A)') '        ./data/output_fields00',i,'.h5:/rho'
	else if ((i.GE.100).AND.(i.LT.1000)) then
		write (100,'(A,I3,A)') '        ./data/output_fields0',i,'.h5:/rho'
	else if ((i.GE.1000).AND.(i.LT.10000)) then
		write (100,'(A,I4,A)') '        ./data/output_fields',i,'.h5:/rho'
	endif
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Attribute>'



		write (100,'(A)') '     <Attribute Name="U" AttributeType="Vector" Center="Cell">'
	write (100,form) '       <DataItem Dimensions="',RBz*Nz,RBy*Ny,RBx*Nx,' 3" NumberType="Float" Precision="8" Format="HDF">'
	if (i.LT.10) then
		write (100,'(A,I1,A)') '        ./data/output_fields000',i,'.h5:/U'
	else if ((i.GE.10).AND.(i.LT.100)) then
		write (100,'(A,I2,A)') '        ./data/output_fields00',i,'.h5:/U'
	else if ((i.GE.100).AND.(i.LT.1000)) then
		write (100,'(A,I3,A)') '        ./data/output_fields0',i,'.h5:/U'
	else if ((i.GE.1000).AND.(i.LT.10000)) then
		write (100,'(A,I4,A)') '        ./data/output_fields',i,'.h5:/U'
	endif
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Attribute>'

			write (100,'(A)') '     <Attribute Name="B" AttributeType="Vector" Center="Cell">'
	write (100,form) '       <DataItem Dimensions="',RBz*Nz,RBy*Ny,RBx*Nx,' 3" NumberType="Float" Precision="8" Format="HDF">'
	if (i.LT.10) then
		write (100,'(A,I1,A)') '        ./data/output_fields000',i,'.h5:/B'
	else if ((i.GE.10).AND.(i.LT.100)) then
		write (100,'(A,I2,A)') '        ./data/output_fields00',i,'.h5:/B'
	else if ((i.GE.100).AND.(i.LT.1000)) then
		write (100,'(A,I3,A)') '        ./data/output_fields0',i,'.h5:/B'
	else if ((i.GE.1000).AND.(i.LT.10000)) then
		write (100,'(A,I4,A)') '        ./data/output_fields',i,'.h5:/B'
	endif
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Attribute>'

	write (100,'(A)') '     <Attribute Name="P" AttributeType="Scalar" Center="Cell">'
	write (100,form) '       <DataItem Dimensions="',RBz*Nz,RBy*Ny,RBx*Nx,'" NumberType="Float" Precision="8" Format="HDF">'
	if (i.LT.10) then
		write (100,'(A,I1,A)') '        ./data/output_fields000',i,'.h5:/P'
	else if ((i.GE.10).AND.(i.LT.100)) then
		write (100,'(A,I2,A)') '        ./data/output_fields00',i,'.h5:/P'
	else if ((i.GE.100).AND.(i.LT.1000)) then
		write (100,'(A,I3,A)') '        ./data/output_fields0',i,'.h5:/P'
	else if ((i.GE.1000).AND.(i.LT.10000)) then
		write (100,'(A,I4,A)') '        output_fields',i,'.h5:/P'
	endif
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Attribute>'

		write (100,'(A)') '     <Attribute Name="div B" AttributeType="Scalar" Center="Cell">'
	write (100,form) '       <DataItem Dimensions="',RBz*Nz,RBy*Ny,RBx*Nx,'" NumberType="Float" Precision="8" Format="HDF">'
	if (i.LT.10) then
		write (100,'(A,I1,A)') '        ./data/output_fields000',i,'.h5:/divB'
	else if ((i.GE.10).AND.(i.LT.100)) then
		write (100,'(A,I2,A)') '        ./data/output_fields00',i,'.h5:/divB'
	else if ((i.GE.100).AND.(i.LT.1000)) then
		write (100,'(A,I3,A)') '        ./data/output_fields0',i,'.h5:/divB'
	else if ((i.GE.1000).AND.(i.LT.10000)) then
		write (100,'(A,I4,A)') '        output_fields',i,'.h5:/divB'
	endif
	write (100,'(A)') '       </DataItem>'
	write (100,'(A)') '     </Attribute>'

	write (100,'(A)') '   </Grid>'

enddo


write (100,'(A)') '</Grid>'
write (100,'(A)') ' </Domain>'
write (100,'(A)') '</Xdmf>'


close(100)

end subroutine
