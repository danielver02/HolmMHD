! HolmMHD is a numerical MHD code
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


module globals

integer :: Nx,Ny,Nz,Nt,outsteps
integer :: RBx,RBy,RBz
integer :: iproc,nproc,ierror

logical :: master

double precision :: dx,dy,dz
double precision :: dt,gamma
double precision :: M_PI=2.d0*asin(1.d0)


integer, allocatable,dimension(:) :: proc_for_block,communi_blocks,communi_procs
integer, allocatable,dimension(:) :: communi_edge_blocks,communi_edge_procs
integer, allocatable,dimension(:) :: communi_corner_blocks,communi_corner_procs
integer, allocatable,dimension(:,:) :: neighbours,edge_neighbours,corner_neighbours


end module




program HolmMHD
use globals
implicit none

include "mpif.h"      !Include MPI library variables

integer :: i,j,k,comp
integer :: Nrun
double precision, allocatable, dimension (:,:,:,:) :: U,Uhalf,continuity,momentumflux1,momentumflux2
double precision, allocatable, dimension (:,:,:,:) :: Bhalf,dB,wvel
double precision, allocatable, dimension (:,:,:,:) :: momentumflux3,B,Energflux
double precision, allocatable, dimension (:,:,:) :: rho, rhohalf,P,Phalf,Energ,Energhalf
double precision, allocatable, dimension (:,:,:) :: mom1,mom2,mom3,divB
double precision, allocatable, dimension (:) :: time_array

double precision :: Bsq,Usq,momentum(3),Binter(3)
double precision :: divF

logical :: run

call set_parameters()



! initialise MPI stuff:
  call mpi_init (ierror)
  call mpi_comm_size (mpi_comm_world, nproc, ierror)
  call mpi_comm_rank (mpi_comm_world, iproc, ierror)


master=.FALSE.
if (iproc.EQ.0) master=.TRUE.

if (master) write (*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
if (master) write (*,*) "                     This is HolmMHD 1.0                   "
if (master) write (*,*) "                      Daniel Verscharen                    "
if (master) write (*,*) "              Mullard Space Science Laboratory             "
if (master) write (*,*) "                  d.verscharen@ucl.ac.uk                   "
if (master) write (*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"
if (master) write (*,*) "Parameters:"
if (master) write (*,'(A5,I7,A15,I8,A15,F10.2)') " RBx:",RBx,"Nx:",Nx,"Lx:",dx*RBx*Nx
if (master) write (*,'(A5,I7,A15,I8,A15,F10.2)') " RBy:",RBy,"Ny:",Ny,"Ly:",dy*RBy*Ny
if (master) write (*,'(A5,I7,A15,I8,A15,F10.2)') " RBz:",RBz,"Nz:",Nz,"Lz:",dz*RBz*Nz
if (master) write (*,'(A4,I8,A15,I8,A15,F10.2)') " Nt:",Nt,"outsteps:",outsteps,"gamma:",gamma
if (master) write (*,*) "---"
if (master) write (*,*) "Splitting processes..."
call split_processes()

if (master) write (*,*) "---"

allocate(U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(Uhalf(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(momentumflux1(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(momentumflux2(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(momentumflux3(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(continuity(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(Energflux(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2))
allocate(Bhalf(3,-2:Nx+2,-2:Ny+2,-2:Nz+2))
allocate(dB(3,Nx,Ny,Nz))
allocate(rho(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(rhohalf(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(P(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(Phalf(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(Energ(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(Energhalf(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(mom1(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(mom2(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(mom3(-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(divB(Nx,Ny,Nz))
allocate(wvel(3,-1:Nx+2,-1:Ny+2,-1:Nz+2))
allocate(time_array(0:Nt+1))



run=.TRUE.
Nrun=0
divB=0.d0
if (master) write (*,*) "Initialising system..."
call initialise(rho,U,B,P)
call communicate(rho,energ,U,B)

rhohalf=rho
Uhalf=U
Bhalf=B
Phalf=P

! Initialise the energy density:
do i=-1,Nx+2
do j=-1,Ny+2
do k=-1,Nz+2
	call B_interpol(B,Binter,i,j,k)
	Bsq=(Binter(1)*Binter(1)+Binter(2)*Binter(2)+Binter(3)*Binter(3))/(8.d0*M_PI)
	Usq=U(1,i,j,k)*U(1,i,j,k)+U(2,i,j,k)*U(2,i,j,k)+U(3,i,j,k)*U(3,i,j,k)
	Energ(i,j,k)=P(i,j,k)/(gamma-1.d0)+Bsq+0.5d0*rho(i,j,k)*Usq
  if ((i.LE.Nx).AND.(i.GE.1).AND.(j.LE.Ny).AND.(j.GE.1).AND.(k.LE.Nz).AND.(k.GE.1)) then
    divB(i,j,k)=(B(1,i,j,k)-B(1,i-1,j,k))/dx+(B(2,i,j,k)-B(2,i,j-1,k))/dy+(B(3,i,j,k)-B(3,i,j,k-1))/dz
  endif
enddo
enddo
enddo
Energhalf=Energ
call communicate(rho,energ,U,B)



if (master) write (*,*) "---"

if (master) write (*,*) "Writing output..."
call output(Nrun,rho,U,B,P,divB)
if (master) call generate_xdmf(Nrun,time_array)
if (master) write (*,*) "---"

call mpi_barrier(MPI_COMM_WORLD,ierror)


if (master) write (*,*) "Start integration..."
if (master) write (*,*) "--------------------"
time_array(0)=0.d0

do while (run)  ! Run over all time steps

call check_dt(U,P,rho,B)

Nrun=Nrun+1
time_array(Nrun)=time_array(Nrun-1)+dt

if (master) write (*,'(A11,I10,A4,F9.3,A5,E12.3e2)') " Time step:", Nrun,", t=",time_array(Nrun),", dt=",dt


! Calculate the fluxes:
do i=-1,Nx+2
do j=-1,Ny+2
do k=-1,Nz+2

! Determine maximum wave speed
call calc_wvel(wvel,i,j,k,U,P,rho,B)

call B_interpol(B,Binter,i,j,k)
Bsq=(Binter(1)*Binter(1)+Binter(2)*Binter(2)+Binter(3)*Binter(3))/(8.d0*M_PI)

do comp=1,3
	continuity(comp,i,j,k)=rho(i,j,k)*U(comp,i,j,k)
	momentumflux1(comp,i,j,k)=rho(i,j,k)*U(comp,i,j,k)*U(1,i,j,k)&
		  -Binter(comp)*Binter(1)/(4.d0*M_PI)
	momentumflux2(comp,i,j,k)=rho(i,j,k)*U(comp,i,j,k)*U(2,i,j,k)&
		  -Binter(comp)*Binter(2)/(4.d0*M_PI)
	momentumflux3(comp,i,j,k)=rho(i,j,k)*U(comp,i,j,k)*U(3,i,j,k)&
		  -Binter(comp)*Binter(3)/(4.d0*M_PI)
	Energflux(comp,i,j,k)=(Energ(i,j,k)+P(i,j,k)+Bsq)*U(comp,i,j,k)&
		  -Binter(comp)*(Binter(1)*U(1,i,j,k)+Binter(2)*U(2,i,j,k)+Binter(3)*U(3,i,j,k))/(4.d0*M_PI)
enddo
	momentumflux1(1,i,j,k)=momentumflux1(1,i,j,k)+P(i,j,k)+Bsq
	momentumflux2(2,i,j,k)=momentumflux2(2,i,j,k)+P(i,j,k)+Bsq
	momentumflux3(3,i,j,k)=momentumflux3(3,i,j,k)+P(i,j,k)+Bsq
enddo
enddo
enddo





! Predictor-corrector method:
! half step in time:
call calc_dB(U,B,dB)
mom1=rho*U(1,:,:,:)
mom2=rho*U(2,:,:,:)
mom3=rho*U(3,:,:,:)

do i=1,Nx
do j=1,Ny
do k=1,Nz

  momentum(1)=rho(i,j,k)*U(1,i,j,k)-0.5d0*dt*divF(momentumflux1,mom1,i,j,k,wvel)
	momentum(2)=rho(i,j,k)*U(2,i,j,k)-0.5d0*dt*divF(momentumflux2,mom2,i,j,k,wvel)
	momentum(3)=rho(i,j,k)*U(3,i,j,k)-0.5d0*dt*divF(momentumflux3,mom3,i,j,k,wvel)

	rhohalf(i,j,k)=rho(i,j,k)-0.5d0*dt*divF(continuity,rho,i,j,k,wvel)
	Energhalf(i,j,k)=Energ(i,j,k)-0.5d0*dt*divF(Energflux,energ,i,j,k,wvel)

	Uhalf(1,i,j,k)=momentum(1)/rhohalf(i,j,k)
	Uhalf(2,i,j,k)=momentum(2)/rhohalf(i,j,k)
	Uhalf(3,i,j,k)=momentum(3)/rhohalf(i,j,k)

	Bhalf(1,i,j,k)=B(1,i,j,k)+0.5d0*dt*dB(1,i,j,k)
	Bhalf(2,i,j,k)=B(2,i,j,k)+0.5d0*dt*dB(2,i,j,k)
	Bhalf(3,i,j,k)=B(3,i,j,k)+0.5d0*dt*dB(3,i,j,k)

enddo
enddo
enddo


call communicate(rhohalf,energhalf,Uhalf,Bhalf)


! Calculate Phalf
do i=-1,Nx+2
do j=-1,Ny+2
do k=-1,Nz+2
	call B_interpol(Bhalf,Binter,i,j,k)
	Bsq=(Binter(1)*Binter(1)+Binter(2)*Binter(2)+Binter(3)*Binter(3))/(8.d0*M_PI)
	Usq=Uhalf(1,i,j,k)*Uhalf(1,i,j,k)+Uhalf(2,i,j,k)*Uhalf(2,i,j,k)+Uhalf(3,i,j,k)*Uhalf(3,i,j,k)
	Phalf(i,j,k)=(Energhalf(i,j,k)-Bsq-0.5d0*rhohalf(i,j,k)*Usq)*(gamma-1.d0)
enddo
enddo
enddo






! Calculate the fluxes in the half-step:
do i=-1,Nx+2
do j=-1,Ny+2
do k=-1,Nz+2

! Determine maximum wave speed
call calc_wvel(wvel,i,j,k,Uhalf,Phalf,rhohalf,Bhalf)

call B_interpol(Bhalf,Binter,i,j,k)
Bsq=(Binter(1)*Binter(1)+Binter(2)*Binter(2)+Binter(3)*Binter(3))/(8.d0*M_PI)


do comp=1,3
	continuity(comp,i,j,k)=rhohalf(i,j,k)*Uhalf(comp,i,j,k)
	momentumflux1(comp,i,j,k)=rhohalf(i,j,k)*Uhalf(comp,i,j,k)*Uhalf(1,i,j,k)&
		  -Binter(comp)*Binter(1)/(4.d0*M_PI)
	momentumflux2(comp,i,j,k)=rhohalf(i,j,k)*Uhalf(comp,i,j,k)*Uhalf(2,i,j,k)&
		  -Binter(comp)*Binter(2)/(4.d0*M_PI)
	momentumflux3(comp,i,j,k)=rhohalf(i,j,k)*Uhalf(comp,i,j,k)*Uhalf(3,i,j,k)&
		  -Binter(comp)*Binter(3)/(4.d0*M_PI)
	Energflux(comp,i,j,k)=(Energhalf(i,j,k)+Phalf(i,j,k)+Bsq)*Uhalf(comp,i,j,k)&
		  -Binter(comp)*(Binter(1)*Uhalf(1,i,j,k)+Binter(2)*Uhalf(2,i,j,k)+Binter(3)*Uhalf(3,i,j,k))/(4.d0*M_PI)
enddo
	momentumflux1(1,i,j,k)=momentumflux1(1,i,j,k)+Phalf(i,j,k)+Bsq
	momentumflux2(2,i,j,k)=momentumflux2(2,i,j,k)+Phalf(i,j,k)+Bsq
	momentumflux3(3,i,j,k)=momentumflux3(3,i,j,k)+Phalf(i,j,k)+Bsq
enddo
enddo
enddo


! second half step in time:
call calc_dB(Uhalf,Bhalf,dB)
mom1=rhohalf*Uhalf(1,:,:,:)
mom2=rhohalf*Uhalf(2,:,:,:)
mom3=rhohalf*Uhalf(3,:,:,:)

do i=1,Nx
do j=1,Ny
do k=1,Nz

	momentum(1)=rho(i,j,k)*U(1,i,j,k)-dt*divF(momentumflux1,mom1,i,j,k,wvel)
	momentum(2)=rho(i,j,k)*U(2,i,j,k)-dt*divF(momentumflux2,mom2,i,j,k,wvel)
	momentum(3)=rho(i,j,k)*U(3,i,j,k)-dt*divF(momentumflux3,mom3,i,j,k,wvel)

	rho(i,j,k)=rho(i,j,k)-dt*divF(continuity,rhohalf,i,j,k,wvel)
	Energ(i,j,k)=Energ(i,j,k)-dt*divF(Energflux,energhalf,i,j,k,wvel)

	U(1,i,j,k)=momentum(1)/rho(i,j,k)
	U(2,i,j,k)=momentum(2)/rho(i,j,k)
	U(3,i,j,k)=momentum(3)/rho(i,j,k)

	B(1,i,j,k)=B(1,i,j,k)+dt*dB(1,i,j,k)
	B(2,i,j,k)=B(2,i,j,k)+dt*dB(2,i,j,k)
	B(3,i,j,k)=B(3,i,j,k)+dt*dB(3,i,j,k)

enddo
enddo
enddo

call communicate(rho,energ,U,B)


! We can only interpolate the field after we know it. Therefore, we need an extra loop. We also require P:
do i=-1,Nx+2
do j=-1,Ny+2
do k=-1,Nz+2
	call B_interpol(B,Binter,i,j,k)
	Bsq=(Binter(1)*Binter(1)+Binter(2)*Binter(2)+Binter(3)*Binter(3))/(8.d0*M_PI)
	Usq=U(1,i,j,k)*U(1,i,j,k)+U(2,i,j,k)*U(2,i,j,k)+U(3,i,j,k)*U(3,i,j,k)
	P(i,j,k)=(Energ(i,j,k)-Bsq-0.5d0*rho(i,j,k)*Usq)*(gamma-1.d0)
	if ((i.LE.Nx).AND.(i.GE.1).AND.(j.LE.Ny).AND.(j.GE.1).AND.(k.LE.Nz).AND.(k.GE.1)) then
		divB(i,j,k)=(B(1,i,j,k)-B(1,i-1,j,k))/dx+(B(2,i,j,k)-B(2,i,j-1,k))/dy+(B(3,i,j,k)-B(3,i,j,k-1))/dz
	endif
enddo
enddo
enddo

if ((Nrun/outsteps).EQ.((1.d0*Nrun)/(1.d0*outsteps))) then
	if (master) write (*,*) "Writing output..."
	call output(Nrun,rho,U,B,P,divB)
  if (master) call generate_xdmf(Nrun,time_array)
	if (master) write (*,*) "---"
endif

if (isnan(P(1,1,1))) then
	if (master) write (*,*) "ERROR: NaN in pressure at origin."
	if (master) write (*,*) "Terminating..."
	run=.FALSE.
endif

if (Nrun.GE.Nt) run=.FALSE.
enddo

if (master) write (*,*) "---"
if (master) write (*,*) "Finalise MPI session..."
call mpi_finalize (ierror)

if (master) write (*,*) "Terminating HolmMHD."
if (master) write (*,*) "*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*"

end program
