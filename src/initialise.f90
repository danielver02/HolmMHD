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

subroutine initialise(rho,U,B,P)
use globals
implicit none
include "mpif.h"      !Include MPI library variables


integer :: i,j,k, RBi,RBj,RBk, RBnumrun, num_RB, ix,iy,iz, error, RBnum
double precision :: rho(-1:Nx+2,-1:Ny+2,-1:Nz+2),U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2),P(-1:Nx+2,-1:Ny+2,-1:Nz+2),r,Uampl,phi
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: x,y,z,xB,yB,zB
double precision :: kx,ky,kz,phase,Bamp
double precision :: Bsq,Usq,Binter(3)
integer,parameter :: seed = 86456


num_RB=RBx*RBy*RBz
RBnumrun=((iproc*num_RB)/nproc)+1
call RBcoord(RBnumrun,RBi,RBj,RBk)


!! Spherical hydro-blast wave:
!! Example on Mac: RBx=2, RBy=2, RBz=1, Nx=60, Ny=90, Nz=2
! dx=1.d0/(Nx*RBx)
! dy=1.5d0/(Ny*RBy)
!
! rho=1.d0
! gamma=5.d0/3.d0
! P=0.1d0
! U=0.d0
! B=0.d0
!
!
! do i=1,Nx
! do j=1,Ny
! 	x=dx*(i-1+(RBi-1)*Nx)
! 	y=dy*(j-1+(RBj-1)*Ny)
!
! 	if (((x-0.5d0)*(x-0.5d0)+(y-0.75d0)*(y-0.75d0)).LE.0.01)  P(i,j,:)=10.d0
!
! enddo
! enddo




! !Orszag-Tang vortex:
!! Example on Mac: RBx=2, RBy=2, RBz=1, Nx=90, Ny=90, Nz=2
! dx=1.d0/(RBx*Nx)
! dy=1.d0/(RBy*Ny)
! dz=0.01d0/(RBz*Nz)
!
! rho=25.d0/(36.d0*M_PI)
! P=5.d0/(12.d0*M_PI)
!
! U=0.d0
! B=0.d0
! ! Redefine step sizes:
! do i=1,Nx
! do j=1,Ny
!	x=dx*(i-1+(RBi-1)*Nx)
!	y=dy*(j-1+(RBj-1)*Ny)
!
! xB=x-0.5*dx
! yB=y-0.5*dy
!
! U(1,i,j,:)=-sin(2.d0*M_PI*y)
! U(2,i,j,:)=sin(2.d0*M_PI*x)
!
! B(1,i,j,:)=-sin(2.d0*M_PI*yB)*sqrt(4.d0*M_PI)
! B(2,i,j,:)=sin(4.d0*M_PI*xB)*sqrt(4.d0*M_PI)
! enddo
! enddo





!! Kelvin-Helmholtz instability:
!! Example on Mac: RBx=2, RBy=2, RBz=1, Nx=150, Ny=150, Nz=2
dx=1.d0/(RBx*Nx)
dy=1.d0/(RBy*Ny)
dz=0.01d0/(RBz*Nz)

P=2.5d0
gamma=1.4d0
rho=2.d0
U=0.d0
B=0.d0
B(1,:,:,:)=0.5d0
U(1,:,:,:)=0.5d0

call srand(seed)

do i=1,Nx
do j=1,Ny
do k=1,Nz
		x=dx*(i-1+(RBi-1)*Nx)
		y=dy*(j-1+(RBj-1)*Ny)
		z=dz*(k-1+(RBk-1)*Nz)

	 if ((y.GT.0.25d0).AND.(y.LT.0.75d0)) then
		 rho(i,j,k)=1.d0
		 U(1,i,j,k)=-0.5d0
	 endif

  ! add a random seed for the instability:

  U(1,i,j,k)=U(1,i,j,k)+0.1d0*(rand()-0.5d0)
  U(2,i,j,k)=U(2,i,j,k)+0.1d0*(rand()-0.5d0)
enddo
enddo
enddo




!! Alfven wave
!! Example on Mac: RBx=1, RBy=1, RBz=4, Nx=2, Ny=2, Nz=50
! dx=0.01d0/(RBx*Nx)
! dy=0.01d0/(RBy*Ny)
! dz=2.d0/(RBz*Nz)
!
! P=1.d0
! rho=1.d0
! U=0.d0
! B=0.d0
! B(3,:,:,:)=1.d0
! kz=2.d0*2.d0*M_PI/(dz*Nz*RBz)
! do k=1,Nz
! 	z=dz*(k-1+(RBk-1)*Nz)
! 	zB=z-0.5*dz
!
! 	B(1,:,:,k)=0.1d0*cos(kz*zB)*sqrt(4.d0*M_PI)
! 	B(2,:,:,k)=0.1d0*sin(kz*zB)*sqrt(4.d0*M_PI)
! 	U(1,:,:,k)=-0.1d0*cos(kz*z)
! 	U(2,:,:,k)=-0.1d0*sin(kz*z)
! enddo



!! Counter-propagating Alfven waves
!! Example on Mac: RBx=1, RBy=1, RBz=4, Nx=40, Ny=40, Nz=20
! dx=1.d0/(RBx*Nx)
! dy=1.d0/(RBy*Ny)
! dz=15.d0/(RBz*Nz)
!
! P=1.d0
! rho=1.d0
! U=0.d0
! B=0.d0
! B(3,:,:,:)=1.d0
! kx=2.d0*2.d0*M_PI/(dx*Nx*RBx)
! kz=2.d0*2.d0*M_PI/(dz*Nz*RBz)
! do i=1,Nx
! do j=1,Ny
! do k=1,Nz
! 	x=dx*(i-1+(RBi-1)*Nx)
! 	y=dy*(j-1+(RBj-1)*Ny)
! 	z=dz*(k-1+(RBk-1)*Nz)
! 	xB=x-0.5d0*dx
! 	yB=y-0.5d0*dy
! 	zB=z-0.5d0*dz
!
! 	B(2,i,j,k)=exp(-(zB-2.d0)**2/2.d0)*0.1d0*cos(kx*xB+kz*zB)*sqrt(4.d0*M_PI)
! 	U(2,i,j,k)=-exp(-(z-2.d0)**2/2.d0)*0.1d0*cos(kx*x+kz*z)
!
! 	B(1,i,j,k)=exp(-(zB-13.d0)**2/2.d0)*0.1d0*cos(kx*yB+kz*zB)*sqrt(4.d0*M_PI)
! 	U(1,i,j,k)=exp(-(z-13.d0)**2/2.d0)*0.1d0*cos(kx*y+kz*z)
! enddo
! enddo
! enddo





!! MHD aligned rotor:
!! Example on Mac: RBx=2, RBy=2, RBz=1, Nx=50, Ny=50, Nz=2
! dx=1.d0/(RBx*Nx)
! dy=1.d0/(RBy*Ny)
! dz=0.01d0/(RBz*Nz)
!
! P=1.d0
! rho=1.d0
! gamma=1.4d0
! B=0.d0
! U=0.d0
! B(1,:,:,:)=5.d0
!
! do i=1,Nx
! do j=1,Ny
! 	x=dx*(i-1+(RBi-1)*Nx)
! 	y=dy*(j-1+(RBj-1)*Ny)
!
! 	r=sqrt((x-0.5d0)*(x-0.5d0)+(y-0.5d0)*(y-0.5d0))
! 	if ((x-0.5d0).EQ.0.d0) then
! 		phi=0.5d0*M_PI
! 		if ((y-0.5d0).LT.0.d0) phi=-0.5d0*M_PI
! 	else
! 		phi=atan((y-0.5d0)/(x-0.5d0))
! 	endif
!
!
! 	if (r.LE.0.1d0) then
! 		rho(i,j,:)=10.d0
! 		U(1,i,j,:)=-2.d0*r*sin(phi)
! 		U(2,i,j,:)=2.d0*r*cos(phi)
! 		if ((x-0.5d0).LT.0.d0) then
! 			U(1,i,j,:)=-U(1,i,j,:)
! 			U(2,i,j,:)=-U(2,i,j,:)
! 		endif
! 	endif
!
! 	if ((r.LE.0.115d0).AND.(r.GT.0.1d0)) then
! 		rho(i,j,:)=10.d0-9.d0*(r-0.1d0)/0.015d0
! 		Uampl=(2.d0-1.d0*(r-0.1d0)/0.015d0)*0.1d0
! 		U(1,i,j,:)=-Uampl*sin(phi)
! 		U(2,i,j,:)=Uampl*cos(phi)
! 		if ((x-0.5d0).LT.0.d0) then
! 			U(1,i,j,:)=-U(1,i,j,:)
! 			U(2,i,j,:)=-U(2,i,j,:)
! 		endif
! 	endif
! enddo
! enddo



!!Isotropic turbulence initialisation:
!!Setup similar to Wan et al., PhPl 23, 042307, 2016
!! Example on Mac: RBx=2, RBy=2, RBz=2, Nx=30, Ny=30, Nz=30
! dx=1.d0/(RBx*Nx)
! dy=1.d0/(RBy*Ny)
! dz=1.d0/(RBz*Nz)
!
! P=1.d0
! rho=1.d0
! U=0.d0
! B=0.d0
! B(3,:,:,:)=1.d0
!
! do ix=-2,2
! do iy=-2,2
! do iz=-2,2
!
!
! kx=(2.d0*ix)*M_PI/(dx*Nx*RBx)
! ky=(2.d0*iy)*M_PI/(dy*Ny*RBy)
! kz=(2.d0*iz)*M_PI/(dz*Nz*RBz)
!
!
! if (master) phase=2.d0*M_PI*rand()
! call mpi_bcast(phase, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, error)
!
!
! if (sqrt(kx*kx+ky*ky+kz*kz).NE.0.d0) then
! 	Bamp=10.d0*(sqrt(kx*kx+ky*ky+kz*kz))**(-10.d0/3.d0)
! else
! 	Bamp=0.d0
! endif
!
! do i=1,Nx
! do j=1,Ny
! do k=1,Nz
! 	x=dx*(i-1+(RBi-1)*Nx)
! 	y=dy*(j-1+(RBj-1)*Ny)
! 	z=dz*(k-1+(RBk-1)*Nz)
! 	xB=x-0.5*dx
! 	yB=y-0.5*dy
!
! 	if (sqrt(kx*kx+ky*ky).NE.0.d0) then
! 		B(1,i,j,k)=B(1,i,j,k)-Bamp*cos(kx*xB+ky*y+kz*z-phase)*(ky/(sqrt(kx*kx+ky*ky)))*sqrt(4.d0*M_PI)
! 		B(2,i,j,k)=B(2,i,j,k)+Bamp*cos(kx*x+ky*yB+kz*z-phase)*(kx/(sqrt(kx*kx+ky*ky)))*sqrt(4.d0*M_PI)
! 		if (kz.GE.0.d0) then
! 			U(1,i,j,k)=U(1,i,j,k)+Bamp*cos(kx*x+ky*y+kz*z-phase)*(ky/(sqrt(kx*kx+ky*ky)))
! 			U(2,i,j,k)=U(2,i,j,k)-Bamp*cos(kx*x+ky*y+kz*z-phase)*(kx/(sqrt(kx*kx+ky*ky)))
! 		else
! 			U(1,i,j,k)=U(1,i,j,k)-Bamp*cos(kx*x+ky*y+kz*z-phase)*(ky/(sqrt(kx*kx+ky*ky)))
! 			U(2,i,j,k)=U(2,i,j,k)+Bamp*cos(kx*x+ky*y+kz*z-phase)*(kx/(sqrt(kx*kx+ky*ky)))
! 		endif
! 	elseif (kz.NE.0.d0) then
! 		B(1,i,j,k)=B(1,i,j,k)+Bamp*cos(kz*z-phase)*sqrt(4.d0*M_PI)
! 		if (kz.GT.0.d0) then
! 			U(1,i,j,k)=U(1,i,j,k)-Bamp*cos(kz*z-phase)
! 		else
! 			U(1,i,j,k)=U(1,i,j,k)+Bamp*cos(kz*z-phase)
! 		endif
! 	endif
!
! enddo
! enddo
! enddo
!
!
! enddo
! enddo
! enddo



end subroutine
