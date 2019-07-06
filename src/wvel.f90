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


subroutine check_dt(U,P,rho,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables


! Make sure that the Courant criterion is always and everywhere fulfilled:
double precision :: wvel(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),dtnew
double precision :: P(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
integer :: i,j,k


! Every process does this:
dtnew=dt
do i=1,Nx
do j=1,Ny
do k=1,Nz
	call calc_wvel(wvel,i,j,k,U,P,rho,B)
	dtnew=min(dtnew,1.d0/(wvel(1,i,j,k)/dx+wvel(2,i,j,k)/dy+wvel(3,i,j,k)/dz))
enddo
enddo
enddo

! Now bring all processes together:
call mpi_allreduce(dtnew,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)

end subroutine



!ftp://ftp.sam.math.ethz.ch/pub/sam-reports/reports/reports2009/2009-33.pdf
subroutine calc_wvel(wvel,i,j,k,U,P,rho,B)
use globals
implicit none

integer :: i,j,k
double precision :: B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2), U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: rho(-1:Nx+2,-1:Ny+2,-1:Nz+2), P(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: wvel(3,-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: cs,vA

cs=sqrt(gamma*P(i,j,k)/rho(i,j,k))
vA=sqrt(B(1,i,j,k)*B(1,i,j,k)+B(2,i,j,k)*B(2,i,j,k)+B(3,i,j,k)*B(3,i,j,k))/sqrt(4.d0*M_PI*rho(i,j,k))

if(isnan(cs)) cs=0.d0

wvel(1,i,j,k)=sqrt((cs*cs+vA*vA)**2-4.d0*cs*cs*B(1,i,j,k)*B(1,i,j,k)/(4.d0*M_PI*rho(i,j,k)))
wvel(1,i,j,k)=sqrt(0.5d0*(wvel(1,i,j,k)+cs*cs+vA*vA))+abs(U(1,i,j,k))

wvel(2,i,j,k)=sqrt((cs*cs+vA*vA)**2-4.d0*cs*cs*B(2,i,j,k)*B(2,i,j,k)/(4.d0*M_PI*rho(i,j,k)))
wvel(2,i,j,k)=sqrt(0.5d0*(wvel(2,i,j,k)+cs*cs+vA*vA))+abs(U(2,i,j,k))

wvel(3,i,j,k)=sqrt((cs*cs+vA*vA)**2-4.d0*cs*cs*B(3,i,j,k)*B(3,i,j,k)/(4.d0*M_PI*rho(i,j,k)))
wvel(3,i,j,k)=sqrt(0.5d0*(wvel(3,i,j,k)+cs*cs+vA*vA))+abs(U(3,i,j,k))


end subroutine
