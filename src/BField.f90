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


subroutine calc_dB(U,B,dB)
use globals
implicit none

double precision :: E(3,0:Nx,0:Ny,0:Nz),U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2),dB(3,Nx,Ny,Nz)
double precision :: Ucalc(3),Bcalc(3)
integer :: i,j,k

! You get U in the middle of the cells and B on the faces of the cells.
! Calculate what is E=-UxB on the edges of the cells!
! This is the staggered grid!
!

! Define E-Field:
do i=0,Nx
do j=0,Ny
do k=0,Nz
	Ucalc(2)=0.25d0*(U(2,i,j,k)+U(2,i,j+1,k)+U(2,i,j+1,k+1)+U(2,i,j,k+1))
	Ucalc(3)=0.25d0*(U(3,i,j,k)+U(3,i,j+1,k)+U(3,i,j+1,k+1)+U(3,i,j,k+1))
	Bcalc(2)=0.5d0*(B(2,i,j,k)+B(2,i,j,k+1))
	Bcalc(3)=0.5d0*(B(3,i,j,k)+B(3,i,j+1,k))
	E(1,i,j,k)=Bcalc(2)*Ucalc(3)-Bcalc(3)*Ucalc(2)


	Ucalc(1)=0.25d0*(U(1,i,j,k)+U(1,i+1,j,k)+U(1,i+1,j,k+1)+U(1,i,j,k+1) )
	Ucalc(3)=0.25d0*(U(3,i,j,k)+U(3,i+1,j,k)+U(3,i+1,j,k+1)+U(3,i,j,k+1) )
	Bcalc(1)=0.5d0*(B(1,i,j,k)+B(1,i,j,k+1))
	Bcalc(3)=0.5d0*(B(3,i,j,k)+B(3,i+1,j,k))
	E(2,i,j,k)=Bcalc(3)*Ucalc(1)-Bcalc(1)*Ucalc(3)


	Ucalc(1)=0.25d0*(U(1,i,j,k)+U(1,i,j+1,k)+U(1,i+1,j+1,k)+U(1,i+1,j,k))
	Ucalc(2)=0.25d0*(U(2,i,j,k)+U(2,i,j+1,k)+U(2,i+1,j+1,k)+U(2,i+1,j,k))
	Bcalc(1)=0.5d0*(B(1,i,j,k)+B(1,i,j+1,k))
	Bcalc(2)=0.5d0*(B(2,i,j,k)+B(2,i+1,j,k))
	E(3,i,j,k)=Bcalc(1)*Ucalc(2)-Bcalc(2)*Ucalc(1)
enddo
enddo
enddo


do i=1,Nx
do j=1,Ny
do k=1,Nz
	dB(1,i,j,k)=(E(2,i,j,k)-E(2,i,j,k-1))/dz-(E(3,i,j,k)-E(3,i,j-1,k))/dy
	dB(2,i,j,k)=(E(3,i,j,k)-E(3,i-1,j,k))/dx-(E(1,i,j,k)-E(1,i,j,k-1))/dz
	dB(3,i,j,k)=(E(1,i,j,k)-E(1,i,j-1,k))/dy-(E(2,i,j,k)-E(2,i-1,j,k))/dx
enddo
enddo
enddo


end subroutine







subroutine B_interpol(B,Binter,i,j,k)
use globals
implicit none

double precision :: B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2),Binter(3)
integer :: i,j,k



Binter(1)=0.5d0*(B(1,i,j,k)+B(1,i-1,j,k))

Binter(2)=0.5d0*(B(2,i,j,k)+B(2,i,j-1,k))

Binter(3)=0.5d0*(B(3,i,j,k)+B(3,i,j,k-1))


end subroutine
