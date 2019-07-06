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


double precision function divF(F,quant,i,j,k,wvel)
! Calculate the divergence of F at the point (i,j,k)
! F is a vector with F(comp,i,j,k)
use globals
implicit none

double precision :: F(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),quant(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: wvel(3,-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: fph,fmh,flimitp,flimitm,r,fpl,fml
double precision :: top,bot
integer :: i,j,k





! Hybrid scheme with flux limiter:


!!! DERIVATIVE IN X-DIRECTION:
! This is the van Leer flux limiter:
top=quant(i,j,k)-quant(i-1,j,k)
bot=quant(i+1,j,k)-quant(i,j,k)
r=top*bot/(bot*bot+1.d-10)
!flimitp=(r+abs(r))/(1.d0+abs(r))	! van Leer limiter
flimitp=max(0.d0,min(1.d0,r)) !min-mod limiter


top=quant(i-1,j,k)-quant(i-2,j,k)
bot=quant(i,j,k)-quant(i-1,j,k)
r=top*bot/(bot*bot+1.d-10)
!flimitm=(r+abs(r))/(1.d0+abs(r))
flimitm=max(0.d0,min(1.d0,r))


! This is the lower order part (Rusanov):
fpl=0.5d0*(F(1,i,j,k)+F(1,i+1,j,k))-0.5d0*max(wvel(1,i,j,k),wvel(1,i+1,j,k))*(quant(i+1,j,k)-quant(i,j,k))
fml=0.5d0*(F(1,i-1,j,k)+F(1,i,j,k))-0.5d0*max(wvel(1,i,j,k),wvel(1,i-1,j,k))*(quant(i,j,k)-quant(i-1,j,k))

! This is the higher order part (fourth-order centered):

fph=(7.d0/12.d0)*(F(1,i,j,k)+F(1,i+1,j,k)) - (1.d0/12.d0)*(F(1,i-1,j,k)+F(1,i+2,j,k))
fmh=(7.d0/12.d0)*(F(1,i-1,j,k)+F(1,i,j,k)) - (1.d0/12.d0)*(F(1,i-2,j,k)+F(1,i+1,j,k))
!divF=((fpl-0.5d0*flimitp*(fpl-fph))-(fml-0.5d0*flimitm*(fml-fmh)))/dx
divF=((fpl-flimitp*(fpl-fph))-(fml-flimitm*(fml-fmh)))/dx





!!! DERIVATIVE IN Y-DIRECTION:
! This is the van Leer flux limiter:
top=quant(i,j,k)-quant(i,j-1,k)
bot=quant(i,j+1,k)-quant(i,j,k)
r=top*bot/(bot*bot+1.d-10)
!flimitp=(r+abs(r))/(1.d0+abs(r))
flimitp=max(0.d0,min(1.d0,r))


top=quant(i,j-1,k)-quant(i,j-2,k)
bot=quant(i,j,k)-quant(i,j-1,k)
r=top*bot/(bot*bot+1.d-10)
!flimitm=(r+abs(r))/(1.d0+abs(r))
flimitm=max(0.d0,min(1.d0,r))


! This is the lower order part (Rusanov):
fpl=0.5d0*(F(2,i,j,k)+F(2,i,j+1,k))-0.5d0*max(wvel(2,i,j,k),wvel(2,i,j+1,k))*(quant(i,j+1,k)-quant(i,j,k))
fml=0.5d0*(F(2,i,j-1,k)+F(2,i,j,k))-0.5d0*max(wvel(2,i,j,k),wvel(2,i,j-1,k))*(quant(i,j,k)-quant(i,j-1,k))


! This is the higher order part(fourth-order centered):
fph=(7.d0/12.d0)*(F(2,i,j,k)+F(2,i,j+1,k)) - (1.d0/12.d0)*(F(2,i,j-1,k)+F(2,i,j+2,k))
fmh=(7.d0/12.d0)*(F(2,i,j-1,k)+F(2,i,j,k)) - (1.d0/12.d0)*(F(2,i,j-2,k)+F(2,i,j+1,k))
!divF=divF+((fpl-0.5d0*flimitp*(fpl-fph))-(fml-0.5d0*flimitm*(fml-fmh)))/dy
divF=divF+((fpl-flimitp*(fpl-fph))-(fml-flimitm*(fml-fmh)))/dy




!!! DERIVATIVE IN Z-DIRECTION:
! This is the van Leer flux limiter:

top=quant(i,j,k)-quant(i,j,k-1)
bot=quant(i,j,k+1)-quant(i,j,k)
r=top*bot/(bot*bot+1.d-10)
!flimitp=(r+abs(r))/(1.d0+abs(r))
flimitp=max(0.d0,min(1.d0,r))

top=quant(i,j,k-1)-quant(i,j,k-2)
bot=quant(i,j,k)-quant(i,j,k-1)
r=top*bot/(bot*bot+1.d-10)
!flimitm=(r+abs(r))/(1.d0+abs(r))
flimitm=max(0.d0,min(1.d0,r))


! This is the lower order part (Rusanov):
! We may possibly need a factor 0.5 instead of 0.25:
fpl=0.5d0*(F(3,i,j,k)+F(3,i,j,k+1))-0.5d0*max(wvel(3,i,j,k),wvel(3,i,j,k+1))*(quant(i,j,k+1)-quant(i,j,k))
fml=0.5d0*(F(3,i,j,k-1)+F(3,i,j,k))-0.5d0*max(wvel(3,i,j,k),wvel(3,i,j,k-1))*(quant(i,j,k)-quant(i,j,k-1))


! This is the higher order part(fourth-order centered):
fph=(7.d0/12.d0)*(F(3,i,j,k)+F(3,i,j,k+1)) - (1.d0/12.d0)*(F(3,i,j,k-1)+F(3,i,j,k+2))
fmh=(7.d0/12.d0)*(F(3,i,j,k-1)+F(3,i,j,k)) - (1.d0/12.d0)*(F(3,i,j,k-2)+F(3,i,j,k+1))
!divF=divF+((fpl-0.5d0*flimitp*(fpl-fph))-(fml-0.5d0*flimitm*(fml-fmh)))/dz
divF=divF+((fpl-flimitp*(fpl-fph))-(fml-flimitm*(fml-fmh)))/dz

end function
