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

subroutine split_processes()
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: num_RB,ideal_splitting,rest
integer :: i,j,k,l,n,RBnum,coord(3),block_ind
logical :: on_list

! Total number of root blocks:
num_RB=RBx*RBy*RBz
allocate(proc_for_block(num_RB))


if (num_RB.NE.nproc) then
	if (master) write (*,*) "Error: Number of root blocks and number of processes are not equal."
	call mpi_finalize (ierror)
	stop
endif

! ideal splitting (leaves a rest that should be handled by the last process):
ideal_splitting=num_RB/nproc
rest=num_RB-ideal_splitting*nproc



if(master) then
	write (*,*) "MPI Splitting:"
	write (*,'(A21,I14)') " Number of processes:",nproc
	write (*,'(A23,I12)') " Number of root blocks:",num_RB
endif


! Associate each block with a process:
do i=1,num_RB
	proc_for_block(i)=(i-1)/ideal_splitting
	if (i.GT.ideal_splitting*nproc) proc_for_block(i)=nproc-1
enddo


! Determine neighbours:
allocate(neighbours(num_RB,6))
call determine_neighbours()

! determine neighbours that this process iproc has to communicate with:
block_ind=1
allocate(communi_blocks(6*ideal_splitting+6))
communi_blocks=-1
do l=1,num_RB
 if (proc_for_block(l).EQ.iproc) then
  do n=1,6
   if(proc_for_block(neighbours(l,n)).NE.iproc) then
	communi_blocks(block_ind)=neighbours(l,n)
	block_ind=block_ind+1
   endif
  enddo
 endif
enddo



! Now find out which processes you as iproc have to talk to:
block_ind=1
allocate(communi_procs(nproc))
communi_procs=-1
do l=1,6*ideal_splitting+6
	if (communi_blocks(l).NE.-1) then
	! you have to communicate with proc_for_block(communi_blocks(l)) for sure.
	! Let's see if it is on the list already:
	on_list=.FALSE.
	do n=1,nproc-1
		if (proc_for_block(communi_blocks(l)).EQ.communi_procs(n)) on_list=.TRUE.
	enddo

	if (.NOT.on_list) then
		communi_procs(block_ind)=proc_for_block(communi_blocks(l))
		block_ind=block_ind+1
	endif
	endif
enddo

! sort both arrays:
call sort_array(communi_blocks,6*ideal_splitting+6)
call sort_array(communi_procs,nproc-1)

! Now we know that the process iproc has to communicate with all processes in the list
! communi_procs(nproc-1) about the blocks listed in communi_blocks(6*ideal_splitting+6).
! These are just the face neighbours of the block. We also have to know the neighbours behind
! the edges and in the corners of the block. That's the following part.




! Determine edge neighbours:
allocate(edge_neighbours(num_RB,12))
call determine_edge_neighbours()

! determine edge neighbours that this process iproc has to communicate with:
block_ind=1
allocate(communi_edge_blocks(12*ideal_splitting+12))
communi_edge_blocks=-1
do l=1,num_RB
 if (proc_for_block(l).EQ.iproc) then
  do n=1,12
   if(proc_for_block(edge_neighbours(l,n)).NE.iproc) then
	communi_edge_blocks(block_ind)=edge_neighbours(l,n)
	block_ind=block_ind+1
   endif
  enddo
 endif
enddo



! Now find out which processes you as iproc have to talk to:
block_ind=1
allocate(communi_edge_procs(nproc))
communi_edge_procs=-1
do l=1,12*ideal_splitting+12
	if (communi_edge_blocks(l).NE.-1) then
	! you have to communicate with proc_for_block(communi_blocks(l)) for sure.
	! Let's see if it is on the list already:
	on_list=.FALSE.
	do n=1,nproc-1
		if (proc_for_block(communi_edge_blocks(l)).EQ.communi_edge_procs(n)) on_list=.TRUE.
	enddo

	if (.NOT.on_list) then
		communi_edge_procs(block_ind)=proc_for_block(communi_edge_blocks(l))
		block_ind=block_ind+1
	endif
	endif
enddo

! sort both arrays:
call sort_array(communi_edge_blocks,12*ideal_splitting+12)
call sort_array(communi_edge_procs,nproc-1)




! Determine corner neighbours:
allocate(corner_neighbours(num_RB,8))
call determine_corner_neighbours()

! determine corner neighbours that this process iproc has to communicate with:
block_ind=1
allocate(communi_corner_blocks(8*ideal_splitting+8))
communi_corner_blocks=-1
do l=1,num_RB
 if (proc_for_block(l).EQ.iproc) then
  do n=1,8
   if(proc_for_block(corner_neighbours(l,n)).NE.iproc) then
	communi_corner_blocks(block_ind)=corner_neighbours(l,n)
	block_ind=block_ind+1
   endif
  enddo
 endif
enddo



! Now find out which processes you as iproc have to talk to:
block_ind=1
allocate(communi_corner_procs(nproc))
communi_corner_procs=-1
do l=1,8*ideal_splitting+8
	if (communi_corner_blocks(l).NE.-1) then
	! you have to communicate with proc_for_block(communi_blocks(l)) for sure.
	! Let's see if it is on the list already:
	on_list=.FALSE.
	do n=1,nproc-1
		if (proc_for_block(communi_corner_blocks(l)).EQ.communi_corner_procs(n)) on_list=.TRUE.
	enddo

	if (.NOT.on_list) then
		communi_corner_procs(block_ind)=proc_for_block(communi_corner_blocks(l))
		block_ind=block_ind+1
	endif
	endif
enddo

! sort both arrays:
call sort_array(communi_corner_blocks,8*ideal_splitting+8)
call sort_array(communi_corner_procs,nproc-1)

end subroutine








subroutine communicate(rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: l,i,n,elements
integer :: status(MPI_STATUS_SIZE)

logical :: on_list
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)


! How many elements are on your communi_procs list?
elements=0
do i=1,nproc-1
	if (communi_procs(i).NE.-1) elements=elements+1
enddo


! Consider the following cases for each root block (for the faces):
! from neighbour 1: you get the fields for (1:2,1:Ny,1:Nz)  and make them your (Nx+1:Nx+2,1:Ny,1:Nz)
! from neighbour 2: you get the fields for (Nx-1:Nx,1:Ny,1:Nz) and make them your (-1:0,1:Ny,1:Nz)
! from neighbour 3: you get the fields for (1:Nx,1:2,1:Nz)  and make them your (1:Nx,Ny+1:Ny+2,1:Nz)
! from neighbour 4: you get the fields for (1:Nx,Ny-1:Ny,1:Nz) and make them your (1:Nx,-1:0,1:Nz)
! from neighbour 5: you get the fields for (1:Nx,1:Ny,1:2)  and make them your (1:Nx,1:Ny,Nz+1:Nz+2)
! from neighbour 6: you get the fields for (1:Nx,1:Ny,Nz-1:Nz) and make them your (1:Nx,1:Ny,-1:0)
! For the B-field, there is one more at the "small" side



! receive from everything that has a lower process number than iproc
! then send to everything that has a higher process number than iproc
do l=0,nproc-1
	if (l.LT.iproc) then ! receive from l if l is in your list communi_procs:
		do i=1,elements
			if (l.EQ.communi_procs(i)) then
				call receive_from(l,rho,energ,U,B)
				call send_to(l,rho,energ,U,B)
			endif
		enddo
	elseif (l.GT.iproc) then ! send to l if l is in your list communi_procs:
		 do i=1,elements
		 	if (l.EQ.communi_procs(i)) then
		 		call send_to(l,rho,energ,U,B)
		 		call receive_from(l,rho,energ,U,B)
			endif
		enddo
	endif
enddo






! Now to the edges:

! How many elements are on your communi_edge_procs list?
elements=0
do i=1,nproc-1
	if (communi_edge_procs(i).NE.-1) elements=elements+1
enddo


! Consider the following cases for each root block (for the edges):
! from neighbour  1: you get the fields for (1,Ny,1:Nz)     and make them your (Nx+1,0,1:Nz)
! from neighbour  2: you get the fields for (1,1,1:Nz)      and make them your (Nx+1,Ny+1,1:Nz)
! from neighbour  3: you get the fields for (Nx,1,1:Nz)     and make them your (0,Ny+1,1:Nz)
! from neighbour  4: you get the fields for (Nx,Ny,1:Nz)    and make them your (0,0,1:Nz)
! from neighbour  5: you get the fields for (1,1:Ny,1)      and make them your (Nx+1,1:Ny,Nz+1)
! from neighbour  6: you get the fields for (1:Nx,1,1)      and make them your (1:Nx,Ny+1,Nz+1)
! from neighbour  7: you get the fields for (Nx,1:Ny,1)     and make them your (0,1:Ny,Nz+1)
! from neighbour  8: you get the fields for (1:Nx,Ny,1)     and make them your (1:Nx,0,Nz+1)
! from neighbour  9: you get the fields for (1,1:Ny,Nz)     and make them your (Nx+1,1:Ny,0)
! from neighbour 10: you get the fields for (1:Nx,1,Nz)     and make them your (1:Nx,Ny+1,0)
! from neighbour 11: you get the fields for (Nx,1:Ny,Nz)    and make them your (0,1:Ny,0)
! from neighbour 12: you get the fields for (1:Nx,Ny,Nz)    and make them your (1:Nx,0,0)
! For the B-field, there is one more at the "small" side



! receive from everything that has a lower process number than iproc
! then send to everything that has a higher process number than iproc
do l=0,nproc-1
	if (l.LT.iproc) then ! receive from l if l is in your list communi_procs:
		do i=1,elements
			if (l.EQ.communi_edge_procs(i)) then
				call receive_edge_from(l,rho,energ,U,B)
				call send_edge_to(l,rho,energ,U,B)
			endif
		enddo
	elseif (l.GT.iproc) then ! send to l if l is in your list communi_procs:
		 do i=1,elements
		 	if (l.EQ.communi_edge_procs(i)) then
		 		call send_edge_to(l,rho,energ,U,B)
		 		call receive_edge_from(l,rho,energ,U,B)
			endif
		enddo
	endif
enddo





! Now to the corners:

! How many elements are on your communi_corner_procs list?
elements=0
do i=1,nproc-1
	if (communi_corner_procs(i).NE.-1) elements=elements+1
enddo


! Consider the following cases for each root block (for the corners):
! from neighbour  1: you get the fields for (1,Ny,1)     and make them your (Nx+1,0,Nz+1)
! from neighbour  2: you get the fields for (1,1,1)      and make them your (Nx+1,Ny+1,Nz+1)
! from neighbour  3: you get the fields for (Nx,1,1)     and make them your (0,Ny+1,Nz+1)
! from neighbour  4: you get the fields for (Nx,Ny,1)    and make them your (0,0,Nz+1)
! from neighbour  5: you get the fields for (1,Ny,Nz)    and make them your (Nx+1,0,0)
! from neighbour  6: you get the fields for (1,1,Nz)     and make them your (Nx+1,Ny+1,0)
! from neighbour  7: you get the fields for (Nx,1,Nz)    and make them your (0,Ny+1,0)
! from neighbour  8: you get the fields for (Nx,Ny,Nz)   and make them your (0,0,0)



! receive from everything that has a lower process number than iproc
! then send to everything that has a higher process number than iproc
do l=0,nproc-1
	if (l.LT.iproc) then ! receive from l if l is in your list communi_procs:
		do i=1,elements
			if (l.EQ.communi_corner_procs(i)) then
				call receive_corner_from(l,rho,energ,U,B)
				call send_corner_to(l,rho,energ,U,B)
			endif
		enddo
	elseif (l.GT.iproc) then ! send to l if l is in your list communi_procs:
		 do i=1,elements
		 	if (l.EQ.communi_corner_procs(i)) then
		 		call send_corner_to(l,rho,energ,U,B)
		 		call receive_corner_from(l,rho,energ,U,B)
			endif
		enddo
	endif
enddo





! Am I my own neighbour (for the faces)?
do i=1,3
 if (neighbours(iproc+1,2*i).EQ.(iproc+1)) then
 	select case (i)
	case (1) ! equivalent to case 2
 			rho(Nx+1:Nx+2,1:Ny,1:Nz)=rho(1:2,1:Ny,1:Nz)
 			rho(-1:0,1:Ny,1:Nz)=rho(Nx-1:Nx,1:Ny,1:Nz)
 			energ(Nx+1:Nx+2,1:Ny,1:Nz)=energ(1:2,1:Ny,1:Nz)
 			energ(-1:0,1:Ny,1:Nz)=energ(Nx-1:Nx,1:Ny,1:Nz)
 			U(:,Nx+1:Nx+2,1:Ny,1:Nz)=U(:,1:2,1:Ny,1:Nz)
 			U(:,-1:0,1:Ny,1:Nz)=U(:,Nx-1:Nx,1:Ny,1:Nz)
 			B(:,Nx+1:Nx+2,1:Ny,1:Nz)=B(:,1:2,1:Ny,1:Nz)
 			B(:,-2:0,1:Ny,1:Nz)=B(:,Nx-2:Nx,1:Ny,1:Nz)
 		case (2) ! equivalent to case 4
 		 	rho(1:Nx,Ny+1:Ny+2,1:Nz)=rho(1:Nx,1:2,1:Nz)
 			rho(1:Nx,-1:0,1:Nz)=rho(1:Nx,Ny-1:Ny,1:Nz)
 			energ(1:Nx,Ny+1:Ny+2,1:Nz)=energ(1:Nx,1:2,1:Nz)
 			energ(1:Nx,-1:0,1:Nz)=energ(1:Nx,Ny-1:Ny,1:Nz)
 			U(:,1:Nx,Ny+1:Ny+2,1:Nz)=U(:,1:Nx,1:2,1:Nz)
 			U(:,1:Nx,-1:0,1:Nz)=U(:,1:Nx,Ny-1:Ny,1:Nz)
 			B(:,1:Nx,Ny+1:Ny+2,1:Nz)=B(:,1:Nx,1:2,1:Nz)
 			B(:,1:Nx,-2:0,1:Nz)=B(:,1:Nx,Ny-2:Ny,1:Nz)
 		case (3) ! equivalent to case 6
 		  rho(1:Nx,1:Ny,Nz+1:Nz+2)=rho(1:Nx,1:Ny,1:2)
 			rho(1:Nx,1:Ny,-1:0)=rho(1:Nx,1:Ny,Nz-1:Nz)
 			energ(1:Nx,1:Ny,Nz+1:Nz+2)=energ(1:Nx,1:Ny,1:2)
 			energ(1:Nx,1:Ny,-1:0)=energ(1:Nx,1:Ny,Nz-1:Nz)
 			U(:,1:Nx,1:Ny,Nz+1:Nz+2)=U(:,1:Nx,1:Ny,1:2)
 			U(:,1:Nx,1:Ny,-1:0)=U(:,1:Nx,1:Ny,Nz-1:Nz)
 			B(:,1:Nx,1:Ny,Nz+1:Nz+2)=B(:,1:Nx,1:Ny,1:2)
 			B(:,1:Nx,1:Ny,-2:0)=B(:,1:Nx,1:Ny,Nz-2:Nz)
 	end select
 endif
enddo



! Am I my own neighbour (for the edges)?
do i=1,8
 if (edge_neighbours(iproc+1,i).EQ.(iproc+1)) then
	select case (i)
		case (1)
			rho(Nx+1,0,1:Nz)=rho(1,Ny,1:Nz)
			energ(Nx+1,0,1:Nz)=energ(1,Ny,1:Nz)
			U(:,Nx+1,0,1:Nz)=U(:,1,Ny,1:Nz)
			B(:,Nx+1,-1:0,1:Nz)=B(:,1,Ny-1:Ny,1:Nz)
			rho(0,Ny+1,1:Nz)=rho(Nx,1,1:Nz) !This is equivalent to case 3
			energ(0,Ny+1,1:Nz)=energ(Nx,1,1:Nz)
			U(:,0,Ny+1,1:Nz)=U(:,Nx,1,1:Nz)
			B(:,-1:0,Ny+1,1:Nz)=B(:,Nx-1:Nx,1,1:Nz)
		case (2)
			rho(Nx+1,Ny+1,1:Nz)=rho(1,1,1:Nz)
			energ(Nx+1,Ny+1,1:Nz)=energ(1,1,1:Nz)
			U(:,Nx+1,Ny+1,1:Nz)=U(:,1,1,1:Nz)
			B(:,Nx+1,Ny+1,1:Nz)=B(:,1,1,1:Nz)
			rho(0,0,1:Nz)=rho(Nx,Ny,1:Nz) ! This is equivalent to case 4
			energ(0,0,1:Nz)=energ(Nx,Ny,1:Nz)
			U(:,0,0,1:Nz)=U(:,Nx,Ny,1:Nz)
			B(:,-1:0,-1:0,1:Nz)=B(:,Nx-1:Nx,Ny-1:Ny,1:Nz)
		case (5)
			rho(Nx+1,1:Ny,Nz+1)=rho(1,1:Ny,1)
			energ(Nx+1,1:Ny,Nz+1)=energ(1,1:Ny,1)
			U(:,Nx+1,1:Ny,Nz+1)=U(:,1,1:Ny,1)
			B(:,Nx+1,1:Ny,Nz+1)=B(:,1,1:Ny,1)
			rho(0,1:Ny,0)=rho(Nx,1:Ny,Nz) ! This is equivalent to case 11
			energ(0,1:Ny,0)=energ(Nx,1:Ny,Nz)
			U(:,0,1:Ny,0)=U(:,Nx,1:Ny,Nz)
			B(:,-1:0,1:Ny,-1:0)=B(:,Nx-1:Nx,1:Ny,Nz-1:Nz)
		case (6)
			rho(1:Nx,Ny+1,Nz+1)=rho(1:Nx,1,1)
			energ(1:Nx,Ny+1,Nz+1)=energ(1:Nx,1,1)
			U(:,1:Nx,Ny+1,Nz+1)=U(:,1:Nx,1,1)
			B(:,1:Nx,Ny+1,Nz+1)=B(:,1:Nx,1,1)
			rho(1:Nx,0,0)=rho(1:Nx,Ny,Nz) ! This is equivalent to case 12
			energ(1:Nx,0,0)=energ(1:Nx,Ny,Nz)
			U(:,1:Nx,0,0)=U(:,1:Nx,Ny,Nz)
			B(:,1:Nx,-1:0,-1:0)=B(:,1:Nx,Ny-1:Ny,Nz-1:Nz)
		case (7)
			rho(0,1:Ny,Nz+1)=rho(Nx,1:Ny,1)
			energ(0,1:Ny,Nz+1)=energ(Nx,1:Ny,1)
			U(:,0,1:Ny,Nz+1)=U(:,Nx,1:Ny,1)
			B(:,-1:0,1:Ny,Nz+1)=B(:,Nx-1:Nx,1:Ny,1)
			rho(Nx+1,1:Ny,0)=rho(1,1:Ny,Nz) ! This is equivalent to case 9
			energ(Nx+1,1:Ny,0)=energ(1,1:Ny,Nz)
			U(:,Nx+1,1:Ny,0)=U(:,1,1:Ny,Nz)
			B(:,Nx+1,1:Ny,-1:0)=B(:,1,1:Ny,-1:Nz)
		case (8)
			rho(1:Nx,0,Nz+1)=rho(1:Nx,Ny,1)
			energ(1:Nx,0,Nz+1)=energ(1:Nx,Ny,1)
			U(:,1:Nx,0,Nz+1)=U(:,1:Nx,Ny,1)
			B(:,1:Nx,-1:0,Nz+1)=B(:,1:Nx,Ny-1:Ny,1)
			rho(1:Nx,Ny+1,0)=rho(1:Nx,1,Nz) ! This is equivalent to case 10
			energ(1:Nx,Ny+1,0)=energ(1:Nx,1,Nz)
			U(:,1:Nx,Ny+1,0)=U(:,1:Nx,1,Nz)
			B(:,1:Nx,Ny+1,-1:0)=B(:,1:Nx,1,Nz-1:Nz)
	end select
 endif
enddo


! Am I my own neighbour (for the corners)?
do i=1,4
 if (corner_neighbours(iproc+1,i).EQ.(iproc+1)) then
	select case (i)
		case (1)
			rho(Nx+1,0,Nz+1)=rho(1,Ny,1)
			energ(Nx+1,0,Nz+1)=energ(1,Ny,1)
			U(:,Nx+1,0,Nz+1)=U(:,1,Ny,1)
			B(:,Nx+1,0,Nz+1)=B(:,1,Ny,1)
			rho(0,Ny+1,0)=rho(Nx,1,Nz) ! This is equivalent to case 7
			energ(0,Ny+1,0)=energ(Nx,1,Nz)
			U(:,0,Ny+1,0)=U(:,Nx,1,Nz)
			B(:,0,Ny+1,0)=B(:,Nx,1,Nz)
		case (2)
			rho(Nx+1,Ny+1,Nz+1)=rho(1,1,1)
			energ(Nx+1,Ny+1,Nz+1)=energ(1,1,1)
			U(:,Nx+1,Ny+1,Nz+1)=U(:,1,1,1)
			B(:,Nx+1,Ny+1,Nz+1)=B(:,1,1,1)
			rho(0,0,0)=rho(Nx,Ny,Nz) ! This is equivalent to case 8
			energ(0,0,0)=energ(Nx,Ny,Nz)
			U(:,0,0,0)=U(:,Nx,Ny,Nz)
			B(:,0,0,0)=B(:,Nx,Ny,Nz)
		case (3)
			rho(0,Ny+1,Nz+1)=rho(Nx,1,1)
			energ(0,Ny+1,Nz+1)=energ(Nx,1,1)
			U(:,0,Ny+1,Nz+1)=U(:,Nx,1,1)
			B(:,0,Ny+1,Nz+1)=B(:,Nx,1,1)
			rho(Nx+1,0,0)=rho(1,Ny,Nz) ! This is equivalent to case 5
			energ(Nx+1,0,0)=energ(1,Ny,Nz)
			U(:,Nx+1,0,0)=U(:,1,Ny,Nz)
			B(:,Nx+1,0,0)=B(:,1,Ny,Nz)
		case (4)
			rho(0,0,Nz+1)=rho(Nx,Ny,1)
			energ(0,0,Nz+1)=energ(Nx,Ny,1)
			U(:,0,0,Nz+1)=U(:,Nx,Ny,1)
			B(:,0,0,Nz+1)=B(:,Nx,Ny,1)
			rho(Nx+1,Ny+1,0)=rho(1,1,Nz) ! This is equivalent to case 6
			energ(Nx+1,Ny+1,0)=energ(1,1,Nz)
			U(:,Nx+1,Ny+1,0)=U(:,1,1,Nz)
			B(:,Nx+1,Ny+1,0)=B(:,1,1,Nz)
	end select
 endif
enddo

call mpi_barrier(MPI_COMM_WORLD,ierror)
end subroutine




subroutine send_to(proc,rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: proc,communicator,i,proc_is_neighbour,neighbour_list(6),tag,j
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision, allocatable, dimension (:,:,:) :: send_sc
double precision, allocatable, dimension (:,:,:,:) :: send_vec



do i=1,6
neighbour_list = (/1,3,5,2,4,6/) ! The order in which the neighbours are handled
j=neighbour_list(i)
proc_is_neighbour=0
if (neighbours(iproc+1,j).EQ.(proc+1)) proc_is_neighbour = j

	!if((iproc.EQ.0).AND.(proc_is_neighbour.ne.0)) write (*,*) "This is",iproc,"sending to",proc,"as neighbour",proc_is_neighbour

! Consider the following cases for each root block:
! to neighbour 1: you send your (Nx-1:Nx,1:Ny,1:Nz) and make them their (-1:0,1:Ny,1:Nz)
! to neighbour 2: you send your (1:2,1:Ny,1:Nz) and make them their (Nx+1:Nx+2,1:Ny,1:Nz)
! to neighbour 3: you send your (1:Nx,Ny-1:Ny,1:Nz) and make them their (1:Nx,-1:0,1:Nz)
! to neighbour 4: you send your (1:Nx,1:2,1:Nz) and make them their (1:Nx,Ny+1:Ny+2,1:Nz)
! to neighbour 5: you send your (1:Nx,1:Ny,Nz-1:Nz) and make them their (1:Nx,1:Ny,-1:0)
! to neighbour 6: you send your (1:Nx,1:Ny,1:2) and make them their (1:Nx,1:Ny,Nz+1:Nz+2)
! For the B-field, there is one more at the "small" side

! tag system:
! tag = TN
! first digit: T: type - 1: rho, 2: energ, 3: U, 4: B
! second digit: N: neighbour number of the target

select case (proc_is_neighbour)
	case (1)
		allocate(send_sc(2,Ny,Nz))
			send_sc=rho(Nx-1:Nx,1:Ny,1:Nz)
			tag=11
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx-1:Nx,1:Ny,1:Nz)
			tag=21
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,2,Ny,Nz))
			send_vec=U(:,Nx-1:Nx,1:Ny,1:Nz)
			tag=31
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec(3,3,Ny,Nz))
			send_vec=B(:,Nx-2:Nx,1:Ny,1:Nz)
			tag=41
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (2)
		allocate(send_sc(2,Ny,Nz))
			send_sc=rho(1:2,1:Ny,1:Nz)
			tag=12
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:2,1:Ny,1:Nz)
			tag=22
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,2,Ny,Nz))
			send_vec=U(:,1:2,1:Ny,1:Nz)
			tag=32
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1:2,1:Ny,1:Nz)
			tag=42
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (3)
		allocate(send_sc(Nx,2,Nz))
			send_sc=rho(1:Nx,Ny-1:Ny,1:Nz)
			tag=13
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,Ny-1:Ny,1:Nz)
			tag=23
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx,2,Nz))
			send_vec=U(:,1:Nx,Ny-1:Ny,1:Nz)
			tag=33
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec(3,Nx,3,Nz))
			send_vec=B(:,1:Nx,Ny-2:Ny,1:Nz)
			tag=43
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (4)
		allocate(send_sc(Nx,2,Nz))
			send_sc=rho(1:Nx,1:2,1:Nz)
			tag=14
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,1:2,1:Nz)
			tag=24
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx,2,Nz))
			send_vec=U(:,1:Nx,1:2,1:Nz)
			tag=34
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1:Nx,1:2,1:Nz)
			tag=44
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (5)
		allocate(send_sc(Nx,Ny,2))
			send_sc=rho(1:Nx,1:Ny,Nz-1:Nz)
			tag=15
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,1:Ny,Nz-1:Nz)
			tag=25
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx,Ny,2))
			send_vec=U(:,1:Nx,1:Ny,Nz-1:Nz)
			tag=35
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec(3,Nx,Ny,3))
			send_vec=B(:,1:Nx,1:Ny,Nz-2:Nz)
			tag=45
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (6)
		allocate(send_sc(Nx,Ny,2))
			send_sc=rho(1:Nx,1:Ny,1:2)
			tag=16
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,1:Ny,1:2)
			tag=26
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx,Ny,2))
			send_vec=U(:,1:Nx,1:Ny,1:2)
			tag=36
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1:Nx,1:Ny,1:2)
			tag=46
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
end select
enddo
end subroutine





subroutine receive_from(proc,rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: proc,communicator,tag,j
integer :: status(MPI_STATUS_SIZE)
integer :: i,proc_is_neighbour,neighbour_list(6)
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision, allocatable, dimension (:,:,:) :: recv_sc
double precision, allocatable, dimension (:,:,:,:) :: recv_vec


do i=1,6

neighbour_list = (/2,4,6,1,3,5/)
j=neighbour_list(i)
proc_is_neighbour=0
if (neighbours(iproc+1,j).EQ.(proc+1)) proc_is_neighbour = j

!if(proc_is_neighbour.ne.0) write (*,*) "This is",iproc,"receiving from",proc,"as neighbour",proc_is_neighbour
!if(proc_is_neighbour.ne.0) write (*,*) "This is",iproc,"Neighbours:",neighbours(iproc+1,1:6)

! Consider the following cases for each root block:
! from neighbour 1: you get the fields for their (1:2,1:Ny,1:Nz)  and make them your (Nx+1:Nx+2,1:Ny,1:Nz)
! from neighbour 2: you get the fields for their (Nx-1:Nx,1:Ny,1:Nz) and make them your (-1:0,1:Ny,1:Nz)
! from neighbour 3: you get the fields for their (1:Nx,1:2,1:Nz)  and make them your (1:Nx,Ny+1:Ny+2,1:Nz)
! from neighbour 4: you get the fields for their (1:Nx,Ny-1:Ny,1:Nz) and make them your (1:Nx,-1:0,1:Nz)
! from neighbour 5: you get the fields for their (1:Nx,1:Ny,1:2)  and make them your (1:Nx,1:Ny,Nz+1:Nz+2)
! from neighbour 6: you get the fields for their (1:Nx,1:Ny,Nz-1:Nz) and make them your (1:Nx,1:Ny,-1:0)
! For the B-field, there is one more at the "small" side

! tag system:
! tag = TN
! first digit: T: type - 1: rho, 2: energ, 3: U, 4: B
! second digit: N: neighbour number of the target


select case (proc_is_neighbour)
	case (1)
		allocate(recv_sc(2,Ny,Nz))
			tag=12
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1:Nx+2,1:Ny,1:Nz)=recv_sc(:,1:Ny,1:Nz)
			tag=22
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1:Nx+2,1:Ny,1:Nz)=recv_sc(:,1:Ny,1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,2,Ny,Nz))
			tag=32
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1:Nx+2,1:Ny,1:Nz)=recv_vec(:,:,1:Ny,1:Nz)
			tag=42
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1:Nx+2,1:Ny,1:Nz)=recv_vec(:,:,1:Ny,1:Nz)
		deallocate(recv_vec)
	case (2)
		allocate(recv_sc(2,Ny,Nz))
			tag=11
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(-1:0,1:Ny,1:Nz)=recv_sc(:,1:Ny,1:Nz)
			tag=21
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(-1:0,1:Ny,1:Nz)=recv_sc(:,1:Ny,1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,2,Ny,Nz))
			tag=31
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,-1:0,1:Ny,1:Nz)=recv_vec(:,:,1:Ny,1:Nz)
			deallocate(recv_vec)
			allocate(recv_vec(3,3,Ny,Nz))
			tag=41
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,-2:0,1:Ny,1:Nz)=recv_vec(:,:,1:Ny,1:Nz)
		deallocate(recv_vec)
	case (3)
		allocate(recv_sc(Nx,2,Nz))
			tag=14
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,Ny+1:Ny+2,1:Nz)=recv_sc(1:Nx,:,1:Nz)
			tag=24
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,Ny+1:Ny+2,1:Nz)=recv_sc(1:Nx,:,1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx,2,Nz))
			tag=34
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,Ny+1:Ny+2,1:Nz)=recv_vec(:,1:Nx,:,1:Nz)
			tag=44
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,Ny+1:Ny+2,1:Nz)=recv_vec(:,1:Nx,:,1:Nz)
		deallocate(recv_vec)
	case (4)
		allocate(recv_sc(Nx,2,Nz))
			tag=13
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,-1:0,1:Nz)=recv_sc(1:Nx,:,1:Nz)
			tag=23
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,-1:0,1:Nz)=recv_sc(1:Nx,:,1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx,2,Nz))
			tag=33
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,-1:0,1:Nz)=recv_vec(:,1:Nx,:,1:Nz)
			deallocate(recv_vec)
			allocate(recv_vec(3,Nx,3,Nz))
			tag=43
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,-2:0,1:Nz)=recv_vec(:,1:Nx,:,1:Nz)
		deallocate(recv_vec)
	case (5)
		allocate(recv_sc(Nx,Ny,2))
			tag=16
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,1:Ny,Nz+1:Nz+2)=recv_sc(1:Nx,1:Ny,:)
			tag=26
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,1:Ny,Nz+1:Nz+2)=recv_sc(1:Nx,1:Ny,:)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx,Ny,2))
			tag=36
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,1:Ny,Nz+1:Nz+2)=recv_vec(:,1:Nx,1:Ny,:)
			tag=46
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,1:Ny,Nz+1:Nz+2)=recv_vec(:,1:Nx,1:Ny,:)
		deallocate(recv_vec)
	case (6)
		allocate(recv_sc(Nx,Ny,2))
			tag=15
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,1:Ny,-1:0)=recv_sc(1:Nx,1:Ny,:)
			tag=25
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,1:Ny,-1:0)=recv_sc(1:Nx,1:Ny,:)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx,Ny,2))
			tag=35
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,1:Ny,-1:0)=recv_vec(:,1:Nx,1:Ny,:)
			deallocate(recv_vec)
			allocate(recv_vec(3,Nx,Ny,3))
			tag=45
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,1:Ny,-2:0)=recv_vec(:,1:Nx,1:Ny,:)
		deallocate(recv_vec)
end select

enddo

end subroutine














subroutine send_edge_to(proc,rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: proc,communicator,i,proc_is_neighbour,neighbour_list(12),tag,j
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision, allocatable, dimension (:) :: send_sc
double precision, allocatable, dimension (:,:) :: send_vec
double precision, allocatable, dimension (:,:,:) :: send_vec2
double precision, allocatable, dimension (:,:,:,:) :: send_vec3

integer :: RBi,RBj,RBk,RBin,RBjn,RBkn

do i=1,12
neighbour_list = (/1,2,3,4,5,6,7,8,9,10,11,12/) ! The order in which the neighbours are handled
j=neighbour_list(i)
proc_is_neighbour=0
if (edge_neighbours(iproc+1,j).EQ.(proc+1)) proc_is_neighbour = j

!call RBcoord(iproc+1,RBi,RBj,RBk)
!call RBcoord(proc+1,RBin,RBjn,RBkn)
!if((iproc.EQ.13).AND.(proc_is_neighbour.ne.0)) write (*,*) "This is",RBi,RBj,RBk,"sending to",RBin,RBjn,RBkn,"as edge neighbour",proc_is_neighbour



! Consider the following cases for each root block (for the edges):
! to neighbour  1: you send the fields for (Nx,1,1:Nz)     and make them their (0,Ny+1,1:Nz)
! to neighbour  2: you send the fields for (Nx,Ny,1:Nz)    and make them their (0,0,1:Nz)
! to neighbour  3: you send the fields for (1,Ny,1:Nz)     and make them their (Nx+1,0,1:Nz)
! to neighbour  4: you send the fields for (1,1,1:Nz)      and make them their (Nx+1,Ny+1,1:Nz)
! to neighbour  5: you send the fields for (Nx,1:Ny,Nz)    and make them their (0,1:Ny,0)
! to neighbour  6: you send the fields for (1:Nx,Ny,Nz)    and make them their (1:Nx,0,0)
! to neighbour  7: you send the fields for (1,1:Ny,Nz)     and make them their (Nx+1,1:Ny,0)
! to neighbour  8: you send the fields for (1:Nx,1,Nz)     and make them their (1:Nx,Ny+1,0)
! to neighbour  9: you send the fields for (Nx,1:Ny,1)     and make them their (0,1:Ny,Nz+1)
! to neighbour 10: you send the fields for (1:Nx,Ny,1)     and make them their (1:Nx,0,Nz+1)
! to neighbour 11: you send the fields for (1,1:Ny,1)      and make them their (Nx+1,1:Ny,Nz+1)
! to neighbour 12: you send the fields for (1:Nx,1,1)      and make them their (1:Nx,Ny+1,Nz+1)
! For the B-field, there is one more at the "small" side



! tag system:
! tag = TN
! first digit: T: type - 1: rho, 2: energ, 3: U, 4: B
! second digit: N: neighbour number of the target

select case (proc_is_neighbour)
	case (1)
		allocate(send_sc(Nz))
			send_sc=rho(Nx,1,1:Nz)
			tag=11
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,1,1:Nz)
			tag=21
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nz))
			send_vec=U(:,Nx,1,1:Nz)
			tag=31
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec2(3,2,Nz))
			send_vec2=B(:,Nx-1:Nx,1,1:Nz)
			tag=41
			call mpi_send(send_vec2,size(send_vec2),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec2)
	case (2)
		allocate(send_sc(Nz))
			send_sc=rho(Nx,Ny,1:Nz)
			tag=12
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,Ny,1:Nz)
			tag=22
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nz))
			send_vec=U(:,Nx,Ny,1:Nz)
			tag=32
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec3(3,2,2,Nz))
			send_vec3=B(:,Nx-1:Nx,Ny-1:Ny,1:Nz)
			tag=42
			call mpi_send(send_vec3,size(send_vec3),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec3)
	case (3)
		allocate(send_sc(Nz))
			send_sc=rho(1,Ny,1:Nz)
			tag=13
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,Ny,1:Nz)
			tag=23
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nz))
			send_vec=U(:,1,Ny,1:Nz)
			tag=33
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec2(3,2,Nz))
			send_vec2=B(:,1,Ny-1:Ny,1:Nz)
			tag=43
			call mpi_send(send_vec2,size(send_vec2),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec2)
	case (4)
		allocate(send_sc(Nz))
			send_sc=rho(1,1,1:Nz)
			tag=14
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,1,1:Nz)
			tag=24
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nz))
			send_vec=U(:,1,1,1:Nz)
			tag=34
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1,1,1:Nz)
			tag=44
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (5)
		allocate(send_sc(Ny))
			send_sc=rho(Nx,1:Ny,Nz)
			tag=15
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,1:Ny,Nz)
			tag=25
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Ny))
			send_vec=U(:,Nx,1:Ny,Nz)
			tag=35
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec3(3,2,Ny,2))
			send_vec3=B(:,Nx-1:Nx,1:Ny,Nz-1:Nz)
			tag=45
			call mpi_send(send_vec3,size(send_vec3),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec3)
	case (6)
		allocate(send_sc(Nx))
			send_sc=rho(1:Nx,Ny,Nz)
			tag=16
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,Ny,Nz)
			tag=26
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx))
			send_vec=U(:,1:Nx,Ny,Nz)
			tag=36
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec3(3,Nx,2,2))
			send_vec3=B(:,1:Nx,Ny-1:Ny,Nz-1:Nz)
			tag=46
			call mpi_send(send_vec3,size(send_vec3),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec3)
	case (7)
		allocate(send_sc(Ny))
			send_sc=rho(1,1:Ny,Nz)
			tag=17
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,1:Ny,Nz)
			tag=27
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Ny))
			send_vec=U(:,1,1:Ny,Nz)
			tag=37
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec2(3,Ny,2))
			send_vec2=B(:,1,1:Ny,Nz-1:Nz)
			tag=47
			call mpi_send(send_vec2,size(send_vec2),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec2)
	case (8)
		allocate(send_sc(Nx))
			send_sc=rho(1:Nx,1,Nz)
			tag=18
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,1,Nz)
			tag=28
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx))
			send_vec=U(:,1:Nx,1,Nz)
			tag=38
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec2(3,Nx,2))
			send_vec2=B(:,1:Nx,1,Nz-1:Nz)
			tag=48
			call mpi_send(send_vec2,size(send_vec2),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec2)
	case (9)
		allocate(send_sc(Ny))
			send_sc=rho(Nx,1:Ny,1)
			tag=19
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,1:Ny,1)
			tag=29
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Ny))
			send_vec=U(:,Nx,1:Ny,1)
			tag=39
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec2(3,2,Ny))
			send_vec2=B(:,Nx-1:Nx,1:Ny,1)
			tag=49
			call mpi_send(send_vec2,size(send_vec2),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec2)
	case (10)
		allocate(send_sc(Nx))
			send_sc=rho(1:Nx,Ny,1)
			tag=110
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,Ny,1)
			tag=210
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx))
			send_vec=U(:,1:Nx,Ny,1)
			tag=310
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			deallocate(send_vec)
			allocate(send_vec2(3,Nx,2))
			send_vec2=B(:,1:Nx,Ny-1:Ny,1)
			tag=410
			call mpi_send(send_vec2,size(send_vec2),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec2)
	case (11)
		allocate(send_sc(Ny))
			send_sc=rho(1,1:Ny,1)
			tag=111
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,1:Ny,1)
			tag=211
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Ny))
			send_vec=U(:,1,1:Ny,1)
			tag=311
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1,1:Ny,1)
			tag=411
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
	case (12)
		allocate(send_sc(Nx))
			send_sc=rho(1:Nx,1,1)
			tag=112
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1:Nx,1,1)
			tag=212
			call mpi_send(send_sc,size(send_sc),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_sc)
		allocate(send_vec(3,Nx))
			send_vec=U(:,1:Nx,1,1)
			tag=312
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1:Nx,1,1)
			tag=412
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
		deallocate(send_vec)
end select
enddo
end subroutine





subroutine receive_edge_from(proc,rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: proc,communicator,tag,j
integer :: status(MPI_STATUS_SIZE)
integer :: i,proc_is_neighbour,neighbour_list(12)
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision, allocatable, dimension (:) :: recv_sc
double precision, allocatable, dimension (:,:) :: recv_vec
double precision, allocatable, dimension (:,:,:) :: recv_vec2
double precision, allocatable, dimension (:,:,:,:) :: recv_vec3

do i=1,12
neighbour_list = (/3,4,1,2,11,12,9,10,7,8,5,6/)
j=neighbour_list(i)
proc_is_neighbour=0
if (edge_neighbours(iproc+1,j).EQ.(proc+1)) proc_is_neighbour = j



! Consider the following cases for each root block (for the edges):
! from neighbour  1: you get the fields for (1,Ny,1:Nz)     and make them your (Nx+1,0,1:Nz)
! from neighbour  2: you get the fields for (1,1,1:Nz)      and make them your (Nx+1,Ny+1,1:Nz)
! from neighbour  3: you get the fields for (Nx,1,1:Nz)     and make them your (0,Ny+1,1:Nz)
! from neighbour  4: you get the fields for (Nx,Ny,1:Nz)    and make them your (0,0,1:Nz)
! from neighbour  5: you get the fields for (1,1:Ny,1)      and make them your (Nx+1,1:Ny,Nz+1)
! from neighbour  6: you get the fields for (1:Nx,1,1)      and make them your (1:Nx,Ny+1,Nz+1)
! from neighbour  7: you get the fields for (Nx,1:Ny,1)     and make them your (0,1:Ny,Nz+1)
! from neighbour  8: you get the fields for (1:Nx,Ny,1)     and make them your (1:Nx,0,Nz+1)
! from neighbour  9: you get the fields for (1,1:Ny,Nz)     and make them your (Nx+1,1:Ny,0)
! from neighbour 10: you get the fields for (1:Nx,1,Nz)     and make them your (1:Nx,Ny+1,0)
! from neighbour 11: you get the fields for (Nx,1:Ny,Nz)    and make them your (0,1:Ny,0)
! from neighbour 12: you get the fields for (1:Nx,Ny,Nz)    and make them your (1:Nx,0,0)
! For the B-field, there is one more at the "small" side


! tag system:
! tag = TN
! first digit: T: type - 1: rho, 2: energ, 3: U, 4: B
! second digit: N: neighbour number of the target


select case (proc_is_neighbour)
	case (1)
		allocate(recv_sc(Nz))
			tag=13
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,0,1:Nz)=recv_sc(1:Nz)
			tag=23
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,0,1:Nz)=recv_sc(1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nz))
			tag=33
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,0,1:Nz)=recv_vec(:,1:Nz)
			deallocate(recv_vec)
			allocate(recv_vec2(3,2,Nz))
			tag=43
			call mpi_recv(recv_vec2, size(recv_vec2), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,-1:0,1:Nz)=recv_vec2(:,:,1:Nz)
		deallocate(recv_vec2)
	case (2)
		allocate(recv_sc(Nz))
			tag=14
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,Ny+1,1:Nz)=recv_sc(1:Nz)
			tag=24
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,Ny+1,1:Nz)=recv_sc(1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nz))
			tag=34
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,Ny+1,1:Nz)=recv_vec(:,1:Nz)
			tag=44
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,Ny+1,1:Nz)=recv_vec(:,1:Nz)
		deallocate(recv_vec)
	case (3)
		allocate(recv_sc(Nz))
			tag=11
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,Ny+1,1:Nz)=recv_sc(1:Nz)
			tag=21
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,Ny+1,1:Nz)=recv_sc(1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nz))
			tag=31
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,Ny+1,1:Nz)=recv_vec(:,1:Nz)
			deallocate(recv_vec)
			allocate(recv_vec2(3,2,Nz))
			tag=41
			call mpi_recv(recv_vec2, size(recv_vec2), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,-1:0,Ny+1,1:Nz)=recv_vec2(:,:,1:Nz)
		deallocate(recv_vec2)
	case (4)
		allocate(recv_sc(Nz))
			tag=12
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,0,1:Nz)=recv_sc(1:Nz)
			tag=22
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,0,1:Nz)=recv_sc(1:Nz)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nz))
			tag=32
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,0,1:Nz)=recv_vec(:,1:Nz)
			deallocate(recv_vec)
			allocate(recv_vec3(3,2,2,Nz))
			tag=42
			call mpi_recv(recv_vec3, size(recv_vec3), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,-1:0,-1:0,1:Nz)=recv_vec3(:,:,:,1:Nz)
		deallocate(recv_vec3)
	case (5)
		allocate(recv_sc(Ny))
			tag=111
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,1:Ny,Nz+1)=recv_sc(1:Ny)
			tag=211
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,1:Ny,Nz+1)=recv_sc(1:Ny)
		deallocate(recv_sc)
		allocate(recv_vec(3,Ny))
			tag=311
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,1:Ny,Nz+1)=recv_vec(:,1:Ny)
			tag=411
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,1:Ny,Nz+1)=recv_vec(:,1:Ny)
		deallocate(recv_vec)
	case (6)
		allocate(recv_sc(Nx))
			tag=112
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,Ny+1,Nz+1)=recv_sc(1:Nx)
			tag=212
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,Ny+1,Nz+1)=recv_sc(1:Nx)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx))
			tag=312
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,Ny+1,Nz+1)=recv_vec(:,1:Nx)
			tag=412
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,Ny+1,Nz+1)=recv_vec(:,1:Nx)
		deallocate(recv_vec)
	case (7)
		allocate(recv_sc(Ny))
			tag=19
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,1:Ny,Nz+1)=recv_sc(1:Ny)
			tag=29
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,1:Ny,Nz+1)=recv_sc(1:Ny)
		deallocate(recv_sc)
		allocate(recv_vec(3,Ny))
			tag=39
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,1:Ny,Nz+1)=recv_vec(:,1:Ny)
			deallocate(recv_vec)
			allocate(recv_vec2(3,2,Ny))
			tag=49
			call mpi_recv(recv_vec2, size(recv_vec2), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,-1:0,1:Ny,Nz+1)=recv_vec2(:,:,1:Ny)
		deallocate(recv_vec2)
	case (8)
		allocate(recv_sc(Nx))
			tag=110
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,0,Nz+1)=recv_sc(1:Nx)
			tag=210
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,0,Nz+1)=recv_sc(1:Nx)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx))
			tag=310
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,0,Nz+1)=recv_vec(:,1:Nx)
			deallocate(recv_vec)
			allocate(recv_vec2(3,Nx,2))
			tag=410
			call mpi_recv(recv_vec2, size(recv_vec2), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,-1:0,Nz+1)=recv_vec2(:,1:Nx,:)
		deallocate(recv_vec2)
	case (9)
		allocate(recv_sc(Ny))
			tag=17
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,1:Ny,0)=recv_sc(1:Ny)
			tag=27
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,1:Ny,0)=recv_sc(1:Ny)
		deallocate(recv_sc)
		allocate(recv_vec(3,Ny))
			tag=37
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,1:Ny,0)=recv_vec(:,1:Ny)
			deallocate(recv_vec)
			allocate(recv_vec2(3,Ny,2))
			tag=47
			call mpi_recv(recv_vec2, size(recv_vec2), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,1:Ny,-1:0)=recv_vec2(:,1:Ny,:)
		deallocate(recv_vec2)
	case (10)
		allocate(recv_sc(Nx))
			tag=18
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,Ny+1,0)=recv_sc(1:Nx)
			tag=28
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,Ny+1,0)=recv_sc(1:Nx)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx))
			tag=38
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,Ny+1,0)=recv_vec(:,1:Nx)
			deallocate(recv_vec)
			allocate(recv_vec2(3,Nx,2))
			tag=48
			call mpi_recv(recv_vec2, size(recv_vec2), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,Ny+1,-1:0)=recv_vec2(:,1:Nx,:)
		deallocate(recv_vec2)
	case (11)
		allocate(recv_sc(Ny))
			tag=15
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,1:Ny,0)=recv_sc(1:Ny)
			tag=25
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,1:Ny,0)=recv_sc(1:Ny)
		deallocate(recv_sc)
		allocate(recv_vec(3,Ny))
			tag=35
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,1:Ny,0)=recv_vec(:,1:Ny)
			deallocate(recv_vec)
			allocate(recv_vec3(3,2,Ny,2))
			tag=45
			call mpi_recv(recv_vec3, size(recv_vec3), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,-1:0,1:Ny,-1:0)=recv_vec3(:,:,1:Ny,:)
		deallocate(recv_vec3)
	case (12)
		allocate(recv_sc(Nx))
			tag=16
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(1:Nx,0,0)=recv_sc(1:Nx)
			tag=26
			call mpi_recv(recv_sc, size(recv_sc), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(1:Nx,0,0)=recv_sc(1:Nx)
		deallocate(recv_sc)
		allocate(recv_vec(3,Nx))
			tag=36
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,1:Nx,0,0)=recv_vec(:,1:Nx)
			deallocate(recv_vec)
			allocate(recv_vec3(3,Nx,2,2))
			tag=46
			call mpi_recv(recv_vec3, size(recv_vec3), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,1:Nx,-1:0,-1:0)=recv_vec3(:,1:Nx,:,:)
		deallocate(recv_vec3)
end select

enddo

end subroutine








subroutine send_corner_to(proc,rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: proc,communicator,i,proc_is_neighbour,neighbour_list(8),tag,j
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: send_sc
double precision :: send_vec(3)



do i=1,8
neighbour_list = (/1,2,3,4,5,6,7,8/) ! The order in which the neighbours are handled
j=neighbour_list(i)
proc_is_neighbour=0
if (corner_neighbours(iproc+1,j).EQ.(proc+1)) proc_is_neighbour = j



! Consider the following cases for each root block (for the corners):
! to neighbour  1: you send the fields for (Nx,1,Nz)    and make them their (0,Ny+1,0)
! to neighbour  2: you send the fields for (Nx,Ny,Nz)   and make them their (0,0,0)
! to neighbour  3: you send the fields for (1,Ny,Nz)    and make them their (Nx+1,0,0)
! to neighbour  4: you send the fields for (1,1,Nz)     and make them their (Nx+1,Ny+1,0)
! to neighbour  5: you send the fields for (Nx,1,1)     and make them their (0,Ny+1,Nz+1)
! to neighbour  6: you send the fields for (Nx,Ny,1)    and make them their (0,0,Nz+1)
! to neighbour  7: you send the fields for (1,Ny,1)     and make them their (Nx+1,0,Nz+1)
! to neighbour  8: you send the fields for (1,1,1)      and make them their (Nx+1,Ny+1,Nz+1)




! tag system:
! tag = TN
! first digit: T: type - 1: rho, 2: energ, 3: U, 4: B
! second digit: N: neighbour number of the target

select case (proc_is_neighbour)
	case (1)
			send_sc=rho(Nx,1,Nz)
			tag=11
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,1,Nz)
			tag=21
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,Nx,1,Nz)
			tag=31
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,Nx,1,Nz)
			tag=41
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (2)
			send_sc=rho(Nx,Ny,Nz)
			tag=12
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,Ny,Nz)
			tag=22
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,Nx,Ny,Nz)
			tag=32
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,Nx,Ny,Nz)
			tag=42
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (3)
			send_sc=rho(1,Ny,Nz)
			tag=13
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,Ny,Nz)
			tag=23
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,1,Ny,Nz)
			tag=33
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1,Ny,Nz)
			tag=43
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (4)
			send_sc=rho(1,1,Nz)
			tag=14
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,1,Nz)
			tag=24
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,1,1,Nz)
			tag=34
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1,1,Nz)
			tag=44
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (5)
			send_sc=rho(Nx,1,1)
			tag=15
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,1,1)
			tag=25
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,Nx,1,1)
			tag=35
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,Nx,1,1)
			tag=45
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (6)
			send_sc=rho(Nx,Ny,1)
			tag=16
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(Nx,Ny,1)
			tag=26
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,Nx,Ny,1)
			tag=36
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,Nx,Ny,1)
			tag=46
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (7)
			send_sc=rho(1,1,1)
			tag=17
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,Ny,1)
			tag=27
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,1,Ny,1)
			tag=37
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1,Ny,1)
			tag=47
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
	case (8)
			send_sc=rho(1,1,1)
			tag=18
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_sc=energ(1,1,1)
			tag=28
			call mpi_send(send_sc,1,MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=U(:,1,1,1)
			tag=38
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
			send_vec=B(:,1,1,1)
			tag=48
			call mpi_send(send_vec,size(send_vec),MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,ierror)
end select
enddo
end subroutine





subroutine receive_corner_from(proc,rho,energ,U,B)
use globals
implicit none
include "mpif.h"      !Include MPI library variables

integer :: proc,communicator,tag,j
integer :: status(MPI_STATUS_SIZE)
integer :: i,proc_is_neighbour,neighbour_list(8)
double precision :: U(3,-1:Nx+2,-1:Ny+2,-1:Nz+2),B(3,-2:Nx+2,-2:Ny+2,-2:Nz+2)
double precision :: energ(-1:Nx+2,-1:Ny+2,-1:Nz+2),rho(-1:Nx+2,-1:Ny+2,-1:Nz+2)
double precision :: recv_sc
double precision :: recv_vec(3)


do i=1,8
neighbour_list = (/7,8,5,6,3,4,1,2/)
j=neighbour_list(i)
proc_is_neighbour=0
if (corner_neighbours(iproc+1,j).EQ.(proc+1)) proc_is_neighbour = j


! Consider the following cases for each root block (for the corners):
! from neighbour  1: you get the fields for (1,Ny,1)     and make them your (Nx+1,0,Nz+1)
! from neighbour  2: you get the fields for (1,1,1)      and make them your (Nx+1,Ny+1,Nz+1)
! from neighbour  3: you get the fields for (Nx,1,1)     and make them your (0,Ny+1,Nz+1)
! from neighbour  4: you get the fields for (Nx,Ny,1)    and make them your (0,0,Nz+1)
! from neighbour  5: you get the fields for (1,Ny,Nz)    and make them your (Nx+1,0,0)
! from neighbour  6: you get the fields for (1,1,Nz)     and make them your (Nx+1,Ny+1,0)
! from neighbour  7: you get the fields for (Nx,1,Nz)    and make them your (0,Ny+1,0)
! from neighbour  8: you get the fields for (Nx,Ny,Nz)   and make them your (0,0,0)




! tag system:
! tag = TN
! first digit: T: type - 1: rho, 2: energ, 3: U, 4: B
! second digit: N: neighbour number of the target


select case (proc_is_neighbour)
	case (1)
			tag=17
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,0,Nz+1)=recv_sc
			tag=27
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,0,Nz+1)=recv_sc
			tag=37
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,0,Nz+1)=recv_vec(:)
			tag=47
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,0,Nz+1)=recv_vec(:)
	case (2)
			tag=18
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,Ny+1,Nz+1)=recv_sc
			tag=28
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,Ny+1,Nz+1)=recv_sc
			tag=38
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,Ny+1,Nz+1)=recv_vec(:)
			tag=48
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,Ny+1,Nz+1)=recv_vec(:)
	case (3)
			tag=15
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,Ny+1,Nz+1)=recv_sc
			tag=25
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,Ny+1,Nz+1)=recv_sc
			tag=35
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,Ny+1,Nz+1)=recv_vec(:)
			tag=45
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,0,Ny+1,Nz+1)=recv_vec(:)
	case (4)
			tag=16
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,0,Nz+1)=recv_sc
			tag=26
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,0,Nz+1)=recv_sc
			tag=36
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,0,Nz+1)=recv_vec(:)
			tag=46
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,0,0,Nz+1)=recv_vec(:)
	case (5)
			tag=13
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,0,0)=recv_sc
			tag=23
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,0,0)=recv_sc
			tag=33
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,0,0)=recv_vec(:)
			tag=43
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,0,0)=recv_vec(:)
	case (6)
			tag=14
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(Nx+1,Ny+1,0)=recv_sc
			tag=24
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(Nx+1,Ny+1,0)=recv_sc
			tag=34
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,Nx+1,Ny+1,0)=recv_vec(:)
			tag=44
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,Nx+1,Ny+1,0)=recv_vec(:)
	case (7)
			tag=11
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,Ny+1,0)=recv_sc
			tag=21
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,Ny+1,0)=recv_sc
			tag=31
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,Ny+1,0)=recv_vec(:)
			tag=41
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,0,Ny+1,0)=recv_vec(:)
	case (8)
			tag=12
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			rho(0,0,0)=recv_sc
			tag=22
			call mpi_recv(recv_sc, 1, MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			energ(0,0,0)=recv_sc
			tag=32
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			U(:,0,0,0)=recv_vec(:)
			tag=42
			call mpi_recv(recv_vec, size(recv_vec), MPI_DOUBLE_PRECISION,proc,tag,MPI_COMM_WORLD,status,ierror)
			B(:,0,0,0)=recv_vec(:)
end select

enddo

end subroutine










! sort an array by its entries (except if the entry is -1):
subroutine sort_array(array,length)
use globals
implicit none
integer :: length,i,swap
integer :: array(length)
logical :: change

change=.TRUE.

do while (change)
change=.FALSE.
do i=1,length-1
 if (array(i+1).NE.-1) then
   if (array(i+1).LT.array(i)) then
   	change=.TRUE.
   	swap=array(i+1)
   	array(i+1)=array(i)
   	array(i)=swap
   endif
 endif
enddo
enddo

end subroutine




! determine the neighbours of block RBnum
subroutine determine_neighbours()
use globals
implicit none
integer :: RBnumrun,RBnum
integer :: i,j,k


do RBnumrun=1,RBx*RBy*RBz
	call RBcoord(RBnumrun,i,j,k)

	neighbours(RBnumrun,1)=RBnum(i+1,j,k)
	neighbours(RBnumrun,2)=RBnum(i-1,j,k)
	neighbours(RBnumrun,3)=RBnum(i,j+1,k)
	neighbours(RBnumrun,4)=RBnum(i,j-1,k)
	neighbours(RBnumrun,5)=RBnum(i,j,k+1)
	neighbours(RBnumrun,6)=RBnum(i,j,k-1)

enddo
end subroutine



! determine the edge neighbours of block RBnum
subroutine determine_edge_neighbours()
use globals
implicit none
integer :: RBnumrun,RBnum
integer :: i,j,k

do RBnumrun=1,RBx*RBy*RBz
	call RBcoord(RBnumrun,i,j,k)

	edge_neighbours(RBnumrun,1)=RBnum(i+1,j-1,k)
	edge_neighbours(RBnumrun,2)=RBnum(i+1,j+1,k)
	edge_neighbours(RBnumrun,3)=RBnum(i-1,j+1,k)
	edge_neighbours(RBnumrun,4)=RBnum(i-1,j-1,k)
	edge_neighbours(RBnumrun,5)=RBnum(i+1,j,k+1)
	edge_neighbours(RBnumrun,6)=RBnum(i,j+1,k+1)
	edge_neighbours(RBnumrun,7)=RBnum(i-1,j,k+1)
	edge_neighbours(RBnumrun,8)=RBnum(i,j-1,k+1)
	edge_neighbours(RBnumrun,9)=RBnum(i+1,j,k-1)
	edge_neighbours(RBnumrun,10)=RBnum(i,j+1,k-1)
	edge_neighbours(RBnumrun,11)=RBnum(i-1,j,k-1)
	edge_neighbours(RBnumrun,12)=RBnum(i,j-1,k-1)


enddo

end subroutine


subroutine determine_corner_neighbours()
use globals
implicit none
integer :: RBnumrun,RBnum
integer :: i,j,k

do RBnumrun=1,RBx*RBy*RBz
	call RBcoord(RBnumrun,i,j,k)

	corner_neighbours(RBnumrun,1)=RBnum(i+1,j-1,k+1)
	corner_neighbours(RBnumrun,2)=RBnum(i+1,j+1,k+1)
	corner_neighbours(RBnumrun,3)=RBnum(i-1,j+1,k+1)
	corner_neighbours(RBnumrun,4)=RBnum(i-1,j-1,k+1)
	corner_neighbours(RBnumrun,5)=RBnum(i+1,j-1,k-1)
	corner_neighbours(RBnumrun,6)=RBnum(i+1,j+1,k-1)
	corner_neighbours(RBnumrun,7)=RBnum(i-1,j+1,k-1)
	corner_neighbours(RBnumrun,8)=RBnum(i-1,j-1,k-1)


enddo

end subroutine



! retrieve the integer coordinates for root block RBnum
subroutine RBcoord(RBnum,i,j,k)
use globals
implicit none
integer :: i,j,k,RBnum

i=1+(RBnum-1)/(RBz*RBy)
j=1+(RBnum-1-(i-1)*RBz*RBy)/RBz
k=RBnum-(i-1)*RBz*RBy-(j-1)*RBz

end subroutine



! retrieve the integer root block number for root block with coordinates i,j,k
integer function RBnum(i,j,k)
use globals
implicit none

integer :: i,j,k
integer :: itest,jtest,ktest

itest=i
jtest=j
ktest=k

if (i.LT.1) itest=i+RBx
if (i.GT.RBx) itest=i-RBx

if (j.LT.1) jtest=j+RBy
if (j.GT.RBy) jtest=j-RBy

if (k.LT.1) ktest=k+RBz
if (k.GT.RBz) ktest=k-RBz

RBnum=ktest+(jtest-1)*RBz+(itest-1)*RBz*RBy



end function
