!-------------------------------------------------------------------------------
!This file reads the binary data and creates an azimuthal averaged density of  -
!the disk                                                                      -
!Taken from readbinaryall.f90 and modified by Oscar Barragán, Oct-22-2014      -
!-------------------------------------------------------------------------------

!Parameters:
#define NX 1800
#define NY 600
#define NZ 1

program AVERAGE_DENSITY

  implicit none
  !INCLUDE THE MPI LIBRARY
  include "mpif.h"

  !MPI Variables
  integer :: ierr, rank, nprocs
  integer status(MPI_STATUS_SIZE)
  !Tags for the communicators
  integer :: tlims, tptype, tind

  !Local Variables
  integer :: i,j,k
  real*8  :: data(NX*NY*NZ)
  character(len=15) :: filenamebinary
  character(len=15) :: filenametxt
  integer, parameter :: itertot=1
  real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: limits(6) ! xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: avden, r, dr
  real*8 :: arup, ardn, area !area up and area down, area
  character(100) :: table
  integer :: iter, ind
  character(1) :: ptype

  !START THE MPI CALL
  call MPI_INIT(ierr)

  !Obtain the number of processors and the actual rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  !Let the first processor do all this work
  if ( rank == 0 ) then

  !It uses the same input file that readbinaryall
  !Read the data form the input.dat file
  call read_input(limits,ptype,ind)

  end if

  !First processor, tell to the other processors what you read
  if (nprocs > 1 ) then
      tlims = 1
      tptype = 2
      tind = 3
      if ( rank > 0 ) then
        call MPI_RECV(limits,6,MPI_REAL8,0,tlims,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(ptype,1,MPI_CHARACTER,0,tptype,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(ind,1,MPI_INTEGER,0,tind,MPI_COMM_WORLD,status,ierr) 
      else if ( rank == 0 ) then
        do i = 1, nprocs - 1 
          call MPI_SEND(limits,6,MPI_REAL8,i,tlims,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(ptype,1,MPI_CHARACTER,i,tptype,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(ind,1,MPI_INTEGER,i,tind,MPI_COMM_WORLD,status,ierr) 
        end do
      end if

  end if
  !Now all the processors know what is going on 

  !put the values in x,y,z to avoid confussion
  xmin = limits(1); xmax = limits(2)
  ymin = limits(3); ymax = limits(4)
  zmin = limits(5); zmax = limits(6)

  !Let's divide all the work between all the processors
  iter = rank + 1
  do while ( iter <= itertot )  

  print*,'iter = ',iter,' in processor ',rank + 1

  call createdensitybinaryfilename(iter,filenamebinary)
  
  !Opening binary file
  open(unit=100,                 &
       status="old",             &
       file=filenamebinary,      &
       form="unformatted",       &
       access="direct",          &
       recl = NX*NY*NZ*8)

  !Reading file
  read(100,rec=1) data
  close(100)


  open(101,file='avden.dat')

  !Lets calculate the azimuthal averaged density from data

 dr = ( ymax - ymin ) / ( NY - 1 ) 
 r = ymin
 !Let's fill the file
   do j=1, NY
     r = r + dr
     avden = 0.0
     do k=1, NZ
        do i=1, NX
          avden = avden + data(i+(j-1)*NX+(k-1)*NX*NY)
        end do
     end do
          !arup = r*r
          !ardn = ( r -dr )*( r - dr )
          !area = arup - ardn
          !now avden is the surface density, divide by area to get mass
          !avden = avden / area !Now this is the mass
          write(101,*) r, avden
   end do

  close(101)

  iter = iter + nprocs

  end do ! end iter do

  !Wait until all the processors finish
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  !END MPI CALL
  call MPI_FINALIZE(ierr)

end program AVERAGE_DENSITY
