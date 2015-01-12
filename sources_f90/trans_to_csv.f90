! ----------------------------------------------------------------
!   This program reads the binary files of fargo3d and convert  
!    them to csv files with all the fields in, the file can be
!           read with Paraview, and others similars.
!             Created by Barragán O.  Nov-Dec-2014                          
! ----------------------------------------------------------------

program TRANS_TO_CSV

  implicit none
  !INCLUDE THE MPI LIBRARY
  include "mpif.h"

  !MPI Variables
  integer :: ierr, rank, nprocs
  integer status(MPI_STATUS_SIZE)
  !Tags for the communicators
  integer :: tlims, tptype, tind, tnxyz, titertot

  !Local Variables
  integer :: i,j,k
  !prepare a vector to allocate the memory
  real*8, allocatable  :: dendata(:), vrdata(:), vtdata(:), vzdata(:), endata(:)
  character(len=30) :: filenamebinary
  character(len=30) :: filename_csv, mydirsetup
  integer :: itertot
  real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: limits(6) ! xmin, xmax, ymin, ymax, zmin, zmax
  character(30) :: table
  integer :: iter, ind, NXYZ(3), NX, NY, NZ
  character(1) :: ptype

  !START THE MPI CALL
  call MPI_INIT(ierr)

  !Obtain the number of processors and the actual rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  !Let the first processor do all this work
  if ( rank == 0 ) then

    !Create new directories cartesiantxt and images
    call system('mkdir csv_files')

    !Read the data form the input.dat file
    call read_input(NXYZ,itertot,limits,ptype,ind,mydirsetup)

  end if

  !First processor, tell to the other processors what you read
  if (nprocs > 1 ) then

    tlims     = 1
    tptype    = 2
    tind      = 3
    tnxyz     = 4
    titertot  = 5

    if ( rank > 0 ) then
      call MPI_RECV(limits,6,MPI_REAL8,0,tlims,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(ptype,1,MPI_CHARACTER,0,tptype,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(ind,1,MPI_INTEGER,0,tind,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(NXYZ,3,MPI_INTEGER,0,tnxyz,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(itertot,1,MPI_INTEGER,0,titertot,MPI_COMM_WORLD,status,ierr) 
    else if ( rank == 0 ) then
      do i = 1, nprocs - 1 
        call MPI_SEND(limits,6,MPI_REAL8,i,tlims,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(ptype,1,MPI_CHARACTER,i,tptype,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(ind,1,MPI_INTEGER,i,tind,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(NXYZ,3,MPI_INTEGER,i,tnxyz,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(itertot,1,MPI_INTEGER,i,titertot,MPI_COMM_WORLD,status,ierr) 
      end do
    end if

  end if
  !Now all the processors know what is going on 

  !Put the values in an easier readable variables
  NX = NXYZ(1)
  NY = NXYZ(2)
  NZ = NXYZ(3)

  !Let's allocate the memory for read the binary data
  allocate( dendata(NX*NY*NZ) )
  allocate( vrdata(NX*NY*NZ)  )
  allocate( vtdata(NX*NY*NZ)  )
  allocate( vzdata(NX*NY*NZ)  )
  allocate( endata(NX*NY*NZ)  )

  !put the values in x,y,z to avoid confussion
  xmin = limits(1)
  xmax = limits(2)
  ymin = limits(3)
  ymax = limits(4)
  zmin = limits(5)
  zmax = limits(6)

  !Let's divide all the work between all the processors
  iter = rank + 1
  !Start the work now!
  do while ( iter <= itertot )  

     print*,'iter = ',iter,' in processor ',rank + 1

     !create denstiy data
    !call createdensitybinaryfilename(iter,filenamebinary)
     call createfieldbinaryfilename(iter,'gasdens',filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading density file
        read(100,rec=1) dendata
     close(100)

     !create vx (vr) data
     !call createvybinaryfilename(iter,filenamebinary)
     call createfieldbinaryfilename(iter,'gasvx',filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading density file
        read(100,rec=1) vrdata
     close(100)

     !create vy data
     !call createvxbinaryfilename(iter,filenamebinary)
     call createfieldbinaryfilename(iter,'gasvy',filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading density file
        read(100,rec=1) vtdata
     close(100)

     !create vz data
     !call createvzbinaryfilename(iter,filenamebinary)
     call createfieldbinaryfilename(iter,'gasvz',filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading density file
        read(100,rec=1) vzdata
     close(100)

     !create energy data
     !call createenergybinaryfilename(iter,filenamebinary)
     call createfieldbinaryfilename(iter,'gasenergy',filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading density file
        read(100,rec=1) endata
     close(100)

     !The binary data has been saved in all the arrays
 
     !Create the name of the txt filename and it is stored in table
     call create_csv_files(iter,filename_csv)
     table='csv_files/'//trim(filename_csv)

!-------------------------------------------------------------------------------
! Modify this to be used in all the coordinates systems
!-------------------------------------------------------------------------------


     !The type of plot is given in the input file
     !The index also is given in the input file
     !call cyl2carte(NX,NY,NZ,dendata,vrdata,vtdata,vzdata,endata,ymin,ymax,xmin,xmax,zmin,zmax,table )
     call sph2carte(NX,NY,NZ,dendata,vrdata,vtdata,vzdata,endata,ymin,ymax,xmin,xmax,zmin,zmax,table )

!-------------------------------------------------------------------------------

     !Make a jump as big as the number of processors
     iter = iter + nprocs

  enddo

  !deallocate the memory
  deallocate(dendata)
  deallocate(vrdata)
  deallocate(vtdata)
  deallocate(vzdata)
  deallocate(endata)

  !Wait until all the processors finish
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  !END MPI CALL
  call MPI_FINALIZE(ierr)

end program TRANS_TO_CSV
