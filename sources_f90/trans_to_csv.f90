!--------------------------------------------------------------------------------------------
!                This program reads the binary files of fargo3d and convert                     
!                them to csv files with all the fields in, the file can be
!                        read with Paraview, and similar sofrware.
!              Created by Barragán O., Rendón F. & Álvarez R.  Nov-Dec-2014                          
!--------------------------------------------------------------------------------------------
!                     Now it can be run in parallel using OPENMP (Jan 2015)
!--------------------------------------------------------------------------------------------

program TRANS_TO_CSV

  implicit none
  !INCLUDE THE MPI LIBRARY
  include "mpif.h"

!START VARIABLE DECLARATION
  !MPI Variables
  integer :: ierr, rank, nprocs
  integer status(MPI_STATUS_SIZE)
  !Tags for the communicators
  integer :: tlims, tptype, tind, tnxyz, titermin, titertot, titerjump
  !Local Variables
  integer :: i,j,k
  !prepare a vector to allocate the memory fields of FARGO3D
  !X ALWAYS is azimuth
  !Y corresponds to the radial component (in cylindrical or spherical)
  !Z to z in cylindrical and colatitude in spherical
  real*8, allocatable  :: dendata(:), vxdata(:), vydata(:), vzdata(:), endata(:)
  character(len=30) :: gasfield
  character(len=30) :: filenamebinary
  character(len=30) :: filename_csv 
  integer :: itermin, itertot, iterjump
  real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: limits(6) ! xmin, xmax, ymin, ymax, zmin, zmax
  character(30) :: table
  integer :: iter, NXYZ(3), NX, NY, NZ
  character(1) :: ptype

  !START THE MPI CALL
  call MPI_INIT(ierr)

  !Obtain the number of processors and the actual rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(iterjump * MPI_COMM_WORLD, rank, ierr)

  !Let the first processor do all this work
  if ( rank == 0 ) then

    !Create new directory csv_files
    call system('mkdir csv_files')

    !Read the data from the input.dat file
    call read_input(NXYZ,itermin,itertot,iterjump,limits,ptype)

  end if

  !First processor, tell to the other processors what you read
  if (nprocs > 1 ) then

    tlims     = 1
    tptype    = 2
    tnxyz     = 4
    titertot  = 5
    titermin  = 6
    titerjump = 7

    if ( rank > 0 ) then
      call MPI_RECV(limits  ,6,MPI_REAL8    ,0,tlims    ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(ptype   ,1,MPI_CHARACTER,0,tptype   ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(NXYZ    ,3,MPI_INTEGER  ,0,tnxyz    ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(itertot ,1,MPI_INTEGER  ,0,titertot ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(itermin ,1,MPI_INTEGER  ,0,titermin ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(iterjump,1,MPI_INTEGER  ,0,titerjump,MPI_COMM_WORLD,status,ierr) 
    else if ( rank == 0 ) then
      do i = 1, nprocs - 1 
        call MPI_SEND(limits  ,6,MPI_REAL8    ,i,tlims    ,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(ptype   ,1,MPI_CHARACTER,i,tptype   ,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(NXYZ    ,3,MPI_INTEGER  ,i,tnxyz    ,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(itertot ,1,MPI_INTEGER  ,i,titertot ,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(itermin ,1,MPI_INTEGER  ,i,titermin ,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(iterjump,1,MPI_INTEGER  ,i,titerjump,MPI_COMM_WORLD,status,ierr) 
      end do
    end if

  end if
  !Now all the processors know what is going on 

  !Put the values in an easier readable variables
  NX = NXYZ(1)
  NY = NXYZ(2)
  NZ = NXYZ(3)

  !Let's allocate the memory to read the binary data
  allocate( dendata(NX*NY*NZ) )
  allocate( vxdata(NX*NY*NZ)  )
  allocate( vydata(NX*NY*NZ)  )
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
  iter = itermin + iterjump * rank 
  !Start the work now!
  do while ( iter <= itertot )  

     print*,'Creating disk',iter,'.csv file by processor ',rank + 1

     !create denstiy data
     gasfield = 'gasdens'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading density file
        read(100,rec=1) dendata
     close(100)

     !create vx data
     gasfield = 'gasvx'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening vx binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading vx file
        read(100,rec=1) vxdata
     close(100)

     !create vy data
     gasfield = 'gasvy'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening vy binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading vy file
        read(100,rec=1) vydata
     close(100)

     !create vz data
     gasfield = 'gasvz'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening vz binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading vz file
        read(100,rec=1) vzdata
     close(100)

     !create energy data
     gasfield = 'gasenergy'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening energy binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", access="direct", recl = NX*NY*NZ*8) 
        !Reading energy file
        read(100,rec=1) endata
     close(100)

     !The binary data has been stored in the arrays
 
     !Create the name of the csv filename and it is stored in table
     call create_csv_files(iter,filename_csv)
     table='csv_files/'//trim(filename_csv)
   
     if ( ptype == 'c' ) then !cylindrical
       call cyl2carte(NX,NY,NZ,dendata,vydata,vxdata,vzdata,endata,ymin,ymax,xmin,xmax,zmin,zmax,table )
     else if ( ptype == 's' ) then !spherical
       call sph2carte(NX,NY,NZ,dendata,vydata,vxdata,vzdata,endata,ymin,ymax,xmin,xmax,zmin,zmax,table )
     else if ( ptype == 'r' ) then !cartesian
       call carte2carte(NX,NY,NZ,dendata,vydata,vxdata,vzdata,endata,ymin,ymax,xmin,xmax,zmin,zmax,table )
     else
       print*,'ERROR! your ptype is not allowed'
       stop
     end if

     !Make a jump as big as the number of processors
     iter = iter + iterjump * nprocs

  enddo

  !deallocate the memory
  deallocate(dendata)
  deallocate(vxdata)
  deallocate(vydata)
  deallocate(vzdata)
  deallocate(endata)

  print*,'The csv files have been created inside csv_files directory'

  !Wait until all the processors finish
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  !END MPI CALL
  call MPI_FINALIZE(ierr)

end program TRANS_TO_CSV
