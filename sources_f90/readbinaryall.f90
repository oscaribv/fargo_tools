! ----------------------------------------------------------------
! This program reads the binary files of fargo3d and convert    --
! them to ascii files, to generate automatically png files and  --
! mpg videos.                                                   --
! Created by Barragán O. and Rendón F. since Oct-15-2014        --
! ----------------------------------------------------------------

program FARGO3D_READER

  implicit none
  !INCLUDE THE MPI LIBRARY
  include "mpif.h"

  !MPI Variables
  integer :: ierr, rank, nprocs
  integer status(MPI_STATUS_SIZE)
  !Tags for the communicators
  integer :: tlims, tptype, tind, tnxyz, titertot, tmydirsetup

  !Local Variables
  integer :: i,j,k
  !prepare a vector to allocate the memory
  real*8, allocatable  :: bindata(:)
  character(len=15) :: filenamebinary
  character(len=99) :: myfilenamebinary
  character(len=15) :: filenametxt
  integer :: itertot
  real*8 :: xmin, xmax, ymin, ymax, zmin, zmax
  real*8 :: limits(6) ! xmin, xmax, ymin, ymax, zmin, zmax
  character(100) :: table
  integer :: iter, ind, NXYZ(3), NX, NY, NZ
  character(1) :: ptype
  character(1),parameter :: symbol='/'
  character(len=99) :: mydirsetup

  !symbol = '/'

  !START THE MPI CALL
  call MPI_INIT(ierr)

  !Obtain the number of processors and the actual rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  !Let the first processor do all this work
  if ( rank == 0 ) then

  !Create new directories cartesiantxt and images
  call system('rm -r cartesiantxt')
  call system('mkdir cartesiantxt')
  call system('rm -r images')
  call system('mkdir images')

  !Read the data form the input.dat file
  call read_input(NXYZ,itertot,limits,ptype,ind,mydirsetup)
  !Create the gnuplot file based on the input
  call create_gpl_file(NX,NY,NZ,limits,ptype,ind)

  end if

  !First processor, tell to the other processors what you read
  if (nprocs > 1 ) then
      tlims     = 1
      tptype    = 2
      tind      = 3
      tnxyz     = 4
      titertot  = 5
      tmydirsetup = 6
      if ( rank > 0 ) then
        call MPI_RECV(limits,6,MPI_REAL8,0,tlims,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(ptype,1,MPI_CHARACTER,0,tptype,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(ind,1,MPI_INTEGER,0,tind,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(NXYZ,3,MPI_INTEGER,0,tnxyz,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(itertot,1,MPI_INTEGER,0,titertot,MPI_COMM_WORLD,status,ierr) 
        call MPI_RECV(mydirsetup,99,MPI_CHARACTER,0,tmydirsetup,MPI_COMM_WORLD,status,ierr) 
      else if ( rank == 0 ) then
        do i = 1, nprocs - 1 
          call MPI_SEND(limits,6,MPI_REAL8,i,tlims,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(ptype,1,MPI_CHARACTER,i,tptype,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(ind,1,MPI_INTEGER,i,tind,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(NXYZ,3,MPI_INTEGER,i,tnxyz,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(itertot,1,MPI_INTEGER,i,titertot,MPI_COMM_WORLD,status,ierr) 
          call MPI_SEND(mydirsetup,99,MPI_CHARACTER,i,tmydirsetup,MPI_COMM_WORLD,status,ierr) 
        end do
      end if

  end if
  !Now all the processors know what is going on 

  !Put the values of the dimenssion in easier readable variables
  NX = NXYZ(1)
  NY = NXYZ(2)
  NZ = NXYZ(3)

  !Let's allocate the memory for read the binary data
  allocate(bindata(NX*NY*NZ))

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

  call createdensitybinaryfilename(iter,filenamebinary)
  !myfilenamebinary = '../outputs/mri_den/'//trim(filenamebinary)
!  myfilenamebinary = '../outputs/'//trim(mydirsetup2)//symbol//trim(filenamebinary)
!  print*, 'mydirsetup: ', mydirsetup
  myfilenamebinary = '../outputs/'//trim(mydirsetup)//symbol//trim(filenamebinary)
  
  !Opening binary file
  open(unit=100,                 &
       status="old",             &
       file=myfilenamebinary,      &
       form="unformatted",       &
       access="direct",          &
       recl = NX*NY*NZ*8)

  !Reading file
  read(100,rec=1) bindata
  close(100)
  !The binary data has been saved in bindata
 
  !Create the name of the txt filename and it is stored in table
  call createdensitytxtfilename(iter,filenametxt)
  table = 'cartesiantxt/'//trim(filenametxt)


  !The type of plot is given in the input file
  !The index also is given in the input file
  if ( ptype == 'z' ) then
    call cut_at_Zk(NX,NY,NZ,ind,bindata,ymin,ymax,xmin,xmax,zmin,zmax,table)
  else if ( ptype == 't' ) then
    call cut_at_Ti(NX,NY,NZ,ind,bindata,ymin,ymax,xmin,xmax,zmin,zmax,table)
  else
    print*,'I do not know what you want to do, bye'
    stop
  end if

  !Make a jump as big as the number of processors
  iter = iter + nprocs

  enddo

  !deallocate the memory
  deallocate(bindata)

  !Wait until all the processors finish
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  !The rest of the work is done by just one processor
  if ( rank == 0 ) then 
 
    print*,'Creating png files'
    call system('ls cartesiantxt/*.txt > cartesiantxt/dummy.dat')
    call system('cut -c 21-24 ./cartesiantxt/dummy.dat  > cartesiantxt/enter')
    call system('./makecartfigs.sh')
  
    print*,'Creating video'
    call system('ffmpeg -f image2 -i images/cartesiandisk%04d.png -qscale 1 images/cartesiandisk.mpg')
    call system('totem images/cartesiandisk.mpg')

  end if

  !END MPI CALL
  call MPI_FINALIZE(ierr)

end program FARGO3D_READER
