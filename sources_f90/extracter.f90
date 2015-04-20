!---------------------------------------------------------------------------------------!
!                                      extracter.f90
!                               Written by O. Barragán (Oct-21-2014)
!       It takes the original NX x NY x NZ binary files, then creates new binary
!       arrays with sizes (NXmax - MXmin) x (NYmax - NYmin) x (NZmax - MZmin),
!       these values are set in the input_ext.dat file. 
!       The new files can be transformed to csv files with t2csv, this allows a
!       faster analysis.
!---------------------------------------------------------------------------------------!

program EXTRACTER

!-----------------------------------------------------------------
!       Variables declaration                                   --
!-----------------------------------------------------------------

  implicit none
  !INCLUDE THE MPI LIBRARY
  include "mpif.h"

  !MPI Variables
  integer :: ierr, rank, nprocs
  integer status(MPI_STATUS_SIZE)
  !Tags for the communicators
  integer :: tlims, tenes, tnxyz, titertot, titermin, titerjump, tptype

  !Local Variables
  integer :: i,j,k
  !variables to read the original bin files
  real*8, allocatable  :: dendata(:), vxdata(:), vydata(:), vzdata(:), endata(:)
  !c for circumplanetary data bin files
  real*8, allocatable  :: cdendata(:), cvxdata(:), cvydata(:), cvzdata(:), cendata(:)
  !names for the binary files
  character(len=30) :: filenametxt, filenamebinary
  !number of iterations, same as number of set of files in FARGO3D
  integer :: itertot, itermin, iterjump
  real*8 :: limits(6) ! xmin, xmax, ymin, ymax, zmin, zmax
  character(100) :: table
  character(30) :: fileextracter, gasfield
  integer :: iter, ind
  character(1) :: ptype
  integer :: ENES(6), NXmin, NXmax, NYmin, NYmax, NZmin, NZmax
  integer :: NXYZ(3), NX, NY, NZ
  integer :: sizex, sizey, sizez
  real :: xmin, xmax, ymin, ymax, zmin, zmax, dx, dy, dz
  real*8, allocatable :: xdata(:)

  !START THE MPI CALL
  call MPI_INIT(ierr)

  !Obtain the number of processors and the actual rank
  call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

  !Let the first processor do all this work
  if ( rank == 0 ) then

   !Creating circumplanetary disk directory
   call system('mkdir extracter_files')

   !Read the data form the input_circumplanetary.dat file
   call read_input_extracter(NXYZ,itermin, itertot, iterjump, limits, ENES, ptype)

  end if

  !First processor, tell to the other processors what you read
  if (nprocs > 1 ) then
    tlims =     1
    tenes =     2
    tnxyz =     3
    titertot =  4
    titermin =  5
    titerjump = 6
    tptype =    7
    if ( rank > 0 ) then
      call MPI_RECV(limits,  6,MPI_REAL8,    0,tlims    ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(ENES,    6,MPI_INTEGER,  0,tenes    ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(NXYZ,    3,MPI_INTEGER,  0,tnxyz    ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(itertot, 1,MPI_INTEGER,  0,titertot ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(itermin, 1,MPI_INTEGER,  0,titermin ,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(iterjump,1,MPI_INTEGER,  0,titerjump,MPI_COMM_WORLD,status,ierr) 
      call MPI_RECV(ptype,   1,MPI_CHARACTER,0,tptype   ,MPI_COMM_WORLD,status,ierr) 
    else if ( rank == 0 ) then
      do i = 1, nprocs - 1 
        call MPI_SEND(limits,  6,MPI_REAL8,    i,tlims,    MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(ENES,    6,MPI_INTEGER,  i,tenes,    MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(NXYZ,    3,MPI_INTEGER,  i,tnxyz,    MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(itertot, 1,MPI_INTEGER,  i,titertot, MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(itermin, 1,MPI_INTEGER,  i,titermin, MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(iterjump,1,MPI_INTEGER,  i,titerjump,MPI_COMM_WORLD,status,ierr) 
        call MPI_SEND(ptype,   1,MPI_CHARACTER,i,tptype,   MPI_COMM_WORLD,status,ierr) 
      end do
    end if

  end if
  !Now all the processors know what is going on 

  !put the values in an easier readable way

  !Taking the dimenssions of the original binary arrray
  NX = NXYZ(1)
  NY = NXYZ(2)
  NZ = NXYZ(3)

  !These are the limits where the original array will be cut
  NXmin = ENES(1)
  NXmax = ENES(2)
  NYmin = ENES(3)
  NYmax = ENES(4)
  NZmin = ENES(5)
  NZmax = ENES(6)

  !These are the dimenssions of the new binary array
  sizex = NXmax - NXmin + 1
  sizey = NYmax - NYmin + 1
  sizez = NZmax - NZmin + 1

  !These are the jumps in each directions, just for equally espaced jumps!
  dx = ( limits(2) - limits(1) ) / (NX - 1)
  dy = ( limits(4) - limits(3) ) / (NY - 1)
  dz = ( limits(6) - limits(5) ) / (NZ - 1)

  !Get the new grid limits taking in account the position where we want to cut
  xmin = limits(1) + (NXmin - 1) * dx  
  xmax = limits(1) + (NXmax - 1) * dx  
  ymin = limits(3) + (NYmin - 1) * dy  
  ymax = limits(3) + (NYmax - 1) * dy  
  zmin = limits(5) + (NZmin - 1) * dz  
  zmax = limits(5) + (NZmax - 1) * dz  


  if ( rank == 0 ) then

    !create input.dat to run the output files of this program
    !with trans_to_csv.f90
    open(13,file='extracter_files/input.dat')
   
      write(13,*),sizex,'    !NX' 
      write(13,*),sizey,'    !NY' 
      write(13,*),sizez,'    !NZ' 
      write(13,*),itermin, '    !itermin'
      write(13,*),itertot, '    !itertot'
      write(13,*),iterjump, '    !iterjump'
      write(13,*),xmin, '    !X min' 
      write(13,*),xmax, '    !X max' 
      write(13,*),ymin, '    !R min' 
      write(13,*),ymax, '    !R max' 
      write(13,*),zmin, '    !Z min' 
      write(13,*),zmax, '    !Z max' 
      write(13,*),ptype,'    !type of plot' 

    close(13)
 
  end if  

  !allocate bindata that contains all the data of all the grid
  !Let's allocate the memory for read the binary data
  allocate( dendata(NX*NY*NZ) )
  allocate( vxdata (NX*NY*NZ) )
  allocate( vydata (NX*NY*NZ) )
  allocate( vzdata (NX*NY*NZ) )
  allocate( endata (NX*NY*NZ) )

  !allocate xdata  of the new Size
  !We call it out the do, sice it will recycled for each processor in each iter
  allocate( cdendata(sizex*sizey*sizez) ) 
  allocate( cvxdata (sizex*sizey*sizez) ) 
  allocate( cvydata (sizex*sizey*sizez) ) 
  allocate( cvzdata (sizex*sizey*sizez) ) 
  allocate( cendata (sizex*sizey*sizez) ) 


  !Let's divide all the work between all the processors
  iter = itermin + iterjump * rank
  !Start the work now!
  do while ( iter <= itertot )  

    print*,'Transforming output = ',iter,' in processor ',rank + 1

!----------------------------------------------------------------------
!       Reading the original binary files
!----------------------------------------------------------------------

     !create denstiy data
     gasfield = 'gasdens'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening density binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", &
     access="direct", recl = NX*NY*NZ*8)
        !Reading density file
        read(100,rec=1) dendata
     close(100)

     !create vx data
     gasfield = 'gasvx'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening vx binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", &
access="direct", recl = NX*NY*NZ*8)
        !Reading vx file
        read(100,rec=1) vxdata
     close(100)


     !create vy data
     gasfield = 'gasvy'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening vy binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", &
access="direct", recl = NX*NY*NZ*8)
        !Reading vy file
        read(100,rec=1) vydata
     close(100)

     !create vz data
     gasfield = 'gasvz'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening vz binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", &
access="direct", recl = NX*NY*NZ*8)
        !Reading vz file
        read(100,rec=1) vzdata
     close(100)

     !create energy data
     gasfield = 'gasenergy'
     call createfieldbinaryfilename(iter,gasfield,filenamebinary)
     !Opening energy binary file
     open(unit=100, status="old", file=filenamebinary, form="unformatted", &
access="direct", recl = NX*NY*NZ*8)
        !Reading energy file
        read(100,rec=1) endata
     close(100)

     !The binary data has been stored in the arrays

!----------------------------------------------------------------------
!       Now all the original binary files are stored in RAM
!----------------------------------------------------------------------
!       Let's take the data that we want, previously we decided this
!       when we choose the NXmin, Mxmax, etc.
!----------------------------------------------------------------------

    !Reading from *data variables to c*data variables
    do k = NZmin, NZmax  
      do j = NYmin, NYmax  
        do i = NXmin, NXmax  

          cdendata((i-NXmin)+(j-NYmin)*sizex+(k-NZmin)*sizex*sizey + 1) & 
          = dendata(i+(j-1)*NX+(k-1)*NX*NY)

          cvxdata((i-NXmin)+(j-NYmin)*sizex+(k-NZmin)*sizex*sizey + 1) & 
          = vxdata(i+(j-1)*NX+(k-1)*NX*NY)

          cvydata((i-NXmin)+(j-NYmin)*sizex+(k-NZmin)*sizex*sizey + 1) & 
          = vydata(i+(j-1)*NX+(k-1)*NX*NY)

          cvzdata((i-NXmin)+(j-NYmin)*sizex+(k-NZmin)*sizex*sizey + 1) & 
          = vzdata(i+(j-1)*NX+(k-1)*NX*NY)

          cendata((i-NXmin)+(j-NYmin)*sizex+(k-NZmin)*sizex*sizey + 1) & 
          = endata(i+(j-1)*NX+(k-1)*NX*NY)

        end do
      end do
    end do

!----------------------------------------------------------------------
!       Now all the data that we want is stored in RAM, lets write it
!       in the hard disk, the names are as the original ones but
!       inside extracter_data directory
!----------------------------------------------------------------------

      !create denstiy data
  gasfield = 'gasdens'
  call createfieldbinaryfilename(iter,gasfield,fileextracter)
  table = 'extracter_files/'//trim(fileextracter)
  open(unit=23, status='unknown',file=table,form="unformatted", &
       access="direct", recl= sizex*sizey*sizez*8 )  
  write(23,rec=1) cdendata
  !Now the circumplanetary data is stored in a binary file
  close(23)

  gasfield = 'gasvx'
  call createfieldbinaryfilename(iter,gasfield,fileextracter)
  table = 'extracter_files/'//trim(fileextracter)
  open(unit=23, status='unknown',file=table,form="unformatted", &
       access="direct", recl= sizex*sizey*sizez*8 )  
  write(23,rec=1) cvxdata
  !Now the circumplanetary data is stored in a binary file
  close(23)

  gasfield = 'gasvy'
  call createfieldbinaryfilename(iter,gasfield,fileextracter)
  table = 'extracter_files/'//trim(fileextracter)
  open(unit=23, status='unknown',file=table,form="unformatted", &
       access="direct", recl= sizex*sizey*sizez*8 )  
  write(23,rec=1) cvydata
  !Now the circumplanetary data is stored in a binary file
  close(23)

  gasfield = 'gasvz'
  call createfieldbinaryfilename(iter,gasfield,fileextracter)
  table = 'extracter_files/'//trim(fileextracter)
  open(unit=23, status='unknown',file=table,form="unformatted", &
       access="direct", recl= sizex*sizey*sizez*8 )  
  write(23,rec=1) cvzdata
  !Now the circumplanetary data is stored in a binary file
  close(23)

  gasfield = 'gasenergy'
  call createfieldbinaryfilename(iter,gasfield,fileextracter)
  table = 'extracter_files/'//trim(fileextracter)
  open(unit=23, status='unknown',file=table,form="unformatted", &
       access="direct", recl= sizex*sizey*sizez*8 )  
  write(23,rec=1) cendata
  close(23)

!----------------------------------------------------------------------
!       All the data is stored in hard disk now        
!----------------------------------------------------------------------

  !Make a jump as big as the number of processors
  iter = iter + iterjump * nprocs

!----------------------------------------------------------------------
!    At this point all the variables that ocupates RAM are recycled
!----------------------------------------------------------------------

  enddo

!----------------------------------------------------------------------
!    Deallocate all the vectors
!----------------------------------------------------------------------

  deallocate(dendata )
  deallocate(vxdata  )
  deallocate(vydata  )
  deallocate(vzdata  )
  deallocate(endata  )
  deallocate(cdendata)
  deallocate(cvxdata )
  deallocate(cvydata )
  deallocate(cvzdata )
  deallocate(cendata )

  !Wait until all the processors finish
  call MPI_BARRIER(MPI_COMM_WORLD,ierr) 

  !END MPI CALL
  call MPI_FINALIZE(ierr)

!----------------------------------------------------------------------
!    THE PROGRAM ENDS HERE
!----------------------------------------------------------------------

end program EXTRACTER
