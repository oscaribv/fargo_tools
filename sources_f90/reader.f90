! Example on how to open a FARGO3D output with fortran. Be
! free to adapt it.
! Compilation line: gfortran -cpp reader.f90

!Parameters:
#define NX 384
#define NY 128
#define NZ 10
#define FIELDNAME "gasdens1.dat"

program FARGO3D_READER
  integer :: i,j,k
  real*8  :: data(NX*NY*NZ)
  
  !Opening binary file
  open(unit=100,            &
       status="old",        &
       file=FIELDNAME,      &
       form="unformatted",  &
       access="direct",     &
       recl = NX*NY*NZ*8)

  !Reading file
  read(100,rec=1) data
  close(100)

  !Saving it as an ascii file (Pixels Vs Field). The format is
  !friendly with gnuplot splot wl function (Note that gnuplot can read
  !binaries, see the script reader.gpl). You should adapt it in order
  !to work with physical coordinates. This version is only for 2D
  !outputs
  
  open(unit=101,         &
       file="table.txt", &
       status="unknown")
 
  do k = 1, NZ 
  do j=1,NY
     do i=1,NX
        write(101,*)i,j,k, log10(data(i+(j-1)*NX+(k-1)*NX*NY))
     end do
     write(101,*)
  end do
  end do

  close(101)
  
  !Now, for display the output: 
  !Inside gnuplot --> sp "table.txt" w l
  
end program FARGO3D_READER
