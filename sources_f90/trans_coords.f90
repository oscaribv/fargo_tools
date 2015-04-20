!----------------------------------------------------------------------------------------------------------
!                                               carte2carte
!                                        Written by O. Barragán
!----------------------------------------------------------------------------------------------------------
!This subroutine writes a file in cartesian coordiantes into a csv file 
!The input vaules are:
!  NX, NY, NZ                                   -> Number of points in x, y and z
!  dendata, vxdata, vydata, vzdata, endata      -> The field arrays, each one of these has a size NX*NY*NZ
!  xmin,xmax,ymin,ymax,zmin,zmax              -> the physical sizes of the "box"
!  filename                                     -> name of the output file (csv format)
!----------------------------------------------------------------------------------------------------------
! The output file is a csv file with the columns 
!        X, Y, Z, log(den), vx, vy, vz, log(energy)
! These files can be open with 3D-visulization programs, such as Paraview
!----------------------------------------------------------------------------------------------------------

!Define the subroutine name carte2carte
subroutine carte2carte(NX,NY,NZ,dendata,vydata,vxdata,vzdata,endata,ymin,ymax,xmin,xmax,zmin,zmax,filename)

!START VARIABLE DEFINITIONS
 implicit none
   !Input variables
   integer :: NX, NY, NZ
   real*8, dimension(NX*NY*NZ) :: dendata, vxdata, vydata, vzdata, endata
   real*8 :: xmin,xmax,ymin,ymax,zmin,zmax
   character(30) :: filename
   !Local variables
   integer :: i,j,k
   real*8 :: x, dx, y, dy, z, dz
   real*8 :: den, en
   real*8 :: vx, vy, vz
!END VARIABLE DEFINITIONS

   !Open the file, this will be the output file
   open(unit=101, file=trim(filename), status="unknown")

   !Calculates the size of the jumps
   dx = (xmax - xmin) / real(NX)
   dy = (ymax - ymin) / real(NY) 
   dz = (zmax - zmin) / real(NZ)

   !Where do we start to jump?
   x = xmin
   y = ymin
   z = zmin

   !This line writes the head of the csv file
   write(101,*),'X, Y, Z, lden, vx, vy, vz, lenergy'

   !Let's fill the file
   do k=1,NZ
     do j=1,NY
        do i=1,NX

          !Whose are the values of the fields in this point i,j,k?
          !density and energy are transformed to log10
          den = log10(dendata(i+(j-1)*NX+(k-1)*NX*NY))
          vx = vxdata(i+(j-1)*NX+(k-1)*NX*NY)
          vy = vydata(i+(j-1)*NX+(k-1)*NX*NY)
          vz = vzdata(i+(j-1)*NX+(k-1)*NX*NY)
          en = log10(endata(i+(j-1)*NX+(k-1)*NX*NY))

          !There are no transformation

          !Now I know the values of the coordinates in cartesian coordinates,
          !and also the vector components, lets write it in the output file
          !before I forger it          
          write(101,13) x,',',  y,',', z,',', den,',', vx,',', vy,',', vz,',', en
          13 format(E14.7,a1,E14.7,a1,E14.7,a1,E14.7,a1,E14.7,a1, &
                    E14.7,a1,E14.7,a1,E14.7,a1)

          !Jump in x
          x = x + dx

        end do
        !Restart x
        x = xmin
        !Jump in y
        y = y + dy
     end do
     !Restart r
     y = ymin
     !Jump in z
     z = z + dz
   end do

   !Now the output file is filled with all the values

  !Close the file
  close(101)

end subroutine carte2carte
!End of the subroutine carte2carte



!----------------------------------------------------------------------------------------------------------
!                                               cyl2carte
!                                        Written by O. Barragán
!----------------------------------------------------------------------------------------------------------
!This subroutine transform cylindrical coordinates to cartesian ones
!The input vaules are:
!  NX, NY, NZ                                   -> Number of points in azimuth, radius and z
!  dendata, vrdata, vtdata, vzdata, endata      -> The field arrays, each one of these has a size NX*NY*NZ
!  rmin,rmax,thmin,thmax,zmin,zmax              -> the physical sizes of the "cylindrical box"
!  filename                                     -> name of the output file (csv format)
!----------------------------------------------------------------------------------------------------------
! The output file is a csv file with the columns 
!        X, Y, Z, log(den), vx, vy, vz, log(energy)
! These files can be open with 3D-visulization programs, such as Paraview
!----------------------------------------------------------------------------------------------------------

!Define the subroutine name cyl2carte
subroutine cyl2carte(NX,NY,NZ,dendata,vrdata,vtdata,vzdata,endata,rmin,rmax,thmin,thmax,zmin,zmax,filename)

!START VARIABLE DEFINITIONS
 implicit none
   !Input variables
   integer :: NX, NY, NZ
   real*8, dimension(NX*NY*NZ) :: dendata, vrdata, vtdata, vzdata, endata
   real*8 :: rmin,rmax,thmin,thmax,zmin,zmax
   character(30) :: filename
   !Local variables
   integer :: i,j,k
   real*8 :: r, dr, th, dth, z, dz
   real*8 :: den, vr, vt, vz, en
   real*8 :: vx, vy 
!END VARIABLE DEFINITIONS

   !Open the file, this will be the output file
   open(unit=101, file=trim(filename), status="unknown")

   !Calculates the size of the jumps
   dth = (thmax - thmin) / real(NX-1) !It must be periodic!
   dr = (rmax - rmin) / real(NY) 
   dz = (zmax - zmin) / real(NZ)

   !Where do we start to jump?
   r = rmin
   th = thmin
   z = zmin

   !This line writes the head of the csv file
   write(101,*),'X, Y, Z, lden, vx, vy, vz, lenergy'

   !Let's fill the file
   do k=1,NZ
     do j=1,NY
        do i=1,NX

          !Whose are the values of the fields in this point i,j,k?
          !density and energy are transformed to log10
          den = log10(dendata(i+(j-1)*NX+(k-1)*NX*NY))
          vr = vrdata(i+(j-1)*NX+(k-1)*NX*NY)
          vt = vtdata(i+(j-1)*NX+(k-1)*NX*NY)
          vz = vzdata(i+(j-1)*NX+(k-1)*NX*NY)
          en = log10(endata(i+(j-1)*NX+(k-1)*NX*NY))

          !transforming the velocity components, these have to be transformed as
          !vectors components, more details in O. Barragán, 2015, Master thesis
          vx = vr * cos(th) - vt * r * sin(th)
          vy = vr * sin(th) + vt * r * cos(th)
          !vz = vz

          !Now I know the values of the coordinates in cartesian coordinates,
          !and also the vector components, lets write it in the output file
          !before I forger it          
          write(101,13) r*cos(th),',', r*sin(th),',', z,',', den,',', vx,',', vy,',', vz,',', en
          13 format(E14.7,a1,E14.7,a1,E14.7,a1,E14.7,a1,E14.7,a1, &
                    E14.7,a1,E14.7,a1,E14.7,a1)

          !Jump in azimuth
          th = th + dth

        end do
        !Restart azimuth
        th = thmin
        !Jump in r
        r = r + dr
     end do
     !Restart r
     r = rmin
     !Jump in z
     z = z + dz
   end do

   !Now the output file is filled with all the values

  !Close the file
  close(101)

end subroutine cyl2carte
!End of the subroutine cyl2carte

!----------------------------------------------------------------------------------------------------------
!                                            sph2carte
!                                       Written by O. Barragán
!----------------------------------------------------------------------------------------------------------
!This subroutine transform spherical coordinates to cartesian ones
!The input vaules are:
!  NX, NY, NZ                                   -> Number of points in azimuth, radius and colatitude
!  dendata, vrdata, vazdata, vcoldata, endata      -> The field arrays, each one of these has a size NX*NY*NZ
!  rmin,rmax,thmin,thmax,zmin,zmax              -> the physical sizes of the "spherical box"
!  filename                                     -> name of the output file (csv format)
!----------------------------------------------------------------------------------------------------------
! The output file is a csv file with the columns 
!        X, Y, Z, log(den), vx, vy, vz, log(energy)
! These files can be open with 3D-visulization programs, such as Paraview
!----------------------------------------------------------------------------------------------------------

!Define the subroutine name sph2carte
subroutine sph2carte(NX,NY,NZ,dendata,vrdata,vazdata,vcoldata,endata,rmin,rmax,azmin,azmax,colmin,colmax,filename)

!START VARIABLE DEFINITIONS
 implicit none
   !Input variables
   integer :: NX, NY, NZ
   real*8, dimension(NX*NY*NZ) :: dendata, vrdata, vazdata, vcoldata, endata
   real*8 :: rmin,rmax,azmin,azmax,colmin,colmax
   character(30) :: filename
   !Local variables
   integer :: i,j,k
   real*8 :: r, dr, az, dtaz, col, dcol
   real*8 :: x,y, z
   real*8 :: den, vr, vaz, vcol, en
   real*8 :: vx, vy, vz
   real*8 :: dxdr, dxdaz, dxdcol
   real*8 :: dydr, dydaz, dydcol
   real*8 :: dzdr, dzdaz, dzdcol
!END VARIABLE DEFINITIONS

   !Open the file, this will be the output file
   open(unit=101, file=trim(filename), status="unknown")

   !Calculates the size of the jumps
   dtaz = (azmax - azmin) / real(NX-1) !It must be periodic!
   dr = (rmax - rmin) / real(NY) 
   dcol = (colmax - colmin) / real(NZ)

   !Where do we start to jump?
   r = rmin
   az = azmin
   col = colmin

   !This line writes the head of the csv file
   write(101,*),'X, Y, Z, lden, vx, vy, vz, lenergy'
   !Let's fill the file
   do k=1,NZ
     do j=1,NY
        do i=1,NX

          !Whose are the values of the fields in this point i,j,k?
          !density and energy are transformed to log10
          den = log10(dendata(i+(j-1)*NX+(k-1)*NX*NY))
          vr = vrdata(i+(j-1)*NX+(k-1)*NX*NY)
          vaz = vazdata(i+(j-1)*NX+(k-1)*NX*NY)
          vcol = vcoldata(i+(j-1)*NX+(k-1)*NX*NY)
          en = log10(endata(i+(j-1)*NX+(k-1)*NX*NY))

          !Obtain x,y,z from spherical coordinates          
          x = r * sin(col)
          y = x * sin(az)
          x = x * cos(az)
          z = r * cos(col)

          !transforming the velocity components, these have to be transformed as
          !vectors components, more details in O. Barragán, 2015, Master thesis

          !dx derivatives
          dxdr = sin(col) * cos(az)
          dxdaz =  - r * sin(col) * sin(az)
          dxdcol = r * cos (col) * cos(az)

          !dy derivatives
          dydr = sin(col) * sin(az)
          dydaz = r * sin(col) * cos(az)
          dydcol = r * cos(col) * sin(az)

          !dz derivatives
          dzdr = cos(col)
          dzdaz = 0.0
          dzdcol = - r * sin(col)

          !Let's transform the velocities
          vx = dxdr * vr + dxdcol * vcol + dxdaz * vaz         
          vy = dydr * vr + dydcol * vcol + dydaz * vaz         
          vz = dzdr * vr + dzdcol * vcol + dzdaz * vaz         

          !DONE!
          !The velocities are transformed now

          !Now I know the values of the coordinates in cartesian coordinates,
          !and also the vector components, lets write it in the output file
          !before I forger it          
          write(101,13) x,',', y,',', z,',', den,',', vx,',', vy,',', vz,',', en
          13 format(E14.7,a1,E14.7,a1,E14.7,a1,E14.7,a1,E14.7,a1, &
                    E14.7,a1,E14.7,a1,E14.7,a1)
          
          !Jump is azimuth
          az = az + dtaz

        end do
        !Restart azimuth
        az = azmin
        !Jump in r
        r = r + dr
     end do
     !Restart r
     r = rmin
     !Jump in colatitude
     col = col + dcol
   end do

  !Now the output file is filled with all the values
  !Close the file
  close(101)

end subroutine sph2carte
!End of subroutine sph2carte


!This subroutine transforms cylindrical to cartesian coordinates in a sheet of
!the Z plane
!The input parameters are
!       NX -> Dimenssion in \theta
!       NY -> Dimenssion in R
!       NZ -> Dimenssion in Z
!       k  -> The index of the Z where we want the cut
!       data(NX*NY*NZ) containing the values
!       minimum and maximum values for r, theta and z
!       the name of the output file
subroutine cut_at_Zk(NX,NY,NZ,k,data,rmin,rmax,thmin,thmax,zmin,zmax,filename)

 implicit none
   !Input variables
   integer :: NX, NY, NZ, k
   real*8, dimension(NX*NY*NZ) :: data
   real*8 :: rmin,rmax,thmin,thmax,zmin,zmax
   character*100 filename
   !Local variables
   integer :: i,j
   real*8 :: r, dr, th, dth, z, dz

   !This if avoids segmentation fault errors
   if ( k > NZ ) then
    print*,' k is greater than NZ, avoiding a segmentation fault'
    print*,' Chage this inside input.dat'
    stop
   end if

   !Calculates the size of the jumps
   dth = (thmax - thmin) / real(NX-1) !It must be periodic!
   dr = (rmax - rmin) / real(NY-1) 
   dz = (zmax - zmin) / real(NZ-1)


   !Where do we start to jump?
   r = rmin
   th = thmin
   z = zmin * (k-1)*dz ! find the z

   if ( NZ == 1 ) z = zmin

   !Open the file
   open(unit=101, file=filename, status="unknown")

   !Let's fill the file
     do j=1,NY
        do i=1,NX
          write(101,*) r*cos(th),',', r*sin(th),',', &
                       log10(data(i+(j-1)*NX+(k-1)*NX*NY))
          th = th + dth
        end do
        th = thmin
        r = r + dr
        write(101,*)
     end do

  !Close the file
  close(101)

end subroutine cut_at_Zk


!This subroutine transforms cylindrical to cartesian coordinates in a sheet of
!a given \theta
!The input parameters are
!       NX -> Dimenssion in \theta
!       NY -> Dimenssion in R
!       NZ -> Dimenssion in Z
!       i  -> The index of the i where we want the cut
!       data(NX*NY*NZ) containing the values
!       minimum and maximum values for r, theta and z
!       the name of the output file
subroutine cut_at_Ti(NX,NY,NZ,i,data,rmin,rmax,thmin,thmax,zmin,zmax,filename)

 implicit none
   !Input variables
   integer :: NX, NY, NZ, i
   real*8, dimension(NX*NY*NZ) :: data
   real*8 :: rmin,rmax,thmin,thmax,zmin,zmax
   character*100 filename
   !Local variables
   integer :: j,k
   real*8 :: r, dr, th, dth, z, dz

   !This if avoids segmentation fault errors
   if ( i > NX ) then
    print*,' i is greater than NX, avoiding a segmentation fault'
    print*,' Chage this inside input.dat'
    stop
   end if

   !Calculates the size of the jumps
   dth = (thmax - thmin) / real(NX-1) !It must be periodic!
   dr = (rmax - rmin) / real(NY-1) 
   dz = (zmax - zmin) / real(NZ-1)

   !Where do we start to jump?
   r = rmin
   th = thmin * (i-1)*dth !find the theta actual
   z = zmin 

   if ( NX == 1 ) th = thmin

   !Open the file
   open(unit=101, file=filename, status="unknown")

   !Let's fill the file
   do k=1, NZ
     do j=1, NY
          write(101,*) r, z, &
                       log10(data(i+(j-1)*NX+(k-1)*NX*NY))
        r = r + dr
     end do
     write(101,*)
     r = rmin
     z = z + dz
   end do

  !Close the file
  close(101)

end subroutine cut_at_Ti
