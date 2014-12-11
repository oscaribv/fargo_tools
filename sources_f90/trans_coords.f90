!This subroutine transforms cylindrical to cartesian coordinates.
!All the data cube!!!!
!The input parameters are
!       NX -> Dimenssion in \theta
!       NY -> Dimenssion in R
!       NZ -> Dimenssion in Z
!       data(NX*NY*NZ) containing the values
!       minimum and maximum values for r, theta and z
!       the name of the output file
subroutine cyl2carte(NX,NY,NZ,dendata,vrdata,vtdata,vzdata,endata,rmin,rmax,thmin,thmax,zmin,zmax,filename)

 implicit none
   !Input variables
   integer :: NX, NY, NZ
   real*8, dimension(NX*NY*NZ) :: dendata, vrdata, vtdata, vzdata, endata
   real*8 :: rmin,rmax,thmin,thmax,zmin,zmax
   character(100) :: filename
   !Local variables
   integer :: i,j,k
   real*8 :: r, dr, th, dth, z, dz
   real*8 :: den, vr, vt, vz, en

   !Open the file
   open(unit=101, file=filename, status="unknown")

   !Calculates the size of the jumps
   dth = (thmax - thmin) / real(NX-1) !It must be periodic!
   dr = (rmax - rmin) / real(NY) 
   dz = (zmax - zmin) / real(NZ)

   !Where do we start to jump?
   r = rmin
   th = thmin
   z = zmin


   write(101,*),'X, Y, Z, log(den), vr, vtheta, vz, log(energy)'
   !Let's fill the file
   do k=1,NZ
     do j=1,NY
        do i=1,NX
          den = log10(dendata(i+(j-1)*NX+(k-1)*NX*NY))
          vr = vrdata(i+(j-1)*NX+(k-1)*NX*NY)
          vt = vtdata(i+(j-1)*NX+(k-1)*NX*NY)
          vz = vzdata(i+(j-1)*NX+(k-1)*NX*NY)
          en = log10(endata(i+(j-1)*NX+(k-1)*NX*NY))
          
          write(101,*) r*cos(th),',', r*sin(th),',', z,',', den,',', vr,',', vt,',', vz,',', en
          th = th + dth
        end do
        th = thmin
        r = r + dr
     end do
     r = rmin
     z = z + dz
   end do

  !Close the file
  close(101)

end subroutine cyl2carte

!This subroutine transforms spherical to cartesian coordinates.
!All the data cube!!!!
!The input parameters are
!       NX -> Dimenssion in azimuth
!       NY -> Dimenssion in radius
!       NZ -> Dimenssion in colatitude
!       data(NX*NY*NZ) containing the values
!       minimum and maximum values for r, azimuth and colatitude
!       the name of the output file
subroutine sph2carte(NX,NY,NZ,dendata,vrdata,vazdata,vcoldata,endata,rmin,rmax,azmin,azmax,colmin,colmax,filename)

 implicit none
   !Input variables
   integer :: NX, NY, NZ
   real*8, dimension(NX*NY*NZ) :: dendata, vrdata, vazdata, vcoldata, endata
   real*8 :: rmin,rmax,azmin,azmax,colmin,colmax
   character(100) :: filename
   !Local variables
   integer :: i,j,k
   real*8 :: r, dr, az, dtaz, col, dcol
   real*8 :: x,y, z
   real*8 :: den, vr, vaz, vcol, en
   real*8 :: vx, vy, vz
   real*8 :: dxdr, dxdaz, dxdcol
   real*8 :: dydr, dydaz, dydcol
   real*8 :: dzdr, dzdaz, dzdcol

   !Open the file
   open(unit=101, file=filename, status="unknown")

   !Calculates the size of the jumps
   dtaz = (azmax - azmin) / real(NX-1) !It must be periodic!
   dr = (rmax - rmin) / real(NY) 
   dcol = (colmax - colmin) / real(NZ)

   !Where do we start to jump?
   r = rmin
   az = azmin
   col = colmin

   write(101,*),'X, Y, Z, log(den), vx, vy, vz, log(energy)'
   !Let's fill the file
   do k=1,NZ
     do j=1,NY
        do i=1,NX
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

          !lets transform the velocities to carseian veolicities to be used in
          !paravew easily and to plot velocity maps
          !creating transformation coefficients
          dxdr = sin(col) * cos(az)
          dxdaz =  - r * sin(col) * sin(az)
          dxdcol = r * cos (col) * cos(az)

          dydr = sin(col) * sin(az)
          dydaz = r * sin(col) * cos(az)
          dydcol = r * cos(col) * sin(az)

          dzdr = cos(col)
          dzdaz = 0.0
          dzdcol = - r * sin(col)

          !Let's transform the velocities
          vx = dxdr * vr + dxdcol * vcol + dxdaz * vaz         
          vy = dydr * vr + dydcol * vcol + dydaz * vaz         
          vz = dzdr * vr + dzdcol * vcol + dzdaz * vaz         
          !DONE!

          !write(101,*) x,',', y,',', z,',', den,',', vr,',', vaz,',', vcol,',', en
          write(101,*) x,',', y,',', z,',', den,',', vx,',', vy,',', vz,',', en
          az = az + dtaz
        end do
        az = azmin
        r = r + dr
     end do
     r = rmin
     col = col + dcol
   end do

  !Close the file
  close(101)

end subroutine sph2carte


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
