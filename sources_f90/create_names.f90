!-------------------------------------------------------------
! Producing name for binary files, the field name as input 
! This subroutina replaces the previous ones for density
! vx, vy, vz, and energy
!-------------------------------------------------------------

     subroutine createfieldbinaryfilename(i,fieldname,filename)
!    input variables
     integer           :: i
!    output variables
        !Field name can be gasdens, gasvx, gasvy, etc
     character(len=15) :: fieldname
     character(len=30) :: filename
!    subroutine variables
     character(len=4)  :: ix

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename=trim(fieldname)//trim(ix)//'.dat'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename=trim(fieldname)//trim(ix)//'.dat'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename=trim(fieldname)//trim(ix)//'.dat'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename=trim(fieldname)//trim(ix)//'.dat'
        endif

     end subroutine

!-------------------------------------------------------------
! Producing name for density as a human readable file .txt
! Outputs files containing X Y Z and density
!-------------------------------------------------------------

     subroutine createdensitytxtfilename(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='gasdens000'//trim(ix)//'.txt'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='gasdens00'//trim(ix)//'.txt'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='gasdens0'//trim(ix)//'.txt'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='gasdens'//trim(ix)//'.txt'
        endif

     end subroutine

!-------------------------------------------------------------
! Producing names of csv files
!-------------------------------------------------------------

     subroutine create_csv_files(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='disk'//trim(ix)//'.csv'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='disk'//trim(ix)//'.csv'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='disk'//trim(ix)//'.csv'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='disk'//trim(ix)//'.csv'
        endif

     end subroutine


!----------------------------------------------------------------------
! Producing name for binary density file for the circumplanetary disk
!-----------------------------------------------------------------------

     subroutine create_circumplanetary_den_binaryfilename(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        call system('cd circumplanetary_data')

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='gasdens'//trim(ix)//'.dat'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='gasdens'//trim(ix)//'.dat'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='gasdens'//trim(ix)//'.dat'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='gasdens'//trim(ix)//'.dat'
        endif

        call system('cd ..')

     end subroutine
!----------------------------------------------------------------------
! Producing name for binary vx file for the circumplanetary disk
!-----------------------------------------------------------------------

     subroutine create_circumplanetary_vx_binaryfilename(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        call system('cd circumplanetary_data')

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='gasvx'//trim(ix)//'.dat'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='gasvx'//trim(ix)//'.dat'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='gasvx'//trim(ix)//'.dat'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='gasvx'//trim(ix)//'.dat'
        endif

        call system('cd ..')

     end subroutine
!----------------------------------------------------------------------
! Producing name for binary vy file for the circumplanetary disk
!-----------------------------------------------------------------------

     subroutine create_circumplanetary_vy_binaryfilename(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        call system('cd circumplanetary_data')

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='gasvy'//trim(ix)//'.dat'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='gasvy'//trim(ix)//'.dat'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='gasvy'//trim(ix)//'.dat'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='gasvy'//trim(ix)//'.dat'
        endif

        call system('cd ..')

     end subroutine

!----------------------------------------------------------------------
! Producing name for binary vz file for the circumplanetary disk
!-----------------------------------------------------------------------

     subroutine create_circumplanetary_vz_binaryfilename(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        call system('cd circumplanetary_data')

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='gasvz'//trim(ix)//'.dat'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='gasvz'//trim(ix)//'.dat'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='gasvz'//trim(ix)//'.dat'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='gasvz'//trim(ix)//'.dat'
        endif

        call system('cd ..')

     end subroutine

!----------------------------------------------------------------------
! Producing name for binary energy file for the circumplanetary disk
!-----------------------------------------------------------------------

     subroutine create_circumplanetary_ene_binaryfilename(i,filename)
!    input variables
     integer           :: i
!    output variables
     character(len=15) :: filename
!    subroutine variables
     character(len=4)  :: ix

        call system('cd circumplanetary_data')

        if (i.lt.10) then
              write(ix,'(I1)') i
        filename='gasenergy'//trim(ix)//'.dat'
        endif
        if ((i.ge.10).and.(i.lt.100)) then
              write(ix,'(I2)') i
        filename='gasenergy'//trim(ix)//'.dat'
        endif
        if ((i.ge.100).and.(i.lt.1000)) then
              write(ix,'(I3)') i
        filename='gasenergy'//trim(ix)//'.dat'
        endif
        if ((i.ge.1000)) then
              write(ix,'(I4)') i
        filename='gasenergy'//trim(ix)//'.dat'
        endif

        call system('cd ..')

     end subroutine


!----------------------------------------------------------------------
!This subroutine will create the p2cplot.gpl authomatically, with the
!aim to adapt it to the type of plot
!----------------------------------------------------------------------

subroutine create_gpl_file(NX,NY,NZ,limits,ptype,ind)

 implicit none

 real*8 :: limits(6)
 integer :: NX,NY,NZ
 character(1) :: ptype
 integer :: ind 
 real*4 :: cut, vmin, vmax, rsize

 !check what cut is
 if ( ptype == 'z' ) then
  vmax = limits(6); vmin = limits(5)
  cut = ( vmax - vmin ) / (NZ - 1.0)
  cut = vmin + (ind-1)*cut
 else if ( ptype == 't' ) then
  vmax = limits(2); vmin = limits(1)
  cut = ( vmax - vmin ) / (NX - 1.0)
  cut = vmin + (ind-1)*cut
 end if

 !open the file
 open(17,file='p2cplot.gpl',status='unknown')

 write(17,*),'#This file was created automathically by the main program'
 write(17,*),'#You can modify it if you want'
 write(17,*),''
 write(17,*),'set terminal pngcairo size 800,800'
 write(17,*),'set output ''images/cartesiandiskxxxx.png'' '
! write(17,*),'set lmargin at screen 0.05'
! write(17,*),'set rmargin at screen 0.85'
! write(17,*),'set bmargin at screen 0.1'
! write(17,*),'set tmargin at screen 0.9'
 write(17,*),'set pm3d map'
 write(17,*),'unset key'
 write(17,*),'set multiplot'
 write(17,*),'set parametric'
 write(17,*),'set isosamples 500'
 write(17,*),'set datafile separator '','''
 !the next depends in z or theta
 if ( ptype == 'z' ) then
   !size of r for z plots
   rsize = real(limits(4))
   write(17,*),'set xr [-',rsize*1.1,':',rsize*1.1,']'
   write(17,*),'set yr [-',rsize*1.1,':',rsize*1.1,']'
   write(17,*),'set xl ''X'''
   write(17,*),'set yl ''Y'''
   write(17,*),'set title ''Cut at Z = ',cut,' '''
 else if ( ptype == 't' ) then
   write(17,*),'set xr [',limits(3),':',limits(4),']'
   write(17,*),'set yr [',limits(5),':',limits(6),']'
   write(17,*),'set xl ''r'''
   write(17,*),'set yl ''z'''
   write(17,*),'set title ''Cut at theta = ',cut,' '''
 end if
 write(17,*),'set palette rgb 34,35,36'
 write(17,*),'filename = ''cartesiantxt/gasdensxxxx.txt'''
 write(17,*),'splot filename u 1:2:3'

 close(17)

end subroutine create_gpl_file
