!This subroutine reads the input file to be called inside the main program, all
!the data must be in the next order
!xmin,
!xmax,
!ymin,
!ymax,
!zmin,
!zmax

subroutine read_input(NXYZ,itermin,itertot,iterjump,limits,ptype,ind,mydirsetup)

 implicit none

 !output parameters
 real*8 :: limits(6) 
 character(1) :: ptype
 character(99) :: mydirsetup
 integer :: ind, NXYZ(3), itermin, itertot, iterjump

 !the input file must be called input.dat
 open(13,file='input.dat',status='old')
 !start to read the parameters
   read(13,*),NXYZ(1)
   read(13,*),NXYZ(2)
   read(13,*),NXYZ(3)
   read(13,*),itermin
   read(13,*),itertot
   read(13,*),iterjump
   read(13,*),limits(1)
   read(13,*),limits(2)
   read(13,*),limits(3)
   read(13,*),limits(4)
   read(13,*),limits(5)
   read(13,*),limits(6)
   read(13,*),ptype
   read(13,*),ind
   read(13,*),mydirsetup
 close(13) 

end subroutine read_input


!This subroutine reads the input file to be called inside the circumplanetary
!program, all the data must be in the next order
!xmin,
!xmax,
!ymin,
!ymax,
!zmin,
!zmax,
!NXmin,
!NXmax,
!NYmin,
!NYmax,
!NZmin,
!NZmax,

subroutine read_input_circumplanetary(NXYZ,itertot,limits,ENES)

 implicit none

 !output parameters
 !real limits of X,Y,Z
 real*8 :: limits(6) 
 !idex where the cut must be done
 integer :: ENES(6), NXYZ(3), itertot

 !the input file must be called input.dat
 open(13,file='input_circumplanetary.dat',status='old')
 !start to read the parameters
   read(13,*),NXYZ(1)
   read(13,*),NXYZ(2)
   read(13,*),NXYZ(3)
   read(13,*),itertot
   read(13,*),limits(1)
   read(13,*),limits(2)
   read(13,*),limits(3)
   read(13,*),limits(4)
   read(13,*),limits(5)
   read(13,*),limits(6)
   read(13,*),ENES(1)
   read(13,*),ENES(2)
   read(13,*),ENES(3)
   read(13,*),ENES(4)
   read(13,*),ENES(5)
   read(13,*),ENES(6)
 close(13) 

end subroutine read_input_circumplanetary
