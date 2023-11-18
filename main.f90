!=======================================
!3D ISOBC MODEL V1.0
!=======================================

PROGRAM full3D

USE param
USE sub
implicit none

x1=int(1500./dx)+1
x2=int(6000./dx)+1
x3=int(15000./dx)+1
x4=int(30000./dx)+1
x5=int(60000./dx)+1

! runtime parameters
ntot = INT(90.*24.*3600./dt)+1
time = 0.0
! output frequency
ntt=360./60.
nout = INT(ntt*3600./dt)
ntt2=120./60.
nout2 = INT(ntt2*3600./dt)

open (10, file = 'melseries.dat' ) 
open (11, file = 'fluxseries.dat' )

!**********
CALL INIT  ! initialisation
!**********
!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot
  time = time + dt  
  write(6,*)"time (hours)", time/(3600.)
  ad = 1.
  adh= MIN(time/(3.*24.*3600.),1.0)
  adu= MIN(time/(3.*24.*3600.),1.0)
  isd=1+int((isdp-1)*adh)
  ! prognostic equations
  CALL dyn
END DO ! end of iteration loop
END PROGRAM full3D