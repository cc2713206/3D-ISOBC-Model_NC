MODULE param
implicit none
real(8) :: tdp
real(8), parameter :: hoch=140.
real(8), parameter :: loch=8000.
real(8), parameter :: pw=2000. 
real(8) :: slip, speed
INTEGER(4), PARAMETER :: nx = 81 ! number of grid points x-direction
INTEGER(4), PARAMETER :: ny = 71 ! number of grid points y-direction
INTEGER(4), PARAMETER :: nz = 50 ! number of grid points z-direction
real(8) :: tbu(0:ny+1,0:nx+1), tbv(0:ny+1,0:nx+1)
real(8) :: depth(0:ny+1,0:nx+1)
integer, parameter :: ni=10 
REAL(8), parameter :: dx = 750.
real(8), parameter :: dy = 400.
real(8), parameter :: dz = 4.
real(8), parameter :: dt = 10.
real(8), PARAMETER :: bta = 158. ! angle between the north and across-slope direction
real(8), PARAMETER :: tta = -75.2 ! latitude
real(8), PARAMETER :: gv = 9.81 ! acceleration due to gravity
real(8), parameter :: pi=3.14159265359
real(8), parameter :: om=7.292e-5
real(8) :: f ! Coriolis parameter
real(8), parameter :: uini=0.060 
real(8), parameter :: fini=1.d-5 
integer, parameter :: isdp=7
real(8), parameter :: dsf = 0.6d-10 
real(8), parameter :: slp=0.0075 
real(8), parameter :: ar=1./50.
real(8), parameter :: rho0=1030.
real(8), parameter :: rhoi=920.
real(8), parameter :: t0=-2.
real(8), parameter :: s0=34.5
real(8), parameter :: bt=3.87e-5
real(8), parameter :: bs=7.86e-4
real(8), parameter :: cmin=0. ! switch of frazil ice module
real(8), parameter :: ib=730.
real(8), parameter :: smix=34.26
real(8), parameter :: slow=34.36
real(8), parameter :: tmix=-0.10
real(8), parameter :: tlow=0.18
real(8), parameter :: tsur=-1.9
real(8), parameter :: ssur=33.8
real(8), parameter :: l1=-0.0573
real(8), parameter :: l2=0.0832
real(8), parameter :: l3=-7.61e-4
real(8), parameter :: c0=3974.
real(8), parameter :: ci=2009.
real(8), parameter :: ti=-15. 
real(8), parameter :: lf=3.35e5
real(8), parameter :: r1=0.025
real(8), parameter :: r2=0.050
real(8), parameter :: r3=0.100
real(8), parameter :: r4=0.4
real(8), parameter :: r5=0.7
real(8), parameter :: r6=1.0
real(8), parameter :: r7=1.3
real(8), parameter :: r8=1.6
real(8), parameter :: r9=1.9
real(8), parameter :: r10=2.2
real(8), parameter :: kt=1.4e-7
real(8), parameter :: ks=8.d-10
real(8), parameter :: theta=0.05 
real(8), parameter :: azmin=1.95d-6
real(8), parameter :: azmin2=1.4d-7
real(8), parameter :: ekmin=1.d-9 
real(8), parameter :: epmin=1.d-12
real(8), parameter :: azmax=0.5
real(8), parameter :: nni=1000. 

! parameters for S.O.R iteration
real(8) :: ustar(0:nz+1,0:ny+1,0:nx+1), vstar(0:nz+1,0:ny+1,0:nx+1), wstar(0:nz+1,0:ny+1,0:nx+1), qstar(nz,ny,nx)
real(8), parameter :: omega=1.4 
real(8), parameter :: peps=0.1e-1 
real(8) :: at(nz,ny,nx), ab(nz,ny,nx), ae(nz,ny,nx), aw(nz,ny,nx)
real(8) :: atot(nz,ny,nx), an(nz,ny,nx), as(nz,ny,nx) 
real(8) :: drdz, drdx, drdy 
real(8) :: ah(0:nz+1,0:ny+1,0:nx+1), az(0:nz+1,0:ny+1,0:nx+1), taz(0:nz+1,0:ny+1,0:nx+1), saz(0:nz+1,0:ny+1,0:nx+1) ! eddy viscosities
INTEGER :: n,i,j,k,ll
real(8) :: uu, vv, ww
real(8) ::  alpha

real(8) :: time
real(8) :: rho(0:nz+1,0:ny+1,0:nx+1)
real(8) :: p(0:nz+1,0:ny+1,0:nx+1) ! hydrostatic pressure
real(8) :: dq(0:nz+1,0:ny+1,0:nx+1) ! pressure correction
real(8) :: q(0:nz+1,0:ny+1,0:nx+1) ! nonhydrostatic pressure
real(8) :: u(0:nz+1,0:ny+1,0:nx+1)
real(8) :: v(0:nz+1,0:ny+1,0:nx+1)
real(8) :: w(0:nz+1,0:ny+1,0:nx+1)
real(8) :: un(0:nz+1,0:ny+1,0:nx+1), vn(0:nz+1,0:ny+1,0:nx+1), wn(0:nz+1,0:ny+1,0:nx+1)
real(8) :: usum(0:ny+1,0:nx+1), vsum(0:ny+1,0:nx+1) ! depth-integrated flow
LOGICAL :: wet(0:nz+1,0:ny+1,0:nx+1), dry(0:nz+1,0:ny+1,0:nx+1), optf ! wet/dry pointer

! arrays for advection subroutine
real(8) :: CuP(0:nz+1,0:ny+1,0:nx+1), CuN(0:nz+1,0:ny+1,0:nx+1)
real(8) :: CvP(0:nz+1,0:ny+1,0:nx+1), CvN(0:nz+1,0:ny+1,0:nx+1)
real(8) :: CwP(0:nz+1,0:ny+1,0:nx+1), CwN(0:nz+1,0:ny+1,0:nx+1)
real(8) :: B(0:nz+1,0:ny+1,0:nx+1), BN(0:nz+1,0:ny+1,0:nx+1)

real(8) :: rhoa, wi(ni), nu(ni), sh(ni), ri(ni), re(ni), uc2(ni), vi(ni), wt(ni), dv(ni)
real(8) :: tfc(0:nz+1,0:ny+1,0:nx+1), fc(ni,0:nz+1,0:ny+1,0:nx+1), t(0:nz+1,0:ny+1,0:nx+1), s(0:nz+1,0:ny+1,0:nx+1)
real(8) :: tmel(0:ny+1,0:nx+1), mel(0:ny+1,0:nx+1), tpre(0:ny+1,0:nx+1), hsc(0:ny+1,0:nx+1), pre(ni,0:ny+1,0:nx+1)
real(8) :: uss(0:nz+1,0:ny+1,0:nx+1), ek(0:nz+1,0:ny+1,0:nx+1), ep(0:nz+1,0:ny+1,0:nx+1)
real(8) :: atop, abot, ps(0:nz+1,0:ny+1,0:nx+1), pb(0:nz+1,0:ny+1,0:nx+1)
real(8) :: ttfc(0:nz+1,0:ny+1,0:nx+1), sfc(0:nz+1,0:ny+1,0:nx+1), senu(ni,0:nz+1,0:ny+1,0:nx+1), wf0(ni,0:nz+1,0:ny+1,0:nx+1), wf(ni,0:nz+1,0:ny+1,0:nx+1)
real(8) :: dzz

real(8) :: cs, sn, dp(0:ny+1,0:nx+1), tfb(0:ny+1,0:nx+1)
real(8) :: csb, snb, cst, snt

LOGICAL :: land(0:ny+1,0:nx+1) ! land pointer

INTEGER :: ntot, nout, nout2, ncl, zot, zom, zob
real(8) :: ntt, ntt2
integer :: ntem
character( len = 3 ) :: cTemp, ccl 

real(8) :: ekn(0:nz+1,0:ny+1,0:nx+1), epn(0:nz+1,0:ny+1,0:nx+1), fcn(ni,0:nz+1,0:ny+1,0:nx+1), tn(0:nz+1,0:ny+1,0:nx+1),ssn(0:nz+1,0:ny+1,0:nx+1)
real(8) :: sqcd, azn(0:ny+1,0:nx+1)
real(8) :: loca(0:ny+1,0:nx+1), drhos, tnn(0:ny+1,0:nx+1)
real(8) :: rhopini
real(8) :: soli(ny), wbb(ny,nx), etst
integer :: isoli(ny), isl
real(8) :: rch, n22(0:nz+1,0:ny+1,0:nx+1), duvz, ttbb(0:ny+1,0:nx+1), ssbb(0:ny+1,0:nx+1), usuc
real(8) :: cii, lff, kzt, kzs, sanj, wujx, shiz, yuel, lingx
integer :: x1, x2, x3, x4, x5
real(8) :: ad, adh, adu, uad, rhoini(0:nz+1,0:ny+1,0:nx+1)
integer :: isp, js1, js2, js3, js4, isd, loci(0:ny+1,0:nx+1)
real(8) :: aup, ach, adw, fflux(1:nx), tflux(1:nx), bflux(1:nx), upmel, chmel, dwmel, smel, aupmel, achmel, adwmel, asmel, affl, atfl, abfl, atq, abq
real(8) :: jif1, jif2, wei(0:nz+1,0:ny+1,0:nx+1)
real(8), parameter :: prh=10.
real(8), parameter :: uscr=0.0010
real(8) :: avbp(0:ny+1,0:nx+1), avup(0:ny+1,0:nx+1), avvp(0:ny+1,0:nx+1), avcp(0:ny+1,0:nx+1), ujif, vjif 
END MODULE param