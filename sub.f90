MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
implicit none
! local parameters
real(8) :: di
snb=sin(bta*pi/180.)
csb=cos(bta*pi/180.)
snt=sin(tta*pi/180.)
cst=cos(tta*pi/180.)
cs=1./sqrt(1.+slp**2)
sn=slp/sqrt(1.+slp**2)
f=2.*om*(cst*snb*sn+snt*cs)
!ri(1)=r1; ri(2)=r2; ri(3)=r3; ri(4)=r4; ri(5)=r5; ri(6)=r6; ri(7)=r7; ri(8)=r8; ri(9)=r9; ri(10)=r10

!do i =1, ni
!  di=2*ri(i)
!  if(di<=1.27) wi(i)=2.025*di**1.621
!  if(di>1.27 .and. di<=7.) wi(i)=-0.103*di**2+4.069*di-2.024
!  wi(i)=wi(i)/1000.*cs
!  if(ri(i)<=1./13.8**0.5)   nu(i)=1.+0.17*ri(i)*13.8**0.5
!  if(ri(i)>1./13.8**0.5)  nu(i)=1.+0.55*ri(i)**(2./3.)*13.8**(1./3.)
!  if(ri(i)<=1./2432.**0.5)   sh(i)=1.+0.17*ri(i)*2432.**0.5
!  if(ri(i)>1./2432.**0.5)  sh(i)=1.+0.55*ri(i)**(2./3.)*2432.**(1./3.)
!  re(i)=(1.5*ar)**(1./3.)*ri(i)/1000.
!  uc2(i)=theta*(rho0-rhoi)*gv*2.*re(i)/rho0
!  vi(i)=pi*ri(i)**3*ar*2.
!  if(di<1.) uc2(i)=(3.724*(di/10.)**0.256/100.)**2
!  if(di>=1.)uc2(i)=(7.656*(di/10.)**0.569/100.)**2
!enddo

!do i=1, ni-1
!  dv(i)= vi(i+1)-vi(i)
!enddo

dzz=sn*dx 
do k=0, nx+1
  do j=0, ny+1
    dp(j,k)=ib-k*dzz
  enddo
enddo

rhopini=rho0*(1.+bs*(smix-s0)-bt*(tmix-t0))
rhoa=rho0*(1.+bs*(slow-s0)-bt*(tlow-t0))

do j=0,ny+1
    soli(j)=hoch
    isoli(j)=int(soli(j)/dz)
  if(j>=int(10.0e3/dy)+1 .and. j<=int(10.0e3/dy)+1+int(loch/dy))then
    soli(j)=hoch-(-0.5*hoch*cos(2.*pi/loch*dy*(j-int(10.0e3/dy)-1))+0.5*hoch)
    isoli(j)=int(soli(j)/dz)
  endif
enddo

aup=0.
ach=0.
adw=0.
js1=0
js2=0
js3=0
do k=0, nx+1 
  do j=0, ny+1
    dp(j,k)=dp(j,k)+soli(j)*cs     
  enddo
enddo

do k=1, nx-1
  do j=1, ny-1
    if(j<=int(10.0e3/dy))then   
      aup=aup+dx*dy   
      js1=js1+1
    elseif(j>=int(10.0e3/dy) .and. j<=int(10.0e3/dy)+int(loch/dy))then
      ach=ach+dx*dy
      js2=js2+1
    else    
      adw=adw+dx*dy 
      js3=js3+1
    endif       
  enddo
enddo

tdp=dz*nz
DO j = 1,ny 
DO k = 1,nx
  depth(j,k) = tdp-soli(j)  
END DO
depth(j,0)=depth(j,1)
depth(j,nx+1)=depth(j,nx)
END DO

DO k = 1,nx
  depth(0,k)=depth(1,k)
  depth(ny+1,k)=depth(ny,k)
END DO

open ( 1, file = 'depth.dat' )
write(1,*) "variables=x,y,dp"
write(1,*) "zone i=", nx, "j=", ny, "f=point"
DO j = 1, ny
  do k= 1, nx
    WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., dp(j,k)
  enddo
END DO  
close(1) 

wet = .true.
land= .false.

DO j = 1,ny
DO k = 1,nx
do i=0, isoli(j)
  wet(i,j,k) = .false.
  if(k==1) wet(i,j,k-1)=wet(i,j,k)
  if(k==nx)wet(i,j,k+1)=wet(i,j,k)
  if(j==1) wet(i,j-1,k)=wet(i,j,k)
  if(j==ny)wet(i,j+1,k)=wet(i,j,k)  
enddo 
END DO
END DO

!DO k = 1, nx
!  depth(0,k)=0.
  !depth(ny+1,k)=0.
  !land(0,k)= .true.
  !land(ny+1,k)= .true.
  !do i=0,nz+1
  !  wet(i,0,k) = .false.
    !wet(i,ny+1,k) = .false.
!  enddo
!END DO

DO j = 1, ny
  depth(j,0)=0.
  land(j,0)= .true.
  do i=0,nz+1
    wet(i,j,0) = .false.
  enddo
  
  !depth(j,nx+1)=0.
  !land(j,nx+1)= .true.
  !do i=0,nz+1
  !  wet(i,j,nx+1) = .false.
  !enddo
END DO

do j=int(8.0d3/dy)+1, int(8.0d3/dy)+1+int(pw/dy)
  do i=isoli(j)+1, isoli(j)+isdp 
    wet(i,j,0) = .true. 
    land(j,0) = .false.
  enddo  
enddo

! set initial arrays
DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
 p(i,j,k) = 0.0
 q(i,j,k) = 0.0
 dq(i,j,k) = 0.0
 u(i,j,k) = 0. 
 un(i,j,k) = 0. 
 ustar(i,j,k) = 0. 
 v(i,j,k) = 0.
 vn(i,j,k) = 0.
 vstar(i,j,k) = 0.
 w(i,j,k) = 0.0
 wn(i,j,k) = 0.0
 wstar(i,j,k) = 0.0
 uss(i,j,k) = 0.
 !do ll= 1, ni
 !  fc(ll,i,j,k) = 0. 
 !  fcn(ll,i,j,k) = 0.
 !enddo 
 !tfc(i,j,k)=0.
 az(i,j,k)=azmin 
 taz(i,j,k)=azmin2
 saz(i,j,k)=azmin2
 ek(i,j,k)=ekmin
 ekn(i,j,k)=ekmin
 ep(i,j,k)=epmin
 epn(i,j,k)=epmin
 ah(i,j,k) = 0.1
 s(i,j,k) = slow
 t(i,j,k) = tlow
 ssn(i,j,k)= s(i,j,k)
 tn(i,j,k)= t(i,j,k)
 rho(i,j,k)=rho0*(1.+bs*(s(i,j,k)-s0)-bt*(t(i,j,k)-t0)) 
 !rho(i,j,k)=rho(i,j,k)+tfc(i,j,k)*(rhoi-rho(i,j,k)) 
 rhoini(i,j,k)=rho(i,j,k)
END DO
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
 usum(j,k) = 0.0
 vsum(j,k) = 0.0
END DO
END DO

! grid parameters
drdx = 1.0/(rho0*dx)
drdy = 1.0/(rho0*dy)
drdz = 1.0/(rho0*dz)

! friction with coastlines
! slip = 0.0
! slip = 1.0
! slip = 2.0

slip = 0.0

DO k = 0,nx+1
DO j = 0,ny+1
DO i = 0,nz+1
  dry(i,j,k)=.NOT.wet(i,j,k)  
END DO
END DO
END DO

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
  at(i,j,k) = dx/dz
  ab(i,j,k) = dx/dz
  ae(i,j,k) = dz/dx
  aw(i,j,k) = dz/dx
  an(i,j,k) = dz*dx/(dy*dy)
  as(i,j,k) = dz*dx/(dy*dy)
  IF(dry(i,j,k+1)) ae(i,j,k) = 0.0
  IF(dry(i,j,k-1)) aw(i,j,k) = 0.0
  IF(dry(i,j+1,k)) an(i,j,k) = 0.0
  IF(dry(i,j-1,k)) as(i,j,k) = 0.0  
  IF(dry(i+1,j,k)) ab(i,j,k) = 0.0
  IF(dry(i-1,j,k)) at(i,j,k) = 0.0
  atot(i,j,k) = ab(i,j,k)+at(i,j,k)+ae(i,j,k)+aw(i,j,k)+an(i,j,k)+as(i,j,k)
END DO
END DO
END DO

azn=azmin
tmel=0.
tpre=0.
loci=0
loca=0.
avup=0.
avvp=0.
avcp=0.
avbp=0.
wbb=0.

zom=int(hoch/2./dz)+1
zob=int(hoch/dz)+int((tdp-hoch)*0.50/dz)+1
END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters
implicit none
real(8) :: pressx, pressy, pressz
real(8) :: perr, q1, q2, term1
INTEGER :: nsor, nstop
real(8) :: advx(0:nz+1,0:ny+1,0:nx+1), advy(0:nz+1,0:ny+1,0:nx+1), advz(0:nz+1,0:ny+1,0:nx+1)
real(8) :: div, divv(0:nz+1,0:ny+1,0:nx+1), div1, div2, div3, um, vm
real(8) :: dif1, dif2, dif3, dif4, difh, difz, diff, diffu, diffv, diffw
real(8) :: drho, speed, fm
real(8) :: aztt(0:nz+1,0:ny+1,0:nx+1), azbb(0:nz+1,0:ny+1,0:nx+1), ahee(0:nz+1,0:ny+1,0:nx+1), ahww(0:nz+1,0:ny+1,0:nx+1), ahnn(0:nz+1,0:ny+1,0:nx+1), ahss(0:nz+1,0:ny+1,0:nx+1)
real(8) :: azt, azb, ahe, ahw, ahn, ahs

real(8) :: dudz, dvdz, drodz, tpp, kkk, ppp, us0 
real(8) :: tat, tamt, tams, rt, rs, ta, tb, tc, sb1, sb2, sb, rtc, rsc, tt1
real(8) :: utop, ubot, vtop, vbot

real(8), parameter :: c1=0.10
real(8) :: c2, term2, term3
real(8) :: uNN, uS, vE, vW

c2 = c1*dx*dy

! south-north boundaries
uad=adu*uini
DO j = 1,ny
  isl=isoli(j) 
DO i = isl+1,nz 
  if(j>=int(8.0d3/dy)+1 .and. j<=int(8.0d3/dy)+1+int(pw/dy) .and. i<=isl+isdp)then
    if(i<=isl+isd)then
      if(i==isl+1)then 
        u(i,j,0)=uad
        t(i,j,0)=tlow-ad*(tlow-tmix)
        s(i,j,0)=slow-ad*(slow-smix)
        !do ll=1,ni
        !  fc(ll,i,j,0) = ad*fini/ni  
        !  fc(ll,i,j,1)=fc(ll,i,j,0)
        !enddo
      else
        u(i,j,0)=u(i-1,j,0)
        t(i,j,0)=t(i-1,j,0)
        s(i,j,0)=s(i-1,j,0)
        !do ll=1,ni
        !  fc(ll,i,j,0)=fc(ll,i-1,j,0)
        !  fc(ll,i,j,1)=fc(ll,i,j,0)
        !enddo        
      endif
      v(i,j,0)=v(i,j,1) 
      w(i,j,0)=w(i,j,1)        
    else
      u(i,j,0)=u(i,j,1)
      v(i,j,0)=v(i,j,1)
      w(i,j,0)=w(i,j,1)
      t(i,j,0) = t(i,j,1)
      s(i,j,0) = s(i,j,1)
      !do ll=1,ni
      !  fc(ll,i,j,0) = fc(ll,i,j,1)
      !enddo  
    endif 
  else
    u(i,j,0)=0.
    v(i,j,0)=0.
    w(i,j,0)=0.
    t(i,j,0) = t(i,j,1)
    s(i,j,0) = s(i,j,1)
    !do ll=1,ni
    !  fc(ll,i,j,0) = fc(ll,i,j,1)
    !enddo       
  endif
  !tfc(i,j,0)=0. 
  !do ll=1,ni
  !  tfc(i,j,0) = tfc(i,j,0)+fc(ll,i,j,0) 
  !enddo
  rho(i,j,0)=rho0*(1.+bs*(s(i,j,0)-s0)-bt*(t(i,j,0)-t0)) 
  !rho(i,j,0)=rho(i,j,0)+tfc(i,j,0)*(rhoi-rho(i,j,0)) 
  rho(i,j,1)=rho0*(1.+bs*(s(i,j,1)-s0)-bt*(t(i,j,1)-t0)) 
  !rho(i,j,1)=rho(i,j,1)+tfc(i,j,1)*(rhoi-rho(i,j,1))   
  ek(i,j,0) = ek(i,j,1)
  ep(i,j,0) = ep(i,j,1)
  az(i,j,0) = az(i,j,1)
  taz(i,j,0) = taz(i,j,1)
  saz(i,j,0) = saz(i,j,1)
  ah(i,j,0) = ah(i,j,1)

  ek(i,j,nx+1) = ek(i,j,nx)
  ep(i,j,nx+1) = ep(i,j,nx)
  az(i,j,nx+1) = az(i,j,nx)
  taz(i,j,nx+1) = taz(i,j,nx)
  saz(i,j,nx+1) = saz(i,j,nx)
  ah(i,j,nx+1) = ah(i,j,nx)
  u(i,j,nx+1)=u(i,j,nx)
  v(i,j,nx+1)=v(i,j,nx)
  w(i,j,nx+1)=w(i,j,nx)
  t(i,j,nx+1) = t(i,j,nx)
  s(i,j,nx+1) = s(i,j,nx)
  !tfc(i,j,nx+1)=0. 
  !do ll=1,ni
  !  fc(ll,i,j,nx+1) = fc(ll,i,j,nx) 
  !  tfc(i,j,nx+1)=tfc(i,j,nx+1)+fc(ll,i,j,nx+1) 
  !enddo 
  rho(i,j,nx+1)=rho0*(1.+bs*(s(i,j,nx+1)-s0)-bt*(t(i,j,nx+1)-t0)) 
  !rho(i,j,nx+1)=rho(i,j,nx+1)+tfc(i,j,nx+1)*(rhoi-rho(i,j,nx+1))  
END DO
END DO

do k=1,nx
do j=1,ny
  isl=isoli(j)    
! upper boundary  
  u(isl,j,k) = 0. 
  v(isl,j,k) = 0.
  sqcd=sqrt(0.0025) 
  um = u(isl+1,j,k)
  vm = 0.25*(v(isl+1,j,k)+v(isl+1,j,k+1)+ &
& v(isl+1,j-1,k)+v(isl+1,j-1,k+1))
  speed = SQRT(um*um+vm*vm)
  tbu(j,k) = sqcd**2*um*speed
  vm = v(isl+1,j,k)
  um = 0.25*(u(isl+1,j,k)+u(isl+1,j,k-1)+ &
& u(isl+1,j+1,k)+u(isl+1,j+1,k-1))  
  speed = SQRT(um*um+vm*vm)
  tbv(j,k) = sqcd**2*vm*speed  
  uss(isl,j,k)=(tbu(j,k)**2+tbv(j,k)**2)**0.25
  az(isl,j,k)=uss(isl,j,k)*sqcd*dz
  if(az(isl,j,k)==0.) az(isl,j,k)=azmin
  if(az(isl,j,k)>azmax) az(isl,j,k)=azmax
  taz(isl,j,k)=az(isl,j,k)
  if(az(isl,j,k)==azmin) taz(isl,j,k)=azmin2
  saz(isl,j,k)=taz(isl,j,k) 
  ek(isl,j,k)=uss(isl,j,k)**2/0.3
  if(ek(isl,j,k)==0.) ek(isl,j,k)=ekmin 
  ep(isl,j,k)=uss(isl,j,k)**3/sqcd/dz 
  if(ep(isl,j,k)==0.) ep(isl,j,k)=epmin 
  ah(isl,j,k)=ah(isl+1,j,k)

! three-equation formulation
  us0=uss(isl,j,k)
  if (us0>=uscr .and. adh==1.) then
    if(wbb(j,k)<0.)then
      etst=1.
    else
      etst=1./(1.+0.052*us0/abs(f)*0.4*wbb(j,k)/us0**3/0.2)**0.5        
    endif
    tat=1/0.4*log(us0**2*0.052* etst**2/5./abs(f)/1.95e-6)+1./2./0.052/etst-1./0.4  
    tamt=12.5*13.8**(2./3.)-6.
    tams=12.5*2432.**(2./3.)-6.
    if(1./(tat+tamt)>=0.0113)then
      rt=us0*0.0113
    else
      rt=us0/(tat+tamt)  
    endif
    if(1./(tat+tams)>=3.8d-4)then
      rs=us0*3.8d-4
    else
      rs=us0/(tat+tams)  
    endif       
    cii=ci/c0 
    lff=lf/c0
    kzt=(taz(isl,j,k)+taz(isl+1,j,k))*0.5
    kzs=(saz(isl,j,k)+saz(isl+1,j,k))*0.5
    sanj=rt*kzt/dz*t(isl+1,j,k)/(rt+kzt/dz)
    wujx=rt**2/(rt+kzt/dz)-rt
    shiz=kzs/dz*s(isl+1,j,k)/(kzs/dz+rs)
    yuel=l2+l3*dp(j,k)-ti
    lingx=kzs/dz+rs
    ta=l1*cii*rs**2/lingx-l1*rs*cii-l1*wujx
    tb=rs**2/lingx*cii*yuel+rs*shiz*l1*cii-rs*lff-rs*cii*yuel+rs**2/lingx*lff-sanj-(l2+l3*dp(j,k))*wujx
    tc=rs*shiz*lff+rs*shiz*cii*yuel
    sb1=(-tb+sqrt(tb**2-4.*ta*tc))/2./ta
    sb2=(-tb-sqrt(tb**2-4.*ta*tc))/2./ta
    if(sb1>0. .and. sb1<slow) sb=sb1
    if(sb2>0. .and. sb2<slow) sb=sb2
    if(sb/=sb1 .and. sb/=sb2)then
      s(isl,j,k)=s(isl+1,j,k) 
      t(isl,j,k)=t(isl+1,j,k) 
      w(isl,j,k)=0.
      tb=-100.
      sb=50.
      ttbb(j,k)=tb
      ssbb(j,k)=sb 
      mel(j,k)=0. 
      wbb(j,k)=0.
      goto 444
    else
      tb=l1*sb+l2+l3*dp(j,k)
      w(isl,j,k)=-rs/sb*(rs*sb+kzs*s(isl+1,j,k)/dz)/(kzs/dz+rs)+rs   
    endif
    t(isl,j,k)=kzt/dz*t(isl+1,j,k)/(rt+kzt/dz)+rt/(rt+kzt/dz)*tb
    s(isl,j,k)=(rs*sb+kzs/dz*s(isl+1,j,k))/(kzs/dz+rs)    
    ttbb(j,k)=tb
    ssbb(j,k)=sb
    tmel(j,k)=tmel(j,k)-w(isl,j,k)*dt*rho0/rhoi
    mel(j,k)=-w(isl,j,k)*3600.*24.*365.*rho0/rhoi   
  else 
    s(isl,j,k)=s(isl+1,j,k) 
    t(isl,j,k)=t(isl+1,j,k) 
    w(isl,j,k)=0.
    tb=-100.
    sb=50.
    ttbb(j,k)=tb
    ssbb(j,k)=sb 
    mel(j,k)=0.  
    wbb(j,k)=0.
  endif
444  tfb(j,k)=t(isl,j,k)-(l1*s(isl,j,k)+l2+l3*dp(j,k))
  
  !tfc(isl,j,k)=0.  
  !atop = 0.5*(az(isl,j,k)+az(isl+1,j,k))
  !do ll=1, ni  !precipitation for shear form
  !  usuc=us0**2/uc2(ll)
  !  if(usuc<1.) then
  !    fc(ll,isl,j,k)=fc(ll,isl+1,j,k)/(1.-wi(ll)*usuc*dz/atop)
  !    pre(ll,j,k)=-wi(ll)*fc(ll,isl,j,k)*(1.-usuc)
  !  else 
  !    fc(ll,isl,j,k)=fc(ll,isl+1,j,k)/(1.-wi(ll)*dz/atop) 
  !    pre(ll,j,k)=0.
  !  endif
  !  tfc(isl,j,k)=tfc(isl,j,k)+fc(ll,isl,j,k)    
  !  tpre(j,k)=tpre(j,k)-pre(ll,j,k)*dt
  !enddo
  rho(isl,j,k)=rho0*(1.+bs*(s(isl,j,k)-s0)-bt*(t(isl,j,k)-t0)) 
  !rho(isl,j,k)=rho(isl,j,k)+tfc(isl,j,k)*(rhoi-rho(isl,j,k))   
  
! lower boundary
  ek(nz+1,j,k)=ek(nz,j,k)
  ep(nz+1,j,k)=ep(nz,j,k)  
  az(nz+1,j,k)=az(nz,j,k) 
  taz(nz+1,j,k)=taz(nz,j,k)
  saz(nz+1,j,k)=saz(nz,j,k)
  ah(nz+1,j,k)=ah(nz,j,k)
  u(nz+1,j,k)=u(nz,j,k)
  v(nz+1,j,k)=v(nz,j,k)
  w(nz+1,j,k)=w(nz,j,k) 
  t(nz,j,k) = tlow
  s(nz,j,k) = slow  
  t(nz+1,j,k) = t(nz,j,k)
  s(nz+1,j,k) = s(nz,j,k)
  !tfc(nz+1,j,k)=0.  
  !do ll=1, ni
  !  fc(ll,nz+1,j,k)= fc(ll,nz,j,k)
  !  tfc(nz+1,j,k)=tfc(nz+1,j,k)+fc(ll,nz+1,j,k)    
  !enddo 
  rho(nz+1,j,k)=rho0*(1.+bs*(s(nz+1,j,k)-s0)-bt*(t(nz+1,j,k)-t0)) 
  !rho(nz+1,j,k)=rho(nz+1,j,k)+tfc(nz+1,j,k)*(rhoi-rho(nz+1,j,k)) 
enddo
enddo

! east-west boundaries
DO k = 1,nx
DO i = 1,nz
 ek(i,0,k) = ek(i,ny,k)
 ep(i,0,k) = ep(i,ny,k)
 az(i,0,k) = az(i,ny,k)
 taz(i,0,k) = taz(i,ny,k)
 saz(i,0,k) = saz(i,ny,k)
 ah(i,0,k) = ah(i,ny,k)
 u(i,0,k)=u(i,ny,k)
 v(i,0,k)=v(i,ny,k)
 w(i,0,k)=w(i,ny,k)
 t(i,0,k) = t(i,ny,k)
 s(i,0,k) = s(i,ny,k)
 !tfc(i,0,k)=0. 
 !do ll=1,ni
 !  fc(ll,i,0,k) = fc(ll,i,1,k)
 !  tfc(i,0,k) = tfc(i,0,k)+fc(ll,i,0,k) 
 !enddo
 rho(i,0,k)=rho0*(1.+bs*(s(i,0,k)-s0)-bt*(t(i,0,k)-t0)) 
 !rho(i,0,k)=rho(i,0,k)+tfc(i,0,k)*(rhoi-rho(i,0,k))  
 
 ek(i,ny+1,k) = ek(i,1,k)
 ep(i,ny+1,k) = ep(i,1,k)
 az(i,ny+1,k) = az(i,1,k)
 taz(i,ny+1,k) = taz(i,1,k)
 saz(i,ny+1,k) = saz(i,1,k)
 ah(i,ny+1,k) = ah(i,1,k)
 u(i,ny+1,k)=u(i,1,k)!0.
 v(i,ny+1,k)=v(i,1,k)!0.
 w(i,ny+1,k)=w(i,1,k)!0.
 t(i,ny+1,k) = t(i,1,k)
 s(i,ny+1,k) = s(i,1,k)
 !tfc(i,ny+1,k)=0.  
 !do ll=1,ni
 !  fc(ll,i,ny+1,k) = fc(ll,i,ny,k)
 !  tfc(i,ny+1,k) = tfc(i,ny+1,k)+fc(ll,i,ny+1,k)  
 !enddo 
 rho(i,ny+1,k)=rho0*(1.+bs*(s(i,ny+1,k)-s0)-bt*(t(i,ny+1,k)-t0)) 
 !rho(i,ny+1,k)=rho(i,ny+1,k)+tfc(i,ny+1,k)*(rhoi-rho(i,ny+1,k))   
END DO
END DO

IF(MOD(n,nout)==0)THEN

  ntem=n/nout 
  write( cTemp,'(i3)' ) ntem   
  open ( 1, file = 'us' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,us"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  uss(isoli(j),j,k)
    enddo
  END DO  
  close(1)
  
  open ( 1, file = 'ut' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,u"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., u(isoli(j)+1,j,k)
    enddo
  END DO  
  close(1)
  
  open ( 1, file = 'vt' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,v"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., v(isoli(j)+1,j,k)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'uvt' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,u,v"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., u(isoli(j)+1,j,k), v(isoli(j)+1,j,k)
    enddo
  END DO  
  close(1)      

  open ( 1, file = 'uvp' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,up,vp"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., avup(j,k), avvp(j,k)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'wt' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,w"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., w(isoli(j)+1,j,k)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'um' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,u"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., u(zom,j,k)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vm' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,v"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., v(zom,j,k)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wm' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,w"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., w(zom,j,k)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'ub' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,u"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., u(zob,j,k)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vb' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,v"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., v(zob,j,k)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wb' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,w"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., w(zob,j,k)
    enddo
  END DO  
  close(1)  

  open ( 1, file = 'azx1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,az"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, az(i,j,x1)
    enddo
  END DO  
  close(1)   
  open ( 1, file = 'azx2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,az"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, az(i,j,x2)
    enddo
  END DO  
  close(1)  
  open ( 1, file = 'azx3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,az"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, az(i,j,x3)
    enddo
  END DO  
  close(1)  
  open ( 1, file = 'azx4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,az"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, az(i,j,x4)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'azx5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,az"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, az(i,j,x5)
    enddo
  END DO  
  close(1)  
  
  open ( 1, file = 'ux1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,u"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, u(i,j,x1)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vx1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,v"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, v(i,j,x1)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wx1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,w"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, w(i,j,x1)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'ux2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,u"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, u(i,j,x2)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vx2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,v"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, v(i,j,x2)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wx2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,w"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, w(i,j,x2)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'ux3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,u"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, u(i,j,x3)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vx3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,v"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, v(i,j,x3)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wx3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,w"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, w(i,j,x3)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'ux4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,u"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, u(i,j,x4)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vx4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,v"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, v(i,j,x4)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wx4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,w"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, w(i,j,x4)
    enddo
  END DO  
  close(1)      

  open ( 1, file = 'ux5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,u"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, u(i,j,x5)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'vx5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,v"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, v(i,j,x5)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'wx5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,w"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, w(i,j,x5)
    enddo
  END DO  
  close(1)   
  
  open ( 1, file = 'tx1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,t"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, t(i,j,x1)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'sx1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,s"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, s(i,j,x1)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'buox1-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,buo"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, (rhoa-rho(i,j,x1))/rho0*gv*cs
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'tx2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,t"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, t(i,j,x2)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'sx2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,s"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, s(i,j,x2)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'buox2-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,buo"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, (rhoa-rho(i,j,x2))/rho0*gv*cs
    enddo
  END DO  
  close(1)     
  
  open ( 1, file = 'tx3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,t"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, t(i,j,x3)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'sx3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,s"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, s(i,j,x3)
    enddo
  END DO  
  close(1)     

  open ( 1, file = 'buox3-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,buo"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, (rhoa-rho(i,j,x3))/rho0*gv*cs
    enddo
  END DO  
  close(1)     
  
  open ( 1, file = 'tx4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,t"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, t(i,j,x4)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'sx4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,s"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, s(i,j,x4)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'buox4-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,buo"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, (rhoa-rho(i,j,x4))/rho0*gv*cs
    enddo
  END DO  
  close(1)     

  open ( 1, file = 'tx5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,t"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, t(i,j,x5)
    enddo
  END DO  
  close(1) 
  
  open ( 1, file = 'sx5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,s"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, s(i,j,x5)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'buox5-' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=y,z,buo"
  write(1,*) "zone i=", ny, "j=", nz+1, "f=point"
  DO i = 0, nz
    do j= 1, ny
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (j-1)*dy/1000., -i*dz, (rhoa-rho(i,j,x5))/rho0*gv*cs
    enddo
  END DO  
  close(1)   
  
  open ( 1, file = 'tscb' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,tscb"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., t(isoli(j)+1,j,k)-(l1*s(isoli(j)+1,j,k)+l2+l3*(dp(j,k)+dz*cs)) !tfb(j,k)  (dp(j,k)+dz*cs)
    enddo
  END DO  
  close(1)  

  open ( 1, file = 'tt' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,t"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  t(isoli(j)+1,j,k)
    enddo
  END DO  
  close(1)  
  
  open ( 1, file = 'st' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,s"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  s(isoli(j)+1,j,k)
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'buot' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,buo"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., (rhoa-rho(isoli(j)+1,j,k))/rho0*gv*cs
    enddo
  END DO  
  close(1)    

  open ( 1, file = 'buop' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,buop"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., avbp(j,k)
    enddo
  END DO  
  close(1)   
  
  !open ( 1, file = 'tfct' // trim(adjustl( cTemp )) // '.dat' )
  !write(1,*) "variables=x,y,tfc"
  !write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !DO j = 1, ny
  !  do k= 1, nx
  !    WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., tfc(isoli(j)+1,j,k)
  !  enddo
  !END DO  
  !close(1)     

  !open ( 1, file = 'tscm' // trim(adjustl( cTemp )) // '.dat' )
  !write(1,*) "variables=x,y,tsc"
  !write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !DO j = 1, ny
  !  do k= 1, nx
  !    WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  t(zom,j,k)-(l1*s(zom,j,k)+l2+l3*(dp(j,k)+(zom-isoli(j))*dz*cs))
  !  enddo
  !END DO  
  !close(1)  

  open ( 1, file = 'tm' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,t"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  t(zom,j,k)
    enddo
  END DO  
  close(1)  
  
  open ( 1, file = 'sm' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,s"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  s(zom,j,k)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'buom' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,buo"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., (rhoa-rho(zom,j,k))/rho0*gv*cs
    enddo
  END DO  
  close(1)    
  
  !open ( 1, file = 'tfcm' // trim(adjustl( cTemp )) // '.dat' )
  !write(1,*) "variables=x,y,tfc"
  !write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !DO j = 1, ny
  !  do k= 1, nx
  !    WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., tfc(zom,j,k)
  !  enddo
  !END DO  
  !close(1)    
  
  open ( 1, file = 'tb' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,t"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  t(zob,j,k)
    enddo
  END DO  
  close(1)  
  
  open ( 1, file = 'sb' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,s"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000.,  s(zob,j,k)
    enddo
  END DO  
  close(1)    
  
  open ( 1, file = 'buob' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,buo"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., (rhoa-rho(zob,j,k))/rho0*gv*cs
    enddo
  END DO  
  close(1)    
  
  !do ll=1, ni  
  !  ncl=ll
  !  write( ccl,'(i3)' ) ncl
  !  open ( 1, file = trim(adjustl( ccl )) // 'Ct' // trim(adjustl( cTemp )) // '.dat' )
  !  write(1,*) "variables=x,y,"//trim(adjustl( ccl ))//"C"
  !  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !  DO j = 1, ny
  !    DO k = 1, nx
  !      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., fc(ll,isoli(j)+1,j,k)
  !    ENDDO
  !  END DO
  !  close(1)
  !
  !  open ( 1, file = trim(adjustl( ccl )) // 'Cm' // trim(adjustl( cTemp )) // '.dat' )
  !  write(1,*) "variables=x,y,"//trim(adjustl( ccl ))//"C"
  !  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !  DO j = 1, ny
  !    DO k = 1, nx
  !      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., fc(ll,zom,j,k)
  !    ENDDO
  !  END DO
  !  close(1)
  !enddo

  open ( 1, file = 'mel' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,mel"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., mel(j,k)
    enddo
  END DO  
  close(1)   
  
  open ( 1, file = 'tmel' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,tmel"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., tmel(j,k)
    enddo
  END DO  
  close(1)       
  
  !open ( 1, file = 'pre' // trim(adjustl( cTemp )) // '.dat' )
  !write(1,*) "variables=x,y,pre"
  !write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !DO j = 1, ny
  !  do k= 1, nx
  !    WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., tpre(j,k)
  !  enddo
  !END DO  
  !close(1)   
  
  !open ( 1, file = 'hsc' // trim(adjustl( cTemp )) // '.dat' )
  !write(1,*) "variables=x,y,hsc"
  !write(1,*) "zone i=", nx, "j=", ny, "f=point"
  !DO j = 1, ny
  !  do k= 1, nx
  !    WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., hsc(j,k)
  !  enddo
  !END DO  
  !close(1)   

  open ( 1, file = 'hplm' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,y,hplm"
  write(1,*) "zone i=", nx, "j=", ny, "f=point"
  DO j = 1, ny
    do k= 1, nx
      WRITE(1,'(f8.2,x,f8.2,x,f20.15)') (k-1)*dx/1000., (j-1)*dy/1000., loca(j,k)
    enddo
  END DO  
  close(1)     

  open ( 1, file = 'buocc' // trim(adjustl( cTemp )) // '.dat' )
  write(1,*) "variables=x,z,buo"
  write(1,*) "zone i=", nx, "j=", nz+1, "f=point"
  DO i = 0, nz
    do k= 1, nx
      WRITE(1,'(f8.2,x,f10.4,x,f20.15)') (k-1)*dx/1000., -i*dz, (rhoa-rho(i,int(14000./dy)+1,k))/rho0*gv*cs !22500.
    enddo
  END DO  
  close(1)     
ENDIF

IF(MOD(n,nout2)==0)THEN
  ntem=n/nout2   
  WRITE(10,'(f10.4,x,f20.15,x,f20.15,x,f20.15,x,f20.15,x,f20.15,x,f20.15,x,f20.15,x,f20.15)') ntem*ntt2, upmel, chmel, dwmel, smel, aupmel, achmel, adwmel, asmel 
  if(optf==.true.) WRITE(11,'(f10.4,x,f20.15,x,f20.15,x,f20.15,x,f20.15,x,f20.15)') ntem*ntt2, affl, atfl*1.d-9, abfl*1.d6, atq*1.d-3, abq
ENDIF

DO k = 0,nx+1
DO j = 0,ny+1
DO i = 0,nz+1
  CuP(i,j,k) = 0.5*(u(i,j,k)+abs(u(i,j,k)))*dt/dx
  CuN(i,j,k) = 0.5*(u(i,j,k)-abs(u(i,j,k)))*dt/dx
  CvP(i,j,k) = 0.5*(v(i,j,k)+abs(v(i,j,k)))*dt/dy
  CvN(i,j,k) = 0.5*(v(i,j,k)-abs(v(i,j,k)))*dt/dy
  CwP(i,j,k) = 0.5*(w(i,j,k)+abs(w(i,j,k)))*dt/dz
  CwN(i,j,k) = 0.5*(w(i,j,k)-abs(w(i,j,k)))*dt/dz
END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
 divv(i,j,k) = (u(i,j,k)-u(i,j,k-1))/dx+(v(i,j,k)-v(i,j-1,k))/dy+(w(i,j,k)-w(i+1,j,k))/dz
END DO
END DO
END DO

! calculation of turbulence kinetic energy and its dissipation rate 
DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = ek(i,j,k)
END DO
END DO
END DO

CALL advect

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
if(wet(i,j,k))then
  div = dt*B(i,j,k)*divv(i,j,k)
  aztt(i,j,k) = 0.5*(az(i,j,k)+az(i-1,j,k))
  dif1 = aztt(i,j,k)*(ek(i-1,j,k)-ek(i,j,k))/dz
  azbb(i,j,k) = 0.5*(az(i,j,k)+az(i+1,j,k))
  dif2 = azbb(i,j,k)*(ek(i,j,k)-ek(i+1,j,k))/dz
  if(i==nz) dif2=0.
  difz = (dif1-dif2)/dz
  diff = dt*difz
  utop = 0.5*(u(i-1,j,k)+u(i-1,j,k-1))
  ubot = 0.5*(u(i+1,j,k)+u(i+1,j,k-1))
  dudz = (utop-ubot)/(2.0*dz)
  vtop = 0.5*(v(i-1,j,k)+v(i-1,j-1,k))
  vbot = 0.5*(v(i+1,j,k)+v(i+1,j-1,k))
  dvdz = (vtop-vbot)/(2.0*dz)
  ps(i,j,k) = (dudz**2+dvdz**2)*az(i,j,k)
  pb(i,j,k) = gv*cs*(bs*saz(i,j,k)*(s(i-1,j,k)-s(i+1,j,k))/(2.*dz)-bt*taz(i,j,k)*(t(i-1,j,k)-t(i+1,j,k))/(2.*dz))
  kkk=ek(i,j,k)
  ekn(i,j,k)=ek(i,j,k)+BN(i,j,k)+div+diff/1.4+dt*(ps(i,j,k)+pb(i,j,k)-ep(i,j,k))
  if(ekn(i,j,k)<0.) ekn(i,j,k)=kkk
 
  ! horizontal eddy viscosity/diffusivity--Smagorinsky, 1963
  term1 = ((u(i,j,k)-u(i,j,k-1))/dx)**2
  term2 = ((v(i,j,k)-v(i,j-1,k))/dy)**2
  uNN = 0.25*(u(i,j,k)+u(i,j,k-1)+u(i,j+1,k)+u(i,j+1,k-1))
  uS = 0.25*(u(i,j,k)+u(i,j,k-1)+u(i,j-1,k)+u(i,j-1,k-1))
  vE = 0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k+1)+v(i,j-1,k+1))
  vW = 0.25*(v(i,j,k)+v(i,j-1,k)+v(i,j,k-1)+v(i,j-1,k-1))
  term3 = 0.5*((uNN-uS)/dy+(vE-vW)/dx)**2  
  ah(i,j,k) = c2*SQRT(term1+term2+term3)!+1.0
endif
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = ep(i,j,k)
END DO
END DO
END DO

CALL advect

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
if(wet(i,j,k))then
  div = dt*B(i,j,k)*divv(i,j,k)
  dif1 = aztt(i,j,k)*(ep(i-1,j,k)-ep(i,j,k))/dz
  dif2 = azbb(i,j,k)*(ep(i,j,k)-ep(i+1,j,k))/dz
  if(i==nz) dif2=0.
  difz = (dif1-dif2)/dz
  diff = dt*difz
  ppp=ep(i,j,k)
  epn(i,j,k)=ep(i,j,k)+BN(i,j,k)+div+diff/1.3+dt*ep(i,j,k)/ek(i,j,k)*(1.44*ps(i,j,k)+0.8*pb(i,j,k)-1.92*ep(i,j,k))
  if(epn(i,j,k)<0.) epn(i,j,k)=ppp
endif
END DO
END DO
END DO

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
if(wet(i,j,k))then  
  ek(i,j,k)=ekn(i,j,k)
  ep(i,j,k)=epn(i,j,k)
endif    
END DO
END DO
END DO

DO i = 1, nz
DO j = 1, ny
DO k = 1, nx
if(wet(i,j,k))then
  az(i,j,k)=0.09*ek(i,j,k)**2/ep(i,j,k)
  if(az(i,j,k)< azmin) az(i,j,k)=azmin
  if(az(i,j,k)> azmax) az(i,j,k)=azmax  
  if(ek(i,j,k)<=ekmin) az(i,j,k)=azmin  
endif
END DO
END DO 
END DO

DO i = 1, nz
DO j = 1, ny
DO k = 1, nx
if(wet(i,j,k))then
  n22(i,j,k)=-gv*cs/rho0*(rho(i-1,j,k)-rho(i+1,j,k))/(2.*dz)
  utop = 0.5*(u(i-1,j,k)+u(i-1,j,k-1))
  ubot = 0.5*(u(i+1,j,k)+u(i+1,j,k-1))
  dudz = (utop-ubot)/(2.0*dz)
  vtop = 0.5*(v(i-1,j,k)+v(i-1,j-1,k))
  vbot = 0.5*(v(i+1,j,k)+v(i+1,j-1,k))
  dvdz = (vtop-vbot)/(2.0*dz)
  duvz = dudz**2+dvdz**2
  if(az(i,j,k)==azmin)then
    taz(i,j,k)=azmin2
    saz(i,j,k)=azmin2
  else
    if(duvz==0.) then
      taz(i,j,k)=azmin2
    else
      rch=n22(i,j,k)/duvz 
      if(rch<=0.)then
        taz(i,j,k)= az(i,j,k)/0.7 
      else
        !taz(i,j,k)=az(i,j,k)/(0.7*exp(-3./0.7*rch)+4.*rch)  ! Venayagamoorthy and Stretch, 2010
        taz(i,j,k)=az(i,j,k)/(0.7*(1.+10.*rch)**(-0.5)/(1.+10./3.*rch)**(-1.5)) ! Munk and Anderson, 1948
      endif
      if(taz(i,j,k)<azmin2) taz(i,j,k)=azmin2
   endif
   saz(i,j,k)=taz(i,j,k)
  endif
endif
END DO
END DO 
END DO

! multiple-sized frazil ice module to be developed soon

if (cmin==0.) goto 888 

! calculation of T and S

888 DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = t(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
if(wet(i,j,k))then
  div = dt*B(i,j,k)*divv(i,j,k)
  dif1 = ahee(i,j,k)/prh*(t(i,j,k+1)-t(i,j,k))/dx
  IF(dry(i,j,k+1)) dif1 = 0.0
  dif2 = ahww(i,j,k)/prh*(t(i,j,k)-t(i,j,k-1))/dx
  IF(dry(i,j,k-1)) dif2 = 0.0
  dif3 = ahnn(i,j,k)/prh*(t(i,j+1,k)-t(i,j,k))/dy
  IF(dry(i,j+1,k)) dif3 = 0.0
  dif4 = ahss(i,j,k)/prh*(t(i,j,k)-t(i,j-1,k))/dy
  IF(dry(i,j-1,k)) dif4 = 0.0
  difh = (dif1-dif2)/dx + (dif3-dif4)/dy 
  aztt(i,j,k) = 0.5*(taz(i,j,k)+taz(i-1,j,k))
  if(i==isoli(j)+1)then
    if(mel(j,k)==0.)then
      dif1 = 0.
    else
      dif1 = aztt(i,j,k)*(t(i-1,j,k)-t(i,j,k))/dz 
    endif   
  else
    dif1 = aztt(i,j,k)*(t(i-1,j,k)-t(i,j,k))/dz 
  endif
  azbb(i,j,k) = 0.5*(taz(i,j,k)+taz(i+1,j,k))       
  dif2 = azbb(i,j,k)*(t(i,j,k)-t(i+1,j,k))/dz 
  if(i==nz) dif2=0.
  difz = (dif1-dif2)/dz
  diff = dt*(difh+difz)  
  tn(i,j,k)= t(i,j,k)+BN(i,j,k)+div+diff! +ttfc(i,j,k)*dt
endif
END DO
END DO
END DO

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
if(wet(i,j,k))then
  t(i,j,k)=tn(i,j,k)
endif    
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = s(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
if(wet(i,j,k))then
  div = dt*B(i,j,k)*divv(i,j,k)
  dif1 = ahee(i,j,k)/prh*(s(i,j,k+1)-s(i,j,k))/dx
  IF(dry(i,j,k+1)) dif1 = 0.0
  dif2 = ahww(i,j,k)/prh*(s(i,j,k)-s(i,j,k-1))/dx
  IF(dry(i,j,k-1)) dif2 = 0.0
  dif3 = ahnn(i,j,k)/prh*(s(i,j+1,k)-s(i,j,k))/dy
  IF(dry(i,j+1,k)) dif3 = 0.0
  dif4 = ahss(i,j,k)/prh*(s(i,j,k)-s(i,j-1,k))/dy
  IF(dry(i,j-1,k)) dif4 = 0.0
  difh = (dif1-dif2)/dx + (dif3-dif4)/dy 
  aztt(i,j,k) = 0.5*(saz(i,j,k)+saz(i-1,j,k))
  if(i==isoli(j)+1)then
    if(mel(j,k)==0.)then
      dif1 = 0.
    else
      dif1 = aztt(i,j,k)*(s(i-1,j,k)-s(i,j,k))/dz 
    endif   
  else
    dif1 = aztt(i,j,k)*(s(i-1,j,k)-s(i,j,k))/dz 
  endif
  azbb(i,j,k) = 0.5*(saz(i,j,k)+saz(i+1,j,k))   
  dif2 = azbb(i,j,k)*(s(i,j,k)-s(i+1,j,k))/dz
  if(i==nz) dif2=0.
  difz = (dif1-dif2)/dz
  diff = dt*(difh+difz)  
  ssn(i,j,k)= s(i,j,k)+BN(i,j,k)+div+diff! +sfc(i,j,k)*dt
endif
END DO
END DO
END DO

DO k = 1,nx
DO j = 1,ny
DO i = 1,nz
if(wet(i,j,k))then
  s(i,j,k)=ssn(i,j,k)
endif    
END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
if(wet(i,j,k))then
  rho(i,j,k) = rho0*(1.+bs*(s(i,j,k)-s0)-bt*(t(i,j,k)-t0)) 
  !rho(i,j,k) = rho(i,j,k)+tfc(i,j,k)*(rhoi-rho(i,j,k)) 
endif  
END DO
END DO
END DO

!DO k = 1, nx
!DO j = 1, ny
!  tt1=tfb(j,k)
!  if(tt1<0.) then
!    do i=isoli(j)+1,nz  
!      tt1=t(i,j,k)-(l1*s(i,j,k)+l2+l3*(dp(j,k)+(i-isoli(j))*dz*cs))
!      if(tt1>=0.) then 
!        hsc(j,k)=(i-isoli(j)-1)*dz
!        exit
!      endif
!    enddo
!  else
!    hsc(j,k)=0.
!  endif  
!ENDDO
!ENDDO

! calculaiton of plume thickness
        loca=0.
        loci=0
        upmel=0.
        chmel=0.
        dwmel=0.      
        fflux=0.
        tflux=0.
        bflux=0.        
        DO k=1, nx
        DO j=1, ny
            if((rhoa-rho(isoli(j)+1,j,k))>=0.15*(rhoa-rhopini)) then            
            jif1=0.
            jif2=0.
            ujif=0.
            vjif=0.
            DO i = isoli(j)+1, nz    
             wei(i,j,k)= (rhoa-0.5*(rho(i,j,k)+rho(i-1,j,k)))/rho0*gv*cs   
             if(i==isoli(j)+1) wei(i,j,k)= (rhoa-rho(i,j,k))/rho0*gv*cs   
             jif1=jif1+wei(i,j,k)*dz
             jif2=jif2+wei(i,j,k)*(i-isoli(j))*dz**2             
            enddo    
            loca(j,k)=jif2*2./jif1
            loci(j,k)=int(loca(j,k)/dz)+1
            endif               
            if(loca(j,k)/=0.)then
              DO i = isoli(j)+1, isoli(j)+loci(j,k) !nz      
                ujif=ujif+u(i,j,k)*dz
                vjif=vjif+v(i,j,k)*dz
              enddo
              avup(j,k)=ujif/loca(j,k)
              avvp(j,k)=vjif/loca(j,k)
              avcp(j,k)=sqrt(avup(j,k)**2+avvp(j,k)**2)
              avbp(j,k)=jif1/loca(j,k)
            else
              avup(j,k)=0.
              avvp(j,k)=0.
              avbp(j,k)=0.
            endif    
        END DO    
        END DO
        
do k=1, nx-1
  do j=1, ny-1
    if(j<=int(10.0e3/dy))then   
      upmel=upmel+tmel(j,k) 
    elseif(j>=int(10.0e3/dy) .and. j<=int(10.0e3/dy)+int(loch/dy))then
      chmel=chmel+tmel(j,k)
    else    
      dwmel=dwmel+tmel(j,k) 
    endif       
  enddo
enddo        

optf=.false. 
affl=0.
atfl=0.
abfl=0.
do k=nx, nx
  do j=int(10.0e3/dy), int(10.0e3/dy)+int(loch/dy)
      if(loci(j,k)/=0)then
      do i=isoli(j)+1, isoli(j)+loci(j,k)  
        fflux(k)=fflux(k)+u(i,j,k)
        tflux(k)=tflux(k)+u(i,j,k)*(t(i,j,k)+1.865)*rho0*c0*1.d-3
        bflux(k)=bflux(k)+u(i,j,k)*(rhoa-rho(i,j,k))/rho0*gv*cs
      enddo
      endif     
  enddo
  if(fflux(k)/=0.)then
    optf=.true.  
    affl=affl+fflux(k)*dy*dz*1.d-6
    atfl=atfl+tflux(k)*dy*dz
    abfl=abfl+bflux(k)*dy*dz*1.d-6
  endif
enddo
if(optf==.true.)then
  atq=atfl/(affl*1.d6)  
  abq=abfl/affl
endif

aupmel=upmel/js1  
achmel=chmel/js2
adwmel=dwmel/js3
upmel=upmel*dx*dy*rhoi*1.d-12 
chmel=chmel*dx*dy*rhoi*1.d-12
dwmel=dwmel*dx*dy*rhoi*1.d-12
smel=upmel+chmel+dwmel 
asmel=smel/(aup+ach+adw)/rhoi*1.d12       
        
! calculation of momentum equations

! calculate hydrostatic pressure
DO j = 0,ny+1
  isl=isoli(j)    
DO k = 0,nx+1
  P(0,j,k) = 0. 
  DO i = 1,nz+1
    if(i==isl)then
      P(i,j,k) = P(i-1,j,k) + rho(i-1,j,k)*gv*dz*cs
    elseif(i==isl+1)then
      P(i,j,k) = P(i-1,j,k) + rho(i,j,k)*gv*dz*cs
    else
      P(i,j,k) = P(i-1,j,k) + 0.5*(rho(i-1,j,k)+rho(i,j,k))*gv*dz*cs
    endif
  END DO
END DO
END DO

! calculate the nonlinear terms for u-momentum equation
DO i = 0,nz+1
DO j = 0,ny
DO k = 0,nx
  CuP(i,j,k) = 0.25*(u(i,j,k)+u(i,j,k+1)+abs(u(i,j,k))+abs(u(i,j,k+1)))*dt/dx
  CuN(i,j,k) = 0.25*(u(i,j,k)+u(i,j,k+1)-abs(u(i,j,k))-abs(u(i,j,k+1)))*dt/dx
  CvP(i,j,k) = 0.25*(v(i,j,k)+v(i,j,k+1)+abs(v(i,j,k))+abs(v(i,j,k+1)))*dt/dy
  CvN(i,j,k) = 0.25*(v(i,j,k)+v(i,j,k+1)-abs(v(i,j,k))-abs(v(i,j,k+1)))*dt/dy
  CwP(i,j,k) = 0.25*(w(i,j,k)+w(i,j,k+1)+abs(w(i,j,k))+abs(w(i,j,k+1)))*dt/dz
  CwN(i,j,k) = 0.25*(w(i,j,k)+w(i,j,k+1)-abs(w(i,j,k))-abs(w(i,j,k+1)))*dt/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = u(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k+1)-u(i,j,k-1))/dx
  div2 = 0.5*(v(i,j,k)+v(i,j,k+1)-v(i,j-1,k)-v(i,j-1,k+1))/dy
  div3 = 0.5*(w(i,j,k)+w(i,j,k+1)-w(i+1,j,k)-w(i+1,j,k+1))/dz
  div = dt*B(i,j,k)*(div1+div2+div3)  
  advx(i,j,k)= BN(i,j,k)+div
END DO
END DO
END DO

! calculate the nonlinear terms for v-momentum equation
DO i = 0,nz+1
DO j = 0,ny
DO k = 0,nx
  CuP(i,j,k) = 0.25*(u(i,j,k)+u(i,j+1,k)+abs(u(i,j,k))+abs(u(i,j+1,k)))*dt/dx
  CuN(i,j,k) = 0.25*(u(i,j,k)+u(i,j+1,k)-abs(u(i,j,k))-abs(u(i,j+1,k)))*dt/dx
  CvP(i,j,k) = 0.25*(v(i,j,k)+v(i,j+1,k)+abs(v(i,j,k))+abs(v(i,j+1,k)))*dt/dy
  CvN(i,j,k) = 0.25*(v(i,j,k)+v(i,j+1,k)-abs(v(i,j,k))-abs(v(i,j+1,k)))*dt/dy
  CwP(i,j,k) = 0.25*(w(i,j,k)+w(i,j+1,k)+abs(w(i,j,k))+abs(w(i,j+1,k)))*dt/dz
  CwN(i,j,k) = 0.25*(w(i,j,k)+w(i,j+1,k)-abs(w(i,j,k))-abs(w(i,j+1,k)))*dt/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = v(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k)+u(i,j+1,k)-u(i,j,k-1)-u(i,j+1,k-1))/dx
  div2 = 0.5*(v(i,j+1,k)-v(i,j-1,k))/dy
  div3 = 0.5*(w(i,j,k)+w(i,j+1,k)-w(i+1,j,k)-w(i+1,j+1,k))/dz
  div = dt*B(i,j,k)*(div1+div2+div3)  
  advy(i,j,k)= BN(i,j,k)+div
END DO
END DO
END DO

! calculate the nonlinear terms for w-momentum equation
DO i = 1,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  CuP(i,j,k) = 0.25*(u(i,j,k)+u(i-1,j,k)+abs(u(i,j,k))+abs(u(i-1,j,k)))*dt/dx
  CuN(i,j,k) = 0.25*(u(i,j,k)+u(i-1,j,k)-abs(u(i,j,k))-abs(u(i-1,j,k)))*dt/dx
  CvP(i,j,k) = 0.25*(v(i,j,k)+v(i-1,j,k)+abs(v(i,j,k))+abs(v(i-1,j,k)))*dt/dy
  CvN(i,j,k) = 0.25*(v(i,j,k)+v(i-1,j,k)-abs(v(i,j,k))-abs(v(i-1,j,k)))*dt/dy
  CwP(i,j,k) = 0.25*(w(i,j,k)+w(i-1,j,k)+abs(w(i,j,k))+abs(w(i-1,j,k)))*dt/dz
  CwN(i,j,k) = 0.25*(w(i,j,k)+w(i-1,j,k)-abs(w(i,j,k))-abs(w(i-1,j,k)))*dt/dz
END DO
END DO
END DO

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  B(i,j,k) = w(i,j,k)
END DO
END DO
END DO

CALL advect

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k)+u(i-1,j,k)-u(i,j,k-1)-u(i-1,j,k-1))/dx
  div2 = 0.5*(v(i,j,k)+v(i-1,j,k)-v(i,j-1,k)-v(i-1,j-1,k))/dy
  div3 = 0.5*(w(i-1,j,k)-w(i+1,j,k))/dz
  div = dt*B(i,j,k)*(div1+div2+div3)  
  advz(i,j,k)= BN(i,j,k) + div
END DO
END DO
END DO

! calculate ustar, vstar and wstar 
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
IF(wet(i,j,k))THEN
! diffusion of u
 ahe = ah(i,j,k+1)
 dif1 = ahe*(u(i,j,k+1)-u(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 ahw = ah(i,j,k)
 dif2 = ahw*(u(i,j,k)-u(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 ahn = 0.25*(ah(i,j,k)+ah(i,j,k+1)+ah(i,j+1,k)+ah(i,j+1,k+1))
 dif3 = ahn*(u(i,j+1,k)-u(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 IF(land(j+1,k))  dif3 = ahn*(u(i,j+1,k)-slip*u(i,j,k))/dy
 ahs = 0.25*(ah(i,j,k)+ah(i,j,k+1)+ah(i,j-1,k)+ah(i,j-1,k+1))
 dif4 = ahs*(u(i,j,k)-u(i,j-1,k))/dy
 IF(dry(i,j-1,k)) dif4 = 0.0
 IF(land(j-1,k))  dif4 = ahs*(slip*u(i,j,k)-u(i,j-1,k))/dy 
 difh = (dif1-dif2)/dx + (dif3-dif4)/dy
 azt = 0.25*(az(i,j,k)+az(i-1,j,k)+az(i,j,k+1)+az(i-1,j,k+1))
 dif1 = azt*(u(i-1,j,k)-u(i,j,k))/dz
 if(i==isoli(j)+1) dif1=-tbu(j,k) 
 azb = 0.25*(az(i,j,k)+az(i+1,j,k)+az(i,j,k+1)+az(i+1,j,k+1))
 dif2 = azb*(u(i,j,k)-u(i+1,j,k))/dz
 IF(i == nz ) dif2 = 0.
 difz = (dif1-dif2)/dz
 diffu = dt*(difh+difz)  

! diffusion of v
 ahe = 0.25*(ah(i,j,k)+ah(i,j+1,k)+ah(i,j,k+1)+ah(i,j+1,k+1))
 dif1 = ahe*(v(i,j,k+1)-v(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 IF(land(j,k+1))  dif1 =  ahe*(v(i,j,k+1)-slip*v(i,j,k))/dx
 ahw = 0.25*(ah(i,j,k)+ah(i,j+1,k)+ah(i,j,k-1)+ah(i,j+1,k-1))
 dif2 = ahw*(v(i,j,k)-v(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 IF(land(j,k-1)) dif2 = ahw*(slip*v(i,j,k)-v(i,j,k-1))/dx
 ahn = ah(i,j+1,k)
 dif3 = ahn*(v(i,j+1,k)-v(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 ahs = ah(i,j,k)
 dif4 = ahs*(v(i,j,k)-v(i,j-1,k))/dx
 IF(dry(i,j-1,k)) dif4 = 0.0
 difh = (dif1-dif2)/dx + (dif3-dif4)/dy 
 azt = 0.25*(az(i,j,k)+az(i,j+1,k)+az(i-1,j,k)+az(i-1,j+1,k))
 dif1 = azt*(v(i-1,j,k)-v(i,j,k))/dz
 if(i==isoli(j)+1) dif1=-tbv(j,k) 
 azb = 0.25*(az(i,j,k)+az(i,j+1,k)+az(i+1,j,k)+az(i+1,j+1,k))
 dif2 = azb*(v(i,j,k)-v(i+1,j,k))/dz
 IF(i == nz ) dif2 = 0.
 difz = (dif1-dif2)/dz
 diffv = dt*(difh+difz)  

! diffusion of w
 ahe = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j,k+1)+ah(i-1,j,k+1))
 dif1 = ahe*(w(i,j,k+1)-w(i,j,k))/dx
 IF(dry(i,j,k+1)) dif1 = 0.0
 IF(land(j,k+1))  dif1 =  ahe*(w(i,j,k+1)-slip*w(i,j,k))/dx
 ahw = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j,k-1)+ah(i-1,j,k-1))
 dif2 = ahw*(w(i,j,k)-w(i,j,k-1))/dx
 IF(dry(i,j,k-1)) dif2 = 0.0
 IF(land(j,k-1)) dif2 = ahw*(slip*w(i,j,k)-w(i,j,k-1))/dx
 ahn = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j+1,k)+ah(i-1,j+1,k))
 dif3 = ahn*(w(i,j+1,k)-w(i,j,k))/dy
 IF(dry(i,j+1,k)) dif3 = 0.0
 IF(land(j+1,k))  dif3 = ahn*(w(i,j+1,k)-slip*w(i,j,k))/dy 
 ahs = 0.25*(ah(i,j,k)+ah(i-1,j,k)+ah(i,j-1,k)+ah(i-1,j-1,k))
 dif4 = ahs*(w(i,j,k)-w(i,j-1,k))/dy
 IF(dry(i,j-1,k)) dif4 = 0.0
 IF(land(j-1,k))  dif4 = ahs*(slip*w(i,j,k)-w(i,j-1,k))/dy 
 difh = (dif1-dif2)/dx + (dif3-dif4)/dy

 azt = 0.5*(az(i,j,k)+az(i-1,j,k))
 dif1 = azt*(w(i-1,j,k)-w(i,j,k))/dz
 azb = 0.5*(az(i,j,k)+az(i+1,j,k)) 
 dif2 = azb*(w(i,j,k)-w(i+1,j,k))/dz
 if(i==nz) dif2=0.
 difz = (dif1-dif2)/dz
 diffw = dt*(difh+difz)  

  IF(wet(i,j,k))THEN
    um = 0.25*(u(i,j,k)+u(i,j,k-1)+u(i,j+1,k)+u(i,j+1,k-1))
    vm = 0.25*(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1))
    pressx = -drdx*(q(i,j,k+1)-q(i,j,k))-drdx*(p(i,j,k+1)-p(i,j,k))   +rhoa/rho0*gv*sn
    pressx = pressx - gv*rho(i,j,k)/rho0*sn 
    fm = f
    alpha = fm*dt
    IF(wet(i,j,k+1)) ustar(i,j,k) = cos(alpha)*u(i,j,k) + sin(alpha)*vm + dt*pressx + advx(i,j,k) + diffu   -2.*om*cst*csb*w(i,j,k)   
    pressy = -drdy*(q(i,j+1,k)-q(i,j,k))-drdy*(p(i,j+1,k)-p(i,j,k))
    IF(wet(i,j+1,k)) vstar(i,j,k) = cos(alpha)*v(i,j,k) - sin(alpha)*um + dt*pressy + advy(i,j,k) + diffv   -2.*om*(cst*snb*cs-snt*sn)*w(i,j,k)   
    pressz = -drdz*(q(i-1,j,k)-q(i,j,k))
    IF(wet(i-1,j,k)) wstar(i,j,k) = w(i,j,k) + dt*pressz + advz(i,j,k) + diffw   -2.*om*((cst*snb*cs+snt*sn)*v(i,j,k)-cst*csb*u(i,j,k)) 
  END IF

END IF    
END DO
END DO
END DO

! boundary conditions
DO j = 1,ny
  isl=isoli(j)
DO i = isl+1,nz  
  if(j>=int(8.0d3/dy)+1 .and. j<=int(8.0d3/dy)+1+int(pw/dy) .and. i<=isl+isdp)then    
    if(i<=isl+isd)then
      if(i==isl+1)then 
        ustar(i,j,0)=uad
      else
        ustar(i,j,0)=ustar(i-1,j,0)   
      endif
      vstar(i,j,0)=vstar(i,j,1)
      wstar(i,j,0)=wstar(i,j,1)  
    else
      ustar(i,j,0)=ustar(i,j,1) 
      vstar(i,j,0)=vstar(i,j,1) 
      wstar(i,j,0)=wstar(i,j,1) 
    endif       
  else
    ustar(i,j,0)=0.
    vstar(i,j,0)=0.
    wstar(i,j,0)=0.
  endif
  ustar(i,j,nx+1)=ustar(i,j,nx)
  vstar(i,j,nx+1)=vstar(i,j,nx)
  wstar(i,j,nx+1)=wstar(i,j,nx)   
END DO
END DO

DO i = 1,nz
DO k = 1,nx
 ustar(i,0,k)=ustar(i,ny,k)
 ustar(i,ny+1,k)=ustar(i,1,k)
 vstar(i,0,k)=vstar(i,ny,k)
 vstar(i,ny+1,k)=vstar(i,1,k)
 wstar(i,0,k)=wstar(i,ny,k)
 wstar(i,ny+1,k)=wstar(i,1,k)
END DO
END DO

! calculate right-hand side of Poisson equation
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
 qstar(i,j,k) = -1.*rho0/dt*(  &
 &  (ustar(i,j,k)-ustar(i,j,k-1))*dz + &
 &  (vstar(i,j,k)-vstar(i,j-1,k))*dx*dz/dy + &
 &  (wstar(i,j,k)-wstar(i+1,j,k))*dx )
END DO
END DO
END DO

! STEP 4: S.O.R. ITERATION

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
 dq(i,j,k) = 0.0
END DO
END DO
END DO

nstop = 8000

!*****************
DO nsor = 1,8000
!*****************

perr = 0.0

! STEP 1: predict pressure correction
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx

 IF(wet(i,j,k))THEN
 q1 = dq(i,j,k)
 term1 = qstar(i,j,k)  + &  
  &      at(i,j,k)*dq(i-1,j,k) + ab(i,j,k)*dq(i+1,j,k) + & 
  &      aw(i,j,k)*dq(i,j,k-1) + ae(i,j,k)*dq(i,j,k+1) + &
  &      as(i,j,k)*dq(i,j-1,k) + an(i,j,k)*dq(i,j+1,k)
 q2 = (1.0-omega)*q1 + omega*term1/atot(i,j,k) 
 dq(i,j,k) = q2
 perr = MAX(ABS(q2-q1),perr)
 END IF

END DO
END DO
END DO

! STEP 2: predict new velocities 
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  IF(wet(i,j,k))THEN
    pressx = -drdx*(dq(i,j,k+1)-dq(i,j,k))
    IF(wet(i,j,k+1)) un(i,j,k) = ustar(i,j,k) + dt*pressx
    pressy = -drdy*(dq(i,j+1,k)-dq(i,j,k))
    IF(wet(i,j+1,k)) vn(i,j,k) = vstar(i,j,k) + dt*pressy
    pressz = -drdz*(dq(i-1,j,k)-dq(i,j,k))
    IF(wet(i-1,j,k)) wn(i,j,k) = wstar(i,j,k) + dt*pressz
  END IF
END DO
END DO
END DO

! STEP 3a: predict depth-integrated flow
DO j = 1,ny
DO k = 1,nx
  usum(j,k) = 0.
  vsum(j,k) = 0.
  do i = isoli(j)+1,nz  
    usum(j,k) = usum(j,k) + dz*un(i,j,k)
    vsum(j,k) = vsum(j,k) + dz*vn(i,j,k)  
  enddo
END DO
END DO

! lateral boundary conditions
DO j = 1,ny 
 usum(j,0)=0.
 usum(j,nx+1) = 0.
 vsum(j,0)=0.
 vsum(j,nx+1) = vsum(j,nx)
 do i = isoli(j)+1,nz 
 if(j>=int(8.0d3/dy)+1 .and. j<=int(8.0d3/dy)+1+int(pw/dy) .and. i<=isoli(j)+isdp)then    
   if(i<=isl+isd)then
     usum(j,0) = usum(j,0) + dz*uad
     vsum(j,0) = vsum(j,0) + dz*vn(i,j,1) 
   else
    usum(j,0) = usum(j,0) + dz*un(i,j,1)   
    vsum(j,0) = vsum(j,0) + dz*vn(i,j,1)
   endif
 else
 endif
 usum(j,nx+1) = usum(j,nx+1) + dz*un(i,j,nx)
 enddo 
END DO

DO k = 1,nx 
 vsum(0,k)=vsum(ny,k)
 vsum(ny+1,k)=vsum(1,k)
 usum(0,k)=usum(ny,k)
 usum(ny+1,k)=usum(1,k)  
END DO

! STEP 3b: predict surface pressure field
DO j = 1,ny
DO k = 1,nx
 dq(isoli(j),j,k) = -dt*rho0*gv*( (usum(j,k)-usum(j,k-1) )/dx + ( vsum(j,k)-vsum(j-1,k) )/dy ) *cs
END DO
END DO

IF(perr <= peps)THEN
  nstop = nsor
  GOTO 33
END IF

!********************
END DO
!********************

    GOTO 34
33  WRITE(*,*) "No. of Interactions =>", nstop, az(isoli(ny/2),ny/2,nx/2), ps(isoli(ny/2)+1,ny/2,nx/2), pb(isoli(ny/2)+1,ny/2,nx/2), taz(isoli(ny/2)+1,ny/2,nx/2), n22(isoli(ny/2)+1,ny/2,nx/2), ssbb(ny/2,nx/2), s(isoli(ny/2),ny/2,nx/2), s(isoli(ny/2)+1,ny/2,nx/2), ttbb(ny/2,nx/2),  t(isoli(ny/2),ny/2,nx/2), t(isoli(ny/2)+1,ny/2,nx/2), rho(isoli(ny/2),ny/2,nx/2), rho(isoli(ny/2)+1,ny/2,nx/2)
    if(nstop>2500) pause
34  CONTINUE

! updating for next time step
DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
if(wet(i,j,k))then
  q(i,j,k) = q(i,j,k)+dq(i,j,k)
  u(i,j,k) = un(i,j,k)
  v(i,j,k) = vn(i,j,k)
  w(i,j,k) = wn(i,j,k)
endif
END DO
END DO
END DO

DO j = 1,ny
DO k = 1,nx
  q(isoli(j),j,k) = q(isoli(j),j,k)+dq(isoli(j),j,k)
END DO
END DO

! lateral boundary conditions
DO j = 1,ny
 DO i = isoli(j),nz  
   q(i,j,0) = q(i,j,1)
   q(i,j,nx+1) = q(i,j,nx)
 END DO
END DO

DO k = 1,nx
 DO i = isoli(1),nz
   q(i,0,k) = q(i,ny,k)
 END DO
 DO i = isoli(ny),nz
   q(i,ny+1,k) = q(i,1,k)
 END DO 
END DO
  
RETURN
END SUBROUTINE dyn

SUBROUTINE advect
implicit none
! local parameters
real(8) :: RxP(0:nz+1,0:ny+1,0:nx+1), RxN(0:nz+1,0:ny+1,0:nx+1)
real(8) :: RyP(0:nz+1,0:ny+1,0:nx+1), RyN(0:nz+1,0:ny+1,0:nx+1)
real(8) :: RzP(0:nz+1,0:ny+1,0:nx+1), RzN(0:nz+1,0:ny+1,0:nx+1)
real(8) :: dB, term1, term2, term3, term4, term5, term6
real(8) :: BwP, BwN, BeP, BeN
real(8) :: BnP, BnN, BsP, BsN
real(8) :: BbP, BbN, BtP, BtN 

DO i = 0,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  RxP(i,j,k) = 0.0
  RxN(i,j,k) = 0.0
  RyP(i,j,k) = 0.0
  RyN(i,j,k) = 0.0
  RzP(i,j,k) = 0.0
  RzN(i,j,k) = 0.0
END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  dB =  B(i,j,k+1)-B(i,j,k)
  IF(ABS(dB) > 0.0) RxP(i,j,k) = (B(i,j,k)-B(i,j,k-1))/dB
  dB =  B(i,j+1,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RyP(i,j,k) = (B(i,j,k)-B(i,j-1,k))/dB
  dB =  B(i-1,j,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RzP(i,j,k) = (B(i,j,k)-B(i+1,j,k))/dB
END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 0,nx-1
  dB =  B(i,j,k+1)-B(i,j,k)
  IF(ABS(dB) > 0.0) RxN(i,j,k) = (B(i,j,k+2)-B(i,j,k+1))/dB
END DO
END DO
END DO

DO i = 1,nz
DO j = 0,ny-1
DO k = 1,nx
  dB =  B(i,j+1,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RyN(i,j,k) = (B(i,j+2,k)-B(i,j+1,k))/dB
END DO
END DO
END DO

DO i = 2,nz+1
DO j = 1,ny
DO k = 1,nx
  dB =  B(i-1,j,k)-B(i,j,k)
  IF(ABS(dB) > 0.0) RzN(i,j,k) = (B(i-2,j,k)-B(i-1,j,k))/dB
END DO
END DO
END DO   

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx

term1 = (1.0-CuP(i,j,k-1))*(B(i,j,k)-B(i,j,k-1))
BwP = B(i,j,k-1)+0.5*PSI(RxP(i,j,k-1))*term1

term1 = (1.0+CuN(i,j,k-1))*(B(i,j,k)-B(i,j,k-1))
BwN = B(i,j,k)-0.5*PSI(RxN(i,j,k-1))*term1

term1 = (1.0-CuP(i,j,k))*(B(i,j,k+1)-B(i,j,k))
BeP = B(i,j,k)+0.5*PSI(RxP(i,j,k))*term1

term1 = (1.0+CuN(i,j,k))*(B(i,j,k+1)-B(i,j,k))  
BeN = B(i,j,k+1)-0.5*PSI(RxN(i,j,k))*term1

term1 = (1.0-CvP(i,j-1,k))*(B(i,j,k)-B(i,j-1,k))
BsP = B(i,j-1,k)+0.5*PSI(RyP(i,j-1,k))*term1

term1 = (1.0+CvN(i,j-1,k))*(B(i,j,k)-B(i,j-1,k))
BsN = B(i,j,k)-0.5*PSI(RyN(i,j-1,k))*term1

term1 = (1.0-CvP(i,j,k))*(B(i,j+1,k)-B(i,j,k))
BnP = B(i,j,k)+0.5*PSI(RyP(i,j,k))*term1

term1 = (1.0+CvN(i,j,k))*(B(i,j+1,k)-B(i,j,k))  
BnN = B(i,j+1,k)-0.5*PSI(RyN(i,j,k))*term1

term1 = (1.0-CwP(i+1,j,k))*(B(i,j,k)-B(i+1,j,k))
BbP = B(i+1,j,k)+0.5*PSI(RzP(i+1,j,k))*term1

term1 = (1.0+CwN(i+1,j,k))*(B(i,j,k)-B(i+1,j,k))  
BbN = B(i,j,k)-0.5*PSI(RzN(i+1,j,k))*term1

term1 = (1.0-CwP(i,j,k))*(B(i-1,j,k)-B(i,j,k)) 
BtP = B(i,j,k)+0.5*PSI(RzP(i,j,k))*term1

term1 = (1.0+CwN(i,j,k))*(B(i-1,j,k)-B(i,j,k)) 
BtN = B(i-1,j,k)-0.5*PSI(RzN(i,j,k))*term1


term1 = CuP(i,j,k-1)*BwP+CuN(i,j,k-1)*BwN
term2 = CuP(i,j,k)*BeP+CuN(i,j,k)*BeN

term3 = CvP(i,j-1,k)*BsP+CvN(i,j-1,k)*BsN
term4 = CvP(i,j,k)*BnP+CvN(i,j,k)*BnN

term5 = CwP(i+1,j,k)*BbP+CwN(i+1,j,k)*BbN
term6 = CwP(i,j,k)*BtP+CwN(i,j,k)*BtN

BN(i,j,k) = term1-term2+term3-term4+term5-term6

END DO
END DO
END DO

RETURN

END SUBROUTINE advect

REAL(8) FUNCTION psi(r)

! input parameters
implicit none
real(8), INTENT(IN) :: r  

! local parameters 
real(8) :: term1, term2, term3

  term1 = MIN(2.0*r,1.0)
  term2 = MIN(r,2.0)
  term3 = MAX(term1,term2)
  psi = MAX(term3,0.0)

RETURN

END FUNCTION psi 

END MODULE sub