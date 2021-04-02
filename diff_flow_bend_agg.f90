!-----------------------------------------------------------------------------------------------------------------
! Curvature driven feedback on aggregation-diffusion of proteins in lipid bilayers
!
! Simulation code by: Arijit Mahapatra, David Saintillan, Padmini Ranganami (UCSD)
!
! For details of the model, refer to:
! "Curvature driven feedback on aggregation-diffusion of proteins in lipid bilayers" by Mahapatra et al.
! (arxiv_link)
!-----------------------------------------------------------------------------------------------------------------
program sumingup
implicit none
integer:: N,M,P,i,j,k,ptc,itn,lamit,pvti,pvtj,velit,ii,jj,ilg,jlg,lgp,Pw,tnn,iisig
double precision:: sig0, lam0, bsig_s, ell, kp, beta,c,err,Tst,Ted,Tstp,stth,flg,dt1
double precision:: pr, L, D1, D2, D3, D4, D5, D6, D7, D8, D9,start,finish,dlr,fc,kBT,D
double precision:: D10,lr,tT,lX,lY,dx,dy,dt,Dmf,xa,xb,ya,yb,dialn,lgsigdx,lgsigdy,dsig
double precision:: cth,lct,tol,sigtol,ur,zr,zrr,sigrr,err_zeta,err_lam,err_z,totsig
double precision:: zer,zerr,lamrr,lamr,rZe,zenu,znu,vrr,diln2,test,D13,D14,D15,D16,D17,D18
double precision:: signu,dzetae,dzetaw,dzetan,dzetas,zcon,zecon,velerr,utol,mf,mff,zzcon
double precision:: zetadx,zetady,azE,azW,azN,azS,azP,uer,ver,xd,yd,etadx,etady,isig,tsig
double precision:: sigdx,sigdy,Lsig,Lsigsq,sigo,uW,uE,uN,uS,urr,uP,sigsqdy,D11,dts,ifl
double precision:: zeo,zo,azeE,azeW,azeN,azeS,sigsqdx,dvy,dvx,duy,dux,d2zy,d2zxy,d2zx
double precision:: dsign,dsige,dsigw,dsigs,umean,vmean,usum,vsum,xt,yt,r,PI,ddx,ddy,wbc 
double precision:: laW,laS,laP,laN,lamo,laE,lamnu,zera,uup,vup,nu,uo,vo,delta,divm,xc,yc
double precision:: zx,zy,wold,wdx,wdy,hxx,hxy,hyy,wzedx,wzedy,Lzeta,Lzew,wxx,wyy,wxy
double precision::Sumean,Svmean,g1mean,g2mean,Slamean,err_sig,D12,ztotal,sigmax
double precision:: gamma,hatB,hatL,hatA,hatS,hatT,Pe
double precision, dimension(0:500):: x,y,siggs,sigge,siggw,siggn,zh,zv,zetah,wh,wv
double precision, dimension(0:500)::zetav,uh,uv,vh,vv,lamh,lamv,sigh,sigv,etah,etav
double precision, dimension(0:5000):: t,Je,Jw,Jn,Js
double precision, dimension(0:20000000)::xLg,yLg
double precision, dimension(0:500,0:500):: sigr,azeP,Sze,d2h,Sz,Ssi,Sla,shst,uold,vold
double precision, dimension(0:500,0:500):: zeta,z,sig,lam,sigp,u,v,Sv,Su,zold,zetaold
double precision, dimension(0:500,0:500)::sigold
double precision, dimension(0:500,0:500)::asiE,asiW,asiN,asiS,g1,g2,fx,fy,div,asiP,Lpzeta
double precision, dimension(0:500,0:500):: oldu,oldv,oldlam,oldsig,w,zp,g0,eta,Leta 
integer nx,ny
double precision ffx(64,64),ffy(64,64),ux(64,64),uy(64,64), rdm(64,64)
double precision pp(64,64),laf(64,64),lamS(64,64)
double precision S(64,64),zetf(64,64),zf(64,64),zetaS(64,64)
!double precision Lx,Ly

!!Specify Grid numbers, change the size of the variable defined in line 36-39 along with 
!!changing grid number 

nx=64
ny=64
!!! size of the domain      
Lx=1.d0
Ly=1.d0

!! generating a 2D array of random number 
CALL RANDOM_NUMBER(rdm)

!! define parameters here or go to line 66 and define all the primary  Dimensionless number
sig0=2d-4 !! sig_s (#/nm^2)
lam0=1d-4  !! (pN/nm)
ell=8d0  !! spontaneous curvature (nm)
kp=84.d0  !! Bending rigidity
kBT=4.14d-2  !! Thermal energy
!gamma=25*kBT  !! aggregation coefficient 
D=1e5     !! diffusion coefficient (nm^2/s)
nu=5d-6   !! viscosity (pN-S/nm)
pr=0.d0   !! pressure  (pN/nm^2)
L=1000.d0  !! size (nm)

!! Primary dimensionless number
!! hatB= kBT/kp
!! hatL= ell/L
!! hatA= gamma/kBT
!! hatS= sig0*L**2
!! hatT= 2.d0*L**2.d0*lam0/kp
!! Pe= lam0*L**2.d0/(nu*D)





!! The internal dimensionless number calculated based on the parameters 
PI=2.D0*DASIN(1.D0)
D1= kp*ell*sig0/(L*lam0)  !! 2*hatL*hatS/hatT
D2=2.d0*(ell)**2.d0*kp*sig0**2.d0/lam0  !! 4*hatL**2*hatS**2/hatT
D3=kBT*sig0/lam0  !!  2*hatB*hatS/hatT
D5=2.d0*L*ell*sig0  !! 2*hatL*hatS
D6=L**2.d0*kBT*sig0/kp !!  hatS*hatB
D7=2.d0*kp*ell**2.d0*sig0/kBT !! 2*hatL**2*hatS/hatB
D8=2.d0*L**2.d0*lam0/kp !! hatT
D9=2.d0*L**3.d0*pr/kp !! ND for pressure (zero)
D10=kp*ell/(kBT*L)  !! hatL/hatB
D11=(lam0*L**2.d0)/(nu*D) !! Pe

D17=25.d0  !! hatA
D16=D17*D3/2.d0
D13=D16/(sig0*L**2.d0)
!! derived ones
D14=D13/D3  !! hatA/(2*hatS)
D15=D13*D8/4.d0

!D17=2.d0*D16/D3
D18=D16*D8/2.d0

ifl=0

N=nx
M=ny
tnn=1001
iisig=50000
P=tnn*iisig
lr=1.d0
tT=3d-1
lX=1.d0
lY=lX/lr
dx=lX/(N)
dy=lY/(M)
delta=(dx**2.d0+dy**2.d0)**0.5d0
!delta=dx
y(1)=-lY/2.d0
do j=2,M
y(j)=y(j-1)+dy
enddo
x(1)=-lX/2.d0
do i=2,N
x(i)=x(i-1)+dx
enddo
dt=tT/(tnn-1)
dts=dt/iisig
t(1)=0.d0
do k=2,tnn
t(k)=t(k-1)+dt
enddo


do i=0,N+1
do j=0,M+1
zeta(i,j)=0.d0
sig(i,j)=0.d0
z(i,j)=0.d0
lam(i,j)=1.d0
d2h(i,j)=0.d0
shst(i,j)=0
u(i,j)=0.d0
v(i,j)=0.d0
fx(i,j)=0.d0
fy(i,j)=0.d0
w(i,j)=0.d0
eta(i,j)=0.d0
enddo
enddo
!Dmf=kp*muph/(2.d0*(alpha**2.d0+kp*(muph)**2.d0)*sig0*L)

do j=1,M
z(1,j)=0.d0
z(N,j)=0.d0
lam(N,j)=1.d0
lam(1,j)=1.d0
enddo
do i=1,N
z(i,1)=0.d0
z(i,M)=0.d0
lam(i,M)=1.d0
lam(i,1)=1.d0
enddo

Tst=tT/10.d0;
Ted=tT/5.d0;
Tstp=tT/8.d0;
stth=10.d0;



do k=2,M
    !Je(k)=-0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))
    !Jw(k)=0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))  
    !Jn(k)=-0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))
    !Js(k)=0.5d0*(tanh(stth*(t(k)-Tst))-tanh(5*stth*(t(k)-Tstp)))  
    !Je(k)=0.d0
    !Jw(k)=0.d0 
    !Jn(k)=0.d0
    !Js(k)=0.d0  
enddo


!! calculation of log
ilg=1

          do i=1,N   ! open 53
          if (i==1 .OR. i==N) then !open 53a
          do j=1,M    ! open 54
        
          do ii=2,N-1    !open 55
          do jj=2,M-1    ! open 56
          !if (ii/=i .OR. jj/=j)then
          xt=(x(i)-x(ii))
          yt=(y(j)-y(jj))
          r=(xt**2.d0+yt**2.d0)**(0.5d0)
          xlg(ilg)=DLOG(r)
          ilg=ilg+1
          !endif
          enddo     ! close 56
          enddo    ! close 55
          enddo    ! close 54
          endif ! close 53a
          enddo    ! close 53
 jlg=1
 
        do j=1,M   ! open 53
          if (j==1 .OR. j==M) then !open 53a
          do i=2,N-1     ! open 54
          do ii=2,N-1    !open 55
          do jj=2,M-1    ! open 56
          !if (ii/=i .OR. jj/=j)then
          xt=(x(i)-x(ii))
          yt=(y(j)-y(jj))
          r=(xt**2.d0+yt**2.d0)**(0.5d0)
          ylg(jlg)=DLOG(r)
          jlg=jlg+1
          !endif
          enddo     ! close 56
          enddo    ! close 55
          enddo    ! close 54
          endif ! close 53a
          enddo    ! close 53       


cth=20.d0
lct=0.89d-1

!xa=0.d0
!ya=0.d0

xa=-0.25d0
ya=-0.25d0

xb=0.25d0
yb=0.25d0

xc=-0.25d0
yc=0.25d0

xd=0.25d0
yd=-0.25d0

!! Giving initial value (perturbation)
do i=1,N
do j=1,M
!sig(i,j)=0.1d0+0.25d0*(1.d0-tanh(cth*(((x(i)-xa)**2.d0+(y(j)-ya)**2.d0)**0.5d0-lct))) + &
!0.25d0*(1.d0-tanh(cth*(((x(i)-xb)**2.d0+(y(j)-yb)**2.d0)**0.5d0-lct)))  + &
!0.25d0*(1.d0-tanh(cth*(((x(i)-xc)**2.d0+(y(j)-yc)**2.d0)**0.5d0-lct))) !+ &
!0.25d0*(1.d0-tanh(cth*(((x(i)-xd)**2.d0+(y(j)-yd)**2.d0)**0.5d0-lct)))
!sig(i,j)=1.d0-tanh(9.d0*((x(i)**2.d0+y(j)**2.d0)**0.5d0-2d-1))
!sig(i,j)=(tanh(cth*(x(i)+lct))-tanh(cth*(x(i)-lct)))*(tanh(cth*(y(j)+lct))- &
!tanh(cth*(y(j)-lct)))
sig(i,j)=0.1d0 +0.0001d0*rdm(i,j)
enddo
enddo

totsig=0.d0
          do i=1,N
          do j=1,M
          totsig=totsig+dx*dy*sig(i,j)
          enddo
          enddo

tol=5d-7
sigtol=1d-7
!utol=1d-10
ur=7d-1
!! coefficients
 		 azeE=1.d0/dx**2.d0
         azeW=1.d0/dx**2.d0
         azeN=1.d0/dy**2.d0
         azeS=1.d0/dy**2.d0
         
          azE=1.d0/dx**2.d0
          azW=1.d0/dx**2.d0
          azN=1.d0/dy**2.d0
          azS=1.d0/dy**2.d0
          azP=-(2.d0/dx**2.d0+2.d0/dy**2.d0)

        
        laE=1.d0/dx**2.d0
          laW=1.d0/dx**2.d0
          laN=1.d0/dy**2.d0
          laS=1.d0/dy**2.d0
          laP=-2.d0*(1.d0/dx**2.d0+1.d0/dy**2.d0)
          
          uE=1.d0/(dx**2.d0)
   		 uW=1.d0/(dx**2.d0)
   		 uN=1.d0/(dy**2.d0)
   		 uS=1.d0/(dy**2.d0)
   		 uP=-2.d0*(1.d0/(dx**2.d0)+1.d0/(dy**2.d0))




  	fc=0.d0

do k=1,P
if(k==1)then 
   isig=ifl
   Pw=1
  else
  isig=MOD((k-1),iisig)
  if(isig==0.d0)then
  Pw=Pw+1
  endif
  endif

if (k>1) then   ! open 007          ! just to write file

  	if(k>1)then
  	do i=1,N
  	do j=1,M
	sigp(i,j)=sig(i,j)
	zp(i,j)=z(i,j)
	fc=1.d0
	 enddo
	 enddo
	endif
	itn=1
	err=5.d0
	!if (k<3)then
	!dt=0.05d0
	!else
	!dt=dt1
	!endif
	
  
	
	
	!asiP=-(2.d0/dx**2.d0+2.d0/dy**2.d0+1.d0/dt)
        do while (err>tol .OR. itn<5)
        err=0.d0
         !do itn=1,20
        itn=itn+1
  !=====================zeta and z iteration=========================================
  !============================================================================
       
    if(isig==ifl)then   !! begining of time step hopping ..../\/\/\/\/\.../\/\/\/\/\../\/\/\ 
       
         zerr=0.d0
        
          
          
         do i=1,N
         do j=1,M
         zetaold(i,j)=zeta(i,j)
         enddo
         enddo
         
          do i=1,N
         do j=1,M
         zold(i,j)=z(i,j)
         enddo
         enddo
         
         
         
         
         
          do i=1,N
         do j=1,M
         azeP(i,j)=-(2.d0/dx**2.d0+2.d0/dy**2.d0)
         rZe=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/dx**2.d0+   &
        (sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/dy**2.d0
         Sze(i,j)=D5*rZe+D8*shst(i,j)+D6*zeta(i,j)*(2.d0*fc*(sig(i,j)*log(sig(i,j))+ &
         (1.d0-sig(i,j))*log(1.d0-sig(i,j)))+  & 
    D7*(sig(i,j))**2.d0-fc*(2.d0*D18/D6)*sig(i,j)*(sig(i,j)-1.d0)+  & 
    fc*(2.d0*D15/D6)*(((sig(i+1,j)-sig(i-1,j))/(2.d0*dx))**2.d0+  &
    ((sig(i,j+1)-sig(i,j-1))/(2.d0*dy))**2.d0)) +  &
    zeta(i,j)*D8*lam(i,j)
   
         enddo
         enddo
          
          
          do i=1,nx
          do j=1,ny
          zetaS(i,j)=Sze(i,j)
          enddo
          enddo
          
          
     call    poissol(Lx,Ly,nx,ny,zetaS,zetf)
     call    poissol(Lx,Ly,nx,ny,zetf,zf)
     
     
          do i=1,nx
          do j=1,ny
          zeta(i,j)=zetf(i,j)
          z(i,j)=zf(i,j)!-zf(1,1)
          enddo
          enddo
          
        
         do j=1,M
          zeta(0,j)=zeta(N,j)
          zeta(N+1,j)=zeta(1,j)
          z(0,j)=z(N,j)
          z(N+1,j)=z(1,j)
          enddo
          
          do i=1,N
          zeta(i,0)=zeta(i,M)
          zeta(i,M+1)=zeta(i,1)
           z(i,0)=z(i,M)
          z(i,M+1)=z(i,1)
          enddo
         

        
         
              ! now check convergence
     err=0.d0
     err_zeta=0.d0
        do i=1,N
        do j=1,M
        if(abs(zeta(i,j)-zetaold(i,j))>err_zeta) then
        err_zeta=abs(zeta(i,j)-zetaold(i,j))
        endif
        
        enddo
        enddo
        err=err_zeta
        err_z=0.d0
          do i=2,N-1
        do j=2,M-1
        if(abs(z(i,j)-zold(i,j))>err_z) then
        err_z=abs(z(i,j)-zold(i,j))
        endif
        
        enddo
        enddo
        
        
        
        if (err_z>err)then
        err=err_z
        endif
        
        !print*, "time step",k, "z iteration error is",zcon,"err",err,"itn",itn,"err_lam",err_lam
        
       if (k>1)then
        do i=1,N
        do j=1,M
        wold=w(i,j)
        w(i,j)=(ur)*w(i,j)+(1.d0-ur)*((z(i,j)-zp(i,j))/dt)/D11
        err=max(abs(wold-w(i,j)),err)
        enddo
        enddo
       endif
       
       
       do i=1,N
  		 do j=1,M
  		 oldlam(i,j)=lam(i,j)
  		 enddo
  		 enddo 
       
  !!!!! Periodic_stokes_solver starts here=================================
  !! ===== This section uses subroutine given below
   
   if(k>1)then !! periodic_enter
   
   do i=1,N
  do j=1,M
  oldu(i,j)=u(i,j)
  oldv(i,j)=v(i,j)
  enddo
  enddo

  
  
          do i=1,N
          do j=1,M 
          zetadx=(zeta(i+1,j)-zeta(i-1,j))/(2.d0*dx)
          zetady=(zeta(i,j+1)-zeta(i,j-1))/(2.d0*dy)
          wdx=(w(i+1,j)-w(i-1,j))/(2.d0*dx)
          wdy=(w(i,j+1)-w(i,j-1))/(2.d0*dy)
          wzedx=(zeta(i+1,j)*w(i+1,j)-zeta(i-1,j)*w(i-1,j))/(2.d0*dx)
          wzedy=(zeta(i,j+1)*w(i,j+1)-zeta(i,j-1)*w(i,j-1))/(2.d0*dy)
          hxx=(z(i+1,j)-2.d0*z(i,j)+z(i-1,j))/(dx**2.d0)
          hyy=(z(i,j+1)-2.d0*z(i,j)+z(i,j-1))/(dy**2.d0)
          hxy=(z(i+1,j+1)-z(i+1,j-1)-z(i-1,j+1)+z(i-1,j-1))/(4.d0*dx*dy)
          sigdx=(sig(i+1,j)-sig(i-1,j))/(2.d0*dx)
   		  sigdy=(sig(i,j+1)-sig(i,j-1))/(2.d0*dy)
          sigsqdx=((sig(i+1,j))**2.d0-(sig(i-1,j))**2.d0)/(2.d0*dx)
          sigsqdy=((sig(i,j+1))**2.d0-(sig(i,j-1))**2.d0)/(2.d0*dy)
fx(i,j)= -(D1*zeta(i,j)-D3*DLOG(abs(sig(i,j)/(1.d0-sig(i,j))))+D16*(2.d0*sig(i,j)-1.d0)+ &
          fc*D13*eta(i,j))*sigdx+5d-1*(D2)*sigsqdx- & 
     2.d0*w(i,j)*zetadx- &
  2.d0*(wdx*hxx+wdy*hxy) +wzedx
fy(i,j)= -(D1*zeta(i,j)-D3*DLOG(abs(sig(i,j)/(1.d0-sig(i,j))))+D16*(2.d0*sig(i,j)-1.d0)+ &
          fc*D13*eta(i,j))*sigdy+5d-1*(D2)*sigsqdy -  & 
   2.d0*w(i,j)*zetady- &
  2.d0*(wdx*hxy+wdy*hyy) +wzedy
          
          enddo
          enddo    
          
          g1mean=0.d0
          g2mean=0.d0
          do i=1,N !open 51
          do j=1,M ! open 52
          
        
          g1(i,j)=fx(i,j)*dx*dy
          g1mean=g1mean+g1(i,j)
          g2(i,j)=fy(i,j)*dx*dy
          g2mean=g2mean+g2(i,j)
          !Sumean=Sumean+Su(i,j)*ddx*ddy
          !Svmean=Svmean+Sv(i,j)*ddx*ddy
          enddo    ! close 52
          enddo    ! close 51 
       
       
       do i=1,N !open 51
          do j=1,M ! open 52
          fx(i,j)=fx(i,j)-g1mean
          fy(i,j)=fy(i,j)-g2mean
          enddo
          enddo
       
       do i=1,nx
       do j=1,ny
       ux(i,j)=u(i,j)
       uy(i,j)=v(i,j)
       ffx(i,j)=fx(i,j)
       ffy(i,j)=fy(i,j)
       S(i,j)=w(i,j)*zeta(i,j)
       enddo
       enddo
       
       call velfield(Lx,Ly,nx,ny,ux,uy,ffx,ffy,S,pp)
       
       
       
       do i=1,nx
       do j=1,ny
       u(i,j)=ux(i,j)
       v(i,j)=uy(i,j)
       lam(i,j)=1.d0-pp(i,j)
       enddo
       enddo
       
       do i=1,nx
       u(i,M+1)=u(i,1)
       v(i,M+1)=v(i,1)
       lam(i,M+1)=lam(i,1)
       u(i,0)=u(i,M)
       v(i,0)=v(i,M)
       lam(i,0)=lam(i,M)
       enddo
       
       do j=1,ny
       u(N+1,j)=u(1,j)
       v(N+1,j)=v(1,j)
       lam(N+1,j)=lam(1,j)
       u(0,j)=u(N,j)
       v(0,j)=v(N,j)
       lam(0,j)=lam(N,j)
       
       enddo
       
       
       do i=1,N
    do j=1,M
    if(abs(oldu(i,j)-u(i,j))>err)then
    !err=abs(oldu(i,j)-u(i,j))
    endif
    enddo
    enddo
    
    do i=1,N
    do j=1,M
    if(abs(oldv(i,j)-v(i,j))>err)then
    !err=abs(oldv(i,j)-v(i,j))
    endif
    enddo
    enddo
    
       
       z(0,0)=z(N,M)
z(N+1,M+1)=z(1,1)
z(0,M+1)=z(N,1)
z(N+1,0)=z(1,M)

do i=1,N
do j=1,M
d2zx=(z(i+1,j)-2.d0*z(i,j)+z(i-1,j))/(dx**2.d0)
d2zy=(z(i,j+1)-2.d0*z(i,j)+z(i,j-1))/(dy**2.d0)
d2zxy=(z(i+1,j+1)-z(i+1,j-1)-z(i-1,j+1)+z(i-1,j-1))/(4.d0*dx*dy)
dux=(u(i+1,j)-u(i-1,j))/(2.d0*dx)
dvy=(v(i,j+1)-v(i,j-1))/(2.d0*dy)
duy=(u(i,j+1)-u(i,j-1))/(2.d0*dy)
dvx=(v(i+1,j)-v(i-1,j))/(2.d0*dx)
shst(i,j)=2.d0*(d2zx*dux+d2zy*dvy+d2zxy*(duy+dvx))
enddo
enddo
        
    
    endif !! periodic_exit
    
 ! lam iteration starts========================================================
 !=============================================================================
 !!===== This loop is usefull to calculate equilibrium value of lambda===========
 !!==when a prescribed IC is given, like finite patches==========================
      if(k==1)then !!! Lam iteration starts, to enter or not to enter (Note this for this simulation this loop is not used, as k>1)
          do i=1,N
  		 do j=1,M
  		 oldlam(i,j)=lam(i,j)
  		 enddo
  		 enddo

         
          do i=1,N
          !wh(i)=w(i,0)
          w(i,0)=w(i,M)
          w(i,M+1)=w(i,1)

          enddo
          
          do j=1,M
          !wv(j)=w(0,j)
           w(0,j)=w(N,j)
          w(N+1,j)=w(1,j)

          enddo
          
          
          do i=1,N
          do j=1,M
          eta(i,j)=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/dx**2.d0+  &
         (sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/dy**2.d0
         !Leta(i,j)=(eta(i+1,j)-2.d0*eta(i,j)+eta(i-1,j))/dx**2.d0+  &
         !(eta(i,j+1)-2.d0*eta(i,j)+eta(i,j-1))/dy**2.d0
          enddo
          enddo
          
          
          do i=1,N
          do j=1,M
 Lzew=(w(i+1,j)*zeta(i+1,j)-2.d0*w(i,j)*zeta(i,j)+w(i-1,j)*zeta(i-1,j))/(dx**2.d0)+  &
 (w(i,j+1)*zeta(i,j+1)-2.d0*w(i,j)*zeta(i,j)+w(i,j-1)*zeta(i,j-1))/(dy**2.d0)
 Lzeta=(zeta(i+1,j)-2.d0*zeta(i,j)+zeta(i-1,j))/(dx**2.d0)+  &
         (zeta(i,j+1)-2.d0*zeta(i,j)+zeta(i,j-1))/(dy**2.d0)
         hxx=(z(i+1,j)-2.d0*z(i,j)+z(i-1,j))/(dx**2.d0)
          hyy=(z(i,j+1)-2.d0*z(i,j)+z(i,j-1))/(dy**2.d0)
          hxy=(z(i+1,j+1)-z(i+1,j-1)-z(i-1,j+1)+z(i-1,j-1))/(4.d0*dx*dy)
          wxx=(w(i+1,j)-2.d0*w(i,j)+w(i-1,j))/(dx**2.d0)
          wyy=(w(i,j+1)-2.d0*w(i,j)+w(i,j-1))/(dy**2.d0)
          wxy=(w(i+1,j+1)-w(i+1,j-1)-w(i-1,j+1)+w(i-1,j-1))/(4.d0*dx*dy)
          wbc=hxx*wxx+2.d0*hxy*wxy+hyy*wyy
          zetadx=(zeta(i+1,j)-zeta(i-1,j))/(2.d0*dx)
          zetady=(zeta(i,j+1)-zeta(i,j-1))/(2.d0*dy)
          wdx=(w(i+1,j)-w(i-1,j))/(2.d0*dx)
          wdy=(w(i,j+1)-w(i,j-1))/(2.d0*dy)
          sigdx=(sig(i+1,j)-sig(i-1,j))/(2.d0*dx)
          sigdy=(sig(i,j+1)-sig(i,j-1))/(2.d0*dy)
          
          etadx=(eta(i+1,j)-eta(i-1,j))/(2.d0*dx)
          etady=(eta(i,j+1)-eta(i,j-1))/(2.d0*dy)
          
          Lsig=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/dx**2.d0+  &
         (sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/dy**2.d0
         !eta(i,j)=Lsig
         
        Lsigsq=((sig(i+1,j))**2.d0-2.d0*(sig(i,j))**2.d0+(sig(i-1,j))**2.d0)/dx**2.d0+ &
         ((sig(i,j+1))**2.d0-2.d0*(sig(i,j))**2.d0+(sig(i,j-1))**2.d0)/dy**2.d0
          Sla(i,j)= (D1*zeta(i,j)-fc*D3*log(sig(i,j)/(1.d0-sig(i,j)))- &
          fc*D16+fc*D13*eta(i,j))*Lsig +   & 
D1*(zetadx*sigdx+zetady*sigdy)- fc*D3*(1.d0/(1.d0-sig(i,j)))*(sigdx*sigdx+sigdy*sigdy) + &
          fc*D13*(sigdx*etadx+sigdy*etady) -  &
    (5d-1)*(D2-2.d0*fc*D16)*Lsigsq-2.d0*Lzew+4.d0*(zetadx*wdx+zetady*wdy)+ &
    2.d0*w(i,j)*Lzeta+ &
        2.d0*wbc
          end do
          end do
          
        do i=1,nx
        do j=1,ny
        lamS(i,j)=Sla(i,j)
        enddo
        enddo  
          
          
          
       call   poissol(Lx,Ly,nx,ny,lamS,laf)
       
       
           do i=1,nx
        do j=1,ny
        lam(i,j)=1.d0+laf(i,j)
        enddo
        enddo  
          
          
          do i=1,N
        !lamh(i)=lam(i,0)
        lam(i,0)=lam(i,M)
        lam(i,M+1)=lam(i,1)
        
          enddo

    do j=1,M
    !lamv(j)=lam(0,j)
    lam(0,j)=lam(N,j)
    lam(N+1,j)=lam(1,j)
     enddo          
    
         
          
        endif !!! lam iteration to enter or not enter

    
    
     err_lam=0.d0
        
         do i=1,N
    	 do j=1,M
    	 if(abs(oldlam(i,j)-lam(i,j))>err_lam)then
    	 err_lam=abs(oldlam(i,j)-lam(i,j))
    	 endif
    	 enddo
    	 enddo  
          
          if(err_lam>err)then
          err=err_lam
          endif
    
    
    
    
    
    
          
       
!-======== div cal
divm=0.d0
   do i=1,N
   do j=1,M
   div(i,j)=(u(i+1,j)-u(i-1,j))/(2.d0*dx)+(v(i,j+1)-v(i,j-1))/(2.d0*dy)-zeta(i,j)*w(i,j)
   if(abs(div(i,j))>=divm)then
   divm=abs(div(i,j))
   endif
   enddo
   enddo
 !-======== div cal
    
diln2=(v(1,M)+v(N,M)-v(1,1)-v(N,1))*dx/2.d0+ (u(N,1)+u(N,M)-u(1,1)-u(1,M))*dy/2.d0
do i=2,N-1
diln2=diln2+v(i,M)*dx-v(i,1)*dx
enddo     
do j=2,M-1
diln2=diln2+u(N,j)*dy-u(1,j)*dy
enddo

         
  
  
  endif !!! endif for time hopping/\/\/\/\/\/\/\/\/\/\...../\/\/\/\/\/\/...../\/\/\/\/\/\
  
     
!!!!===================Aggregation-Diffusion_equation starts here============== 
!sig iteration starts==========================================================
!==============================================================================
if (k>1)then
 
 D12=D11

         do i=1,N
  		 do j=1,M
  		 oldsig(i,j)=sig(i,j)
  		 enddo
  		 enddo
          
         
          
           err_sig=10.d0
          do while (err_sig>sigtol)
          
          err_sig=0.d0
          
           
           do i=1,N
  		 do j=1,M
  		 sigold(i,j)=sig(i,j)
  		 enddo
  		 enddo
          
          do i=1,N
          do j=1,M
     !mf= (4.d0*sig(i,j)+ 16.d0*sig(i,j)*(sig(i,j)-0.5d0)**2.d0+ D7*sig(i,j)-D17*sig(i,j))
        ! mff=(4.d0+ 32.d0*sig(i,j)*(sig(i,j)-0.5d0)+16.d0*(sig(i,j)-0.5d0)**2.d0)
          
          
          
          
           eta(i,j)=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/dx**2.d0+  &
         (sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/dy**2.d0
         
         
         Lpzeta(i,j)= (zeta(i+1,j)-2.d0*zeta(i,j)+zeta(i-1,j))/dx**2.d0+  &
         (zeta(i,j+1)-2.d0*zeta(i,j)+zeta(i,j-1))/dy**2.d0
         
             
          enddo
          enddo
         
         do i=1,N
          do j=1,M
         
         Leta(i,j)=(eta(i+1,j)-2.d0*eta(i,j)+eta(i-1,j))/dx**2.d0+  &
         (eta(i,j+1)-2.d0*eta(i,j)+eta(i,j-1))/dy**2.d0
          
           mf= 1.d0/(1.d0-sig(i,j))+ D7*sig(i,j)-D17*sig(i,j)
         mff=1.d0/((1.d0-sig(i,j))**2.d0)
         
         
      asiP(i,j)=-(mf*2.d0/dx**2.d0+2.d0*mf/dy**2.d0+1.d0/dts) - &
      D12*w(i,j)*zeta(i,j) -D14*Leta(i,j) - D10*Lpzeta(i,j)
        
        asiE(i,j)=1.d0*mf/dx**2.d0-u(i,j)*D12/(2.d0*dx)+  &
        (mff+D7-D17)*(sig(i+1,j)-sig(i-1,j))/(4.d0*dx**2.d0)-  &
        D10*(zeta(i+1,j)-zeta(i-1,j))/(4.d0*dx**2.d0)- &
        D14*(eta(i+1,j)-eta(i-1,j))/(4.d0*dx**2.d0)
        
          asiW(i,j)=1.d0*mf/dx**2.d0+u(i,j)*D12/(2.d0*dx)- &
          (mff+D7-D17)*(sig(i+1,j)-sig(i-1,j))/(4.d0*dx**2.d0)+ &
          D10*(zeta(i+1,j)-zeta(i-1,j))/(4.d0*dx**2.d0)+ &
          D14*(eta(i+1,j)-eta(i-1,j))/(4.d0*dx**2.d0)
          
          asiN(i,j)=1.d0*mf/dy**2.d0-v(i,j)*D12/(2.d0*dy)+&
          (mff+D7-D17)*(sig(i,j+1)-sig(i,j-1))/(4.d0*dy**2.d0)- &
          D10*(zeta(i,j+1)-zeta(i,j-1))/(4.d0*dy**2.d0)- &
          D14*(eta(i,j+1)-eta(i,j-1))/(4.d0*dy**2.d0)
          
          asiS(i,j)=1.d0*mf/dy**2.d0+v(i,j)*D12/(2.d0*dy)- &
          (mff+D7-D17)*(sig(i,j+1)-sig(i,j-1))/(4.d0*dy**2.d0)+ &
       D10*(zeta(i,j+1)-zeta(i,j-1))/(4.d0*dy**2.d0)+  &
       D14*(eta(i,j+1)-eta(i,j-1))/(4.d0*dy**2.d0)
          
   
    Ssi(i,j)=-sigp(i,j)/dts 
          enddo
          enddo
         
          
         
          
          
        !sigbc now here   

          do j=1,M
        
         
          sig(0,j)=sig(N,j)

          eta(0,j)=eta(N,j)
          enddo



          do i=1,N
        sig(i,0)=sig(i,M)
        eta(i,0)=eta(i,M)
          enddo





          do j=1,M
          sig(N+1,j)=sig(1,j)
          !sig(N+1,j)=sig(N-1,j)
          eta(N+1,j)=eta(1,j)
          enddo


          do i=1,N
          
          sig(i,M+1)=sig(i,1)
          eta(i,M+1)=eta(i,1)
          enddo
 		
 		
          
          !sigrr=10.d0
		!do while (sigrr>sigtol)
          !sigrr=0.d0
          
         !! Ssi was here
         
          
         !! sig bc was here
         
        

         
  		  do i=1,N
          do j=1,M
          sigo=sig(i,j)
          sig(i,j)=(1.d0-ur)*sig(i,j)+ur*(Ssi(i,j)-asiE(i,j)*sig(i+1,j)-      &
          asiW(i,j)*sig(i-1,j)-asiN(i,j)*sig(i,j+1)-asiS(i,j)*sig(i,j-1))/asiP(i,j)
          if(sig(i,j)>0.999d0)then
          sig(i,j)=0.999d0
          elseif (sig(i,j)<0.001d0)then
          sig(i,j)=0.001d0
          endif
          
         ! sigr(i,j)=abs(sig(i,j)-sigo)
          !if (sigr(i,j)>sigrr)then
         ! sigrr=sigr(i,j)
         ! endif
          enddo
          enddo

!print*, "time step",k, "sig iteration error is",sigrr,"sig",sig(1,10),"sigd",dsig
      ! enddo ! end of sigtol 
       
        do i=1,N
    	 do j=1,M
    	 if(abs(sigold(i,j)-sig(i,j))>err_sig)then
    	 err_sig=abs(sigold(i,j)-sig(i,j))
    	 endif
    	 enddo
    	 enddo
       
       
        enddo         !! end of 2nd sig loop
       
       do i=1,N
    	 do j=1,M
    	 if(abs(oldsig(i,j)-sig(i,j))>err)then
    	 !err=abs(oldsig(i,j)-sig(i,j))
    	 endif
    	 enddo
    	 enddo  
endif ! end of  checking for solving for diffusion equation
       
		 
      
       

         pvti=1+(N-1)/2
         pvtj=1+(M-1)/2
!print*,"time step",k, "itn no",itn,"error is",err,"z",z(pvti,pvtj),"v",u(1,1),"div",divm,&
!"w",w(pvti,pvtj),"Pw",Pw,"tsig",tsig,"dsig",dsig
          !itn=itn+1
 

tsig=0.d0
          do i=1,N
          do j=1,M

          tsig=tsig+dx*dy*sig(i,j)
          enddo
          enddo


if(k>1)then

do i=1,N
do j=1,M
sig(i,j)=sig(i,j)+0.1d0*(totsig-tsig)
if(sig(i,j)>0.999d0)then
sig(i,j)=0.999d0
elseif (sig(i,j)<0.001d0)then
sig(i,j)=0.001d0
endif
enddo
enddo

endif




enddo ! for do while ending

 d2h(1,1)=((z(2,1)-z(1,1))/dx)**2.d0+  &
         ((z(1,2)-z(1,1))/dy)**2.d0
          d2h(N,1)=((z(N,1)-z(N-1,1))/dx)**2.d0+   &
        ((z(N,2)-z(N,1))/dy)**2.d0
          d2h(N,M)=((z(N,M)-z(N-1,M))/dx)**2.d0+   &
         ((z(N,M)-z(N,M-1))/dy)**2.d0
          d2h(1,M)=((z(2,M)-z(1,M))/dx)**2.d0+  &
         ((z(1,M)-z(1,M-1))/dy)**2.d0
          do i=1,N-1
          d2h(i,1)=((z(i+1,1)-z(i-1,1))/(2.d0*dx))**2.d0+ &
         ((z(i,2)-z(i,1))/dy)**2.d0
          d2h(i,M)=((z(i+1,M)-z(i-1,M))/(2*dx))**2.d0+ &
         ((z(i,M)-z(i,M-1))/dy)**2.d0
          enddo
          do j=1,M-1
          d2h(1,j)=((z(2,j)-z(1,j))/(dx))**2.d0+ &
        ((z(1,j+1)-z(1,j-1))/(2.d0*dy))**2.d0
          d2h(N,j)=((z(N,j)-z(N-1,j))/(dx))**2.d0+ &
         ((z(N,j+1)-z(N,j-1))/(2.d0*dy))**2.d0
          enddo
          do i=2,N-1
          do j=2,M-1
          d2h(i,j)=((z(i+1,j)-z(i-1,j))/(2.d0*dx))**2.d0+ &
          ((z(i,j+1)-z(i,j-1))/(2.d0*dy))**2.d0
          enddo
          enddo

 !tsig=0.d0
         




endif  ! close 007

tsig=0.d0
 do i=1,N
 do j=1,M

    tsig=tsig+dx*dy*sig(i,j)
enddo
enddo

dsig=0.d0

do i=1,N
do j=1,M
dsig=dsig+ dx*dy*(sig(i,j)-tsig)**2.d0
enddo
enddo




if(isig==ifl)then !writing file

sigmax=0.d0
do i=1,N
do j=1,M
if(sig(i,j)>sigmax)then
sigmax=sig(i,j)
endif
enddo
enddo

print*,"time step",k, "itn no",itn,"error is",err,"z",z(pvti,pvtj),"v",u(1,1),"div",divm,&
"w",w(pvti,pvtj),"Pw",Pw,"tsig",tsig,"dsig",dsig
   
open(unit=27,file="zeta.txt")
         do i=1,N
         write(27,32)(zeta(i,j),j=1,M)
         enddo
  32     format(1x,5000f16.8)

      open(unit=28,file="sig.txt")
         do i=1,N
         write(28,33)(sig(i,j),j=1,M)
         enddo
  33     format(1x,5000f16.8)

  open(unit=29,file="z.txt")
         do i=1,N
         write(29,34)(z(i,j),j=1,M)
         enddo
  34     format(1x,5000f16.8)

  open(unit=30,file="lam.txt")
         do i=1,N
         write(30,35)(lam(i,j),j=1,M)
         enddo
  35     format(1x,5000f16.8)
  open(unit=31,file="u.txt")
         do i=1,N
         write(31,36)(u(i,j),j=1,M)
         enddo
  36     format(1x,5000f16.8)
  open(unit=32,file="v.txt")
         do i=1,N
         write(32,37)(v(i,j),j=1,M)
         enddo
  37     format(1x,5000f16.8)
  
  open(unit=34,file="fx.txt")
         do i=1,N
         write(34,39)(fx(i,j),j=1,M)
         enddo
  39     format(1x,5000f16.8)
  open(unit=35,file="fy.txt")
         do i=1,N
         write(35,40)(fy(i,j),j=1,M)
         enddo
  40     format(1x,5000f16.8)
  open(unit=36,file="div.txt")
         do i=1,N
         write(36,41)(div(i,j),j=1,M)
         enddo
  41     format(1x,5000f16.8)
  open(unit=56,file="w.txt")
         do i=1,N
         write(56,61)(w(i,j),j=1,M)
         enddo
  61     format(1x,5000f16.8)
  open(unit=57,file="Sla.txt")
         do i=1,N
         write(57,62)(Sla(i,j),j=1,M)
         enddo
  62     format(1x,5000f16.8)
open(unit=69,file="dsig.txt")
         write(69,74)dsig,sigmax 
  74     format(1x,5000f18.12)
  open(unit=33,file="totalsig.txt")
         !do i=1,tnn
         write(33,38)(tsig)
         !enddo
  38     format(1x,5000f16.8)
  open(unit=37,file="t.txt")
         !do i=1,tnn
         write(37,42)(dt*(Pw-1))
         !enddo
  42     format(1x,5000f16.8)
  
  

endif  !! writing file

enddo  ! time loop ends


open(unit=38,file="x.txt")
         do i=1,N
         write(38,43)(x(i))
         enddo
  43     format(1x,5000f16.8)
   open(unit=39,file="y.txt")
         do j=1,M
         write(39,44)(y(j))
         enddo
  44     format(1x,5000f16.8)


    

!print*,"lamP=",laP,"z=",z(10,10),"divm",divm

end program sumingup


!C**********************************************************************
!C Subroutine velfield:
!C     Obtain velocity field spectrally (Hasimoto solution)
!C**********************************************************************

      
      
      


      subroutine velfield(Lx,Ly,nx,ny,ux,uy,fx,fy,S,pp)

      implicit none

      integer nx,ny,nph,ix,iy,nn(2)
      integer ixp,ixm,iyp,iym,i,j

      double precision Lx,Ly,dx,dy
      double precision fx(64,64),fy(64,64),phi(64,64),Scomp(2,64,64)
      double precision fcompx(2,64,64),fcompy(2,64,64),phik(2,64,64),pp(64,64)
      double precision ukx(2,64,64),uky(2,64,64),Scom(2,64,64),pk(2,64,64)
       double precision pkx(2,64,64),pky(2,64,64)
      double precision kx,ky,ksq,fk(2,2),dot(2)
      double precision ux(64,64),uy(64,64),S(64,64),px(64,64),py(64,64)
      double precision pi,cst
      double precision bEW(2,2,64,64),E(2,2),W(2,2),dU(2,2)
      double precision alpha,D,dr,gamma

      

      pi = 4.d0*atan(1.d0)
      dx = Lx/real(nx)
      dy = Ly/real(ny)

!C Copy forces to complex arrays
      do ix = 1,nx
         do iy = 1,ny
            fcompx(1,ix,iy) = fx(ix,iy)
            fcompx(2,ix,iy) = 0.d0
            fcompy(1,ix,iy) = fy(ix,iy)
            fcompy(2,ix,iy) = 0.d0
            Scomp(1,ix,iy) = S(ix,iy)
            Scomp(2,ix,iy) = 0.d0
            
         enddo
      enddo
      nn(1) = nx
      nn(2) = ny

!C Take the Fourier transform of force distribution 
      call fourn(fcompx,nn,2,1)
      call fourn(fcompy,nn,2,1)
      call fourn(Scomp,nn,2,1)

!C Get Fourier coefficients of the velocity field (Hasimoto)
      do iy = 1,ny
         ky = iy-1
         if (iy.gt.(ny/2)) ky = iy-1-ny
         ky = ky/Ly

         do ix = 1,nx
            kx = ix-1
            if (ix.gt.(nx/2)) kx = ix-1-nx
            kx = kx/Lx
            
            if (ix.eq.1.and.iy.eq.1) then
                  
               ukx(1,ix,iy) = 0.d0
               ukx(2,ix,iy) = 0.d0
               uky(1,ix,iy) = 0.d0
               uky(2,ix,iy) = 0.d0
               phik(1,ix,iy)=0.d0
               phik(2,ix,iy)=0.d0
               pkx(1,ix,iy) = 0.d0
               pkx(2,ix,iy) = 0.d0
               pky(1,ix,iy) = 0.d0
               pky(2,ix,iy) = 0.d0
               pk(1,ix,iy) = 0.d0
               pk(2,ix,iy) = 0.d0
               
            else
               
               ksq = kx**2.d0+ky**2.d0

               fk(1,1) = fcompx(1,ix,iy)
               fk(1,2) = fcompy(1,ix,iy)
               fk(2,1) = fcompx(2,ix,iy)
               fk(2,2) = fcompy(2,ix,iy)
                  
               dot(1) = fk(1,1)*kx+fk(1,2)*ky+Scomp(2,ix,iy)*(2.d0*pi*ksq)
               dot(2) = fk(2,1)*kx+fk(2,2)*ky-Scomp(1,ix,iy)*(2.d0*pi*ksq)
                  
               ukx(1,ix,iy) = (fk(1,1)-dot(1)*kx/ksq)/ksq/(4*pi**2.d0)
               ukx(2,ix,iy) = (fk(2,1)-dot(2)*kx/ksq)/ksq/(4*pi**2.d0)
               uky(1,ix,iy) = (fk(1,2)-dot(1)*ky/ksq)/ksq/(4*pi**2.d0)
               uky(2,ix,iy) = (fk(2,2)-dot(2)*ky/ksq)/ksq/(4*pi**2.d0)
              ! phik(1,ix,iy)=(Scomp(1,ix,iy))/(4*pi**2.d0*ksq)
              ! phik(2,ix,iy)=(Scomp(2,ix,iy))/(4*pi**2.d0*ksq)
               pkx(1,ix,iy) = (dot(1)*kx/ksq)
               pkx(2,ix,iy) = (dot(2)*kx/ksq)
               pky(1,ix,iy) = (dot(1)*ky/ksq)
               pky(2,ix,iy) = (dot(2)*ky/ksq)
               pk(1,ix,iy)=-(pkx(2,ix,iy)*kx+pky(2,ix,iy)*ky)/(2*pi*ksq)
               pk(2,ix,iy)=(pkx(1,ix,iy)*kx+pky(1,ix,iy)*ky)/(2*pi*ksq)
               
            endif
            
         enddo
      enddo

!C Call the inverse Fourier transform
      call fourn(ukx,nn,2,-1)
      call fourn(uky,nn,2,-1)
      !call fourn(phik,nn,2,-1)
       call fourn(pkx,nn,2,-1)
      call fourn(pky,nn,2,-1)
      call fourn(pk,nn,2,-1)

!C Rescale arrays
      cst = real(nx*ny)
      do ix = 1,nx
         do iy = 1,ny
            ux(ix,iy) = ukx(1,ix,iy)/cst
            uy(ix,iy) = uky(1,ix,iy)/cst
            !phi(ix,iy)=phik(1,ix,iy)/cst
            px(ix,iy) = pkx(1,ix,iy)/cst
            py(ix,iy) = pky(1,ix,iy)/cst
            pp(ix,iy) = pk(1,ix,iy)/cst
         enddo
      enddo

      return
      end

!C**********************************************************************
!C Subroutine poisson solver:
!C     Obtain solution for Poisson equation (Hasimoto solution)
!C**********************************************************************


subroutine poissol(Lx,Ly,nx,ny,S,phi)

      implicit none

      integer nx,ny,nph,ix,iy,nn(2)
      integer ixp,ixm,iyp,iym,i,j

      double precision Lx,Ly,dx,dy
      double precision Scomp(2,64,64)
      double precision phik(2,64,64)
      double precision kx,ky,ksq,dot(2)
      double precision phi(64,64),S(64,64)
      double precision pi,cst
      

      

      pi = 4.d0*atan(1.d0)
      dx = Lx/real(nx)
      dy = Ly/real(ny)

!C Copy forces to complex arrays
      do ix = 1,nx
         do iy = 1,ny
            Scomp(1,ix,iy) = S(ix,iy)
            Scomp(2,ix,iy) = 0.d0
            
         enddo
      enddo
      nn(1) = nx
      nn(2) = ny

!C Take the Fourier transform of force distribution 
     ! call fourn(fcompx,nn,2,1)
      !call fourn(fcompy,nn,2,1)
      call fourn(Scomp,nn,2,1)

!C Get Fourier coefficients of the velocity field (Hasimoto)
      do iy = 1,ny
         ky = iy-1
         if (iy.gt.(ny/2)) ky = iy-1-ny
         ky = ky/Ly

         do ix = 1,nx
            kx = ix-1
            if (ix.gt.(nx/2)) kx = ix-1-nx
            kx = kx/Lx
            
            if (ix.eq.1.and.iy.eq.1) then
                  
               phik(1,ix,iy)=0.d0
               phik(2,ix,iy)=0.d0
               
            else
               
               ksq = kx**2.d0+ky**2.d0

               
               phik(1,ix,iy)=-(Scomp(1,ix,iy))/(4*pi**2.d0*ksq)
               phik(2,ix,iy)=-(Scomp(2,ix,iy))/(4*pi**2.d0*ksq)
               
               
            endif
            
         enddo
      enddo

!C Call the inverse Fourier transform
      
      call fourn(phik,nn,2,-1)
       

!C Rescale arrays
      cst = real(nx*ny)
      do ix = 1,nx
         do iy = 1,ny
            
            phi(ix,iy)=phik(1,ix,iy)/cst
            
         enddo
      enddo

      return
      end



!C*******************************************************************
!C Function fourn:
!C     Computes the fast Fourier transform of a multidimensional
!C     array of complex numbers.
!C     From Numerical Recipes, Press et al., 2001.
!C*******************************************************************

      SUBROUTINE fourn(data,nn,ndim,isign) 

      INTEGER isign,ndim,nn(ndim) 
      double precision data(*) 
      INTEGER i1,i2,i2rev,i3,i3rev,ibit,idim,ifp1,ifp2,ip1,ip2 
      INTEGER ip3,k1,k2,n,nprev,nrem,ntot 
      double precision tempi,tempr 
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp 

      ntot=1 
      do idim=1,ndim 
         ntot=ntot*nn(idim) 
      enddo 
      nprev=1 
      do idim=1,ndim 
         n=nn(idim) 
         nrem=ntot/(n*nprev) 
         ip1=2*nprev 
         ip2=ip1*n 
         ip3=ip2*nrem 
         i2rev=1 
         do i2=1,ip2,ip1 
            if(i2.lt.i2rev) then 
               do i1=i2,i2+ip1-2,2 
                  do i3=i1,ip3,ip2 
                     i3rev=i2rev+i3-i2 
                     tempr=data(i3) 
                     tempi=data(i3+1) 
                     data(i3)=data(i3rev) 
                     data(i3+1)=data(i3rev+1) 
                     data(i3rev)=tempr 
                     data(i3rev+1)=tempi 
                  enddo 
               enddo 
            endif 
            ibit=ip2/2 
 1          if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then 
               i2rev=i2rev-ibit 
               ibit=ibit/2 
               goto 1 
            endif 
            i2rev=i2rev+ibit 
         enddo 
         ifp1=ip1  
 2       if(ifp1.lt.ip2) then 
            ifp2=2*ifp1 
            theta=isign*6.28318530717959d0/(ifp2/ip1) 
            wpr=-2.d0*sin(0.5d0*theta)**2 
            wpi=sin(theta) 
            wr=1.d0 
            wi=0.d0 
            do i3=1,ifp1,ip1 
               do i1=i3,i3+ip1-2,2 
                  do i2=i1,ip3,ifp2 
                     k1=i2 
                     k2=k1+ifp1 
                     tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1) 
                     tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2) 
                     data(k2)=data(k1)-tempr 
                     data(k2+1)=data(k1+1)-tempi
                     data(k1)=data(k1)+tempr 
                     data(k1+1)=data(k1+1)+tempi 
                  enddo 
               enddo 
               wtemp=wr  
               wr=wr*wpr-wi*wpi+wr 
               wi=wi*wpr+wtemp*wpi+wi 
            enddo 
            ifp1=ifp2 
            goto 2 
         endif 
         nprev=n*nprev 
      enddo 
      return 
      END


