!-----------------------------------------------------------------------------------------------------------------
! Cahn-Hilliard model of aggregation diffusion of proteins on flat membrane
!
! Simulation code by: Arijit Mahapatra, David Saintillan, Padmini Ranganami (UCSD)
!
! For details of the model, refer to:
! "Curvature driven feedback on aggregation-diffusion of proteins in lipid bilayers" by Mahapatra et al.
! (arxiv_link)
!-----------------------------------------------------------------------------------------------------------------
program sumingup
implicit none
integer:: N,M,P,i,j,k,ptc,itn,lamit,pvti,pvtj,velit,ii,jj,isig,iisig,Pw,tnn
double precision:: sig0, lam0, alpha, muph, kp, beta,c,err,Tst,Ted,Tstp,stth
double precision:: pr, L, D1, D2, D3, D4, D5, D6, D7, D8, D9
double precision:: D10,lr,tT,lX,lY,dx,dy,dt,Dmf,xa,xb,ya,yb
double precision:: cth,lct,tol,sigtol,ur,zr,zrr,sigrr
double precision:: zer,zerr,lamrr,lamr,rZe,zenu,znu,vrr,sigmax
double precision:: signu,dzetae,dzetaw,dzetan,dzetas,zcon,zecon,velerr,utol
double precision:: zetadx,zetady,azE,azW,azN,azS,azP,uer,ver,xd,yd,totsig,dsig
double precision:: sigdx,sigdy,Lsig,Lsigsq,sigo,uW,uE,uN,uS,urr,uP,sigsqdy,D11
double precision:: zeo,zo,azeE,azeW,azeN,azeS,sigsqdx,dvy,dvx,duy,dux,d2zy,d2zxy,d2zx
double precision:: dsign,dsige,dsigw,dsigs,umean,vmean,usum,vsum,xt,yt,r,PI,ddx,ddy 
double precision:: laW,laS,laP,laN,lamo,laE,lamnu,zera,uup,vup,nu,uo,vo,delta,divm,xc,yc
double precision:: gamma, sigs, D12, Dag, sigrt,zetarr,azeP,D14,D17,mf,mff,kBT,D19
double precision, dimension(0:500):: x,y,siggs,sigge,siggw,siggn,sigh,sigv,zetah,zetav
double precision, dimension(0:5000):: t,tsig,Je,Jw,Jn,Js,sigd,max_sig
double precision, dimension(0:500,0:500):: sigr,Sze,d2h,Sz,Ssi,Sla,shst,uold,vold
double precision, dimension(0:500,0:500):: zeta,z,sig,lam,sigp,u,v,Sv,Su,zold,zetaold
double precision, dimension(0:500,0:500)::asiE,asiW,asiN,asiS,g1,g2,fx,fy,div,asiP
double precision, dimension(0:500,0:500)::dzetax,dzetay,dsigx,dsigy,Lpzeta
double precision, dimension(0:500,0:500):: oldu,oldv,oldlam,oldsig ,oldzeta
double precision rdm(64,64)
!sig0=5d0
!lam0=5d-2
!alpha=1.d0
!muph=1d-3
!kp=500.d0
CALL RANDOM_NUMBER(rdm)




PI=2.D0*DASIN(1.D0)

!!Specify Grid numbers, change the size of the variable defined in line 32 along with 
!!changing grid number 

N=64
M=64
tnn=2001  !! timestep to print result
iisig=30000  !! internal timestep
P=(tnn-1)*iisig+1
lr=1.d0
tT=3d-1  !! time of the simulation
lX=1.d0
lY=lX/lr
dx=lX/(N)
dy=lY/(M)

!kBT=4.14d0

Pw=0
!D14=0.d0
!! Parameters====================>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>///\\\\>>>>>>>>>>>>>>>>>>
dt=tT/(P-1)
D17=25.d0!  hatA
D19=200  ! hatS

!! derived ND number
D14= D17/(2.d0*D19) !.001d0   !D17/(2.d0*sigs*L**2.d0)


y(1)=-lY/2.d0
do j=2,M
y(j)=y(j-1)+dy
enddo
x(1)=-lX/2.d0
do i=2,N
x(i)=x(i-1)+dx
enddo
!dt=tT/(P-1)

t(1)=0.d0
do k=2,P
!t(k)=t(k-1)+dt
enddo


do i=0,N+1
do j=0,M+1
sig(i,j)=0.d0
enddo
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
    Je(k)=0.d0
    Jw(k)=0.d0 
    Jn(k)=0.d0
    Js(k)=0.d0  
enddo



cth=20.d0
lct=5d-2
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

!!! ==========initial_condition_for_sigma==============================
do i=1,N
do j=1,M
!sig(i,j)=sig0*(0.5d0*(1.d0-tanh(cth*(((x(i)-xa)**2.d0+(y(j)-ya)**2.d0)**0.5d0-lct))) + &
!0.5d0*(1.d0-tanh(cth*(((x(i)-xb)**2.d0+(y(j)-yb)**2.d0)**0.5d0-lct)))  + &
!0.5d0*(1.d0-tanh(cth*(((x(i)-xc)**2.d0+(y(j)-yc)**2.d0)**0.5d0-lct))) + &
!0.5d0*(1.d0-tanh(cth*(((x(i)-xd)**2.d0+(y(j)-yd)**2.d0)**0.5d0-lct))))

!sig(i,j)=1.d0-tanh(9.d0*((x(i)**2.d0+y(j)**2.d0)**0.5d0-2d-1))
!sig(i,j)=(tanh(cth*(x(i)+lct))-tanh(cth*(x(i)-lct)))*(tanh(cth*(y(j)+lct))- &
!tanh(cth*(y(j)-lct)))
!sig(i,j)=sig0+ 0.001*rand()
sig(i,j)=0.1d0 +0.0001d0*rdm(i,j)
enddo
enddo

!! to start with 
do i=1,N
sig(i,0)=sig(i,M)
sig(i,M+1)=sig(i,1)
enddo
do j=1,M
sig(0,j)=sig(N,j)
sig(N+1,j)=sig(1,j)
enddo

do i=1,N
do j=1,M
zeta(i,j)=(sig(i+1,j)-2*sig(i,j)+sig(i-1,j))/(dx**2.d0) +  &
(sig(i,j+1)-2*sig(i,j)+sig(i,j-1))/(dy**2.d0) 
enddo
enddo
totsig=0.d0
 dsig=0.d0
          do i=1,N
          do j=1,M
          

          
          totsig=totsig+dx*dy*sig(i,j)
          
          enddo
          enddo
          
          do i=1,N
          do j=1,M
          dsig=dsig+((totsig-sig(i,j)))**2.d0*dx*dy
          enddo
          enddo

tol=1d-8
sigtol=5d-8
!utol=1d-10
ur=7d-1
!! coefficients
 		
		
        azeP=-(2.d0/dx**2.d0+2.d0/dy**2.d0)

isig=iisig
do k=1,P
if (k>1) then   ! open 007          ! just to write file
  	if(k>1)then
  	do i=1,N
  	do j=1,M
	sigp(i,j)=sig(i,j)
	 enddo
	 enddo
	endif
	itn=1
	err=5.d0
        do while (err>tol .OR. itn<5)
         !do itn=1,20
        itn=itn+1
  err=0.d0
 !! =====================Iteration for protein density starts here===============
!sig iteration starts==========================================================
!==============================================================================
if (k>1)then



         do i=1,N
  		 do j=1,M
  		 oldsig(i,j)=sig(i,j)
  		 oldzeta(i,j)=zeta(i,j)
  		 enddo
  		 enddo
          
          
         
		
		
		
		zetarr=0.d0
          
          
           do i=1,N
          do j=1,M
          
       mf= 1.d0/(1.d0-sig(i,j))-D17*sig(i,j)
        mff=1.d0/(1.d0-sig(i,j))**2.d0
          
          dsigx(i,j)=(sig(i+1,j)-sig(i-1,j))/(2.d0*dx)
          dsigy(i,j)=(sig(i,j+1)-sig(i,j-1))/(2.d0*dx)
          
          dzetax(i,j)=(zeta(i+1,j)-zeta(i-1,j))/(2.d0*dx)
          dzetay(i,j)=(zeta(i,j+1)-zeta(i,j-1))/(2.d0*dx)
          
          Lpzeta(i,j)=(zeta(i+1,j)-2*zeta(i,j)+zeta(i-1,j))/(dx**2.d0) +  &
        (zeta(i,j+1)-2*zeta(i,j)+zeta(i,j-1))/(dy**2.d0) 
          
          asiE(i,j)=mf*1.d0/dx**2.d0 + ((mff-D17)*dsigx(i,j)-D14*dzetax(i,j))/(2.d0*dx)
          asiW(i,j)=mf*1.d0/dx**2.d0- ((mff-D17)*dsigx(i,j)-D14*dzetax(i,j))/(2.d0*dx)
          asiN(i,j)=mf*1.d0/dy**2.d0+ ((mff-D17)*dsigy(i,j)-D14*dzetay(i,j))/(2.d0*dy)
          asiS(i,j)=mf*1.d0/dy**2.d0-((mff-D17)*dsigy(i,j)-D14*dzetay(i,j))/(2.d0*dy)
      
           asiP(i,j)=-1.d0/dt-(2.d0/dx**2.d0+2.d0/dy**2.d0)*mf -D14*Lpzeta(i,j)
                     
          enddo
          enddo
          
           do i=1,N
          do j=1,M
        Ssi(i,j)=-sigp(i,j)/dt
          enddo
          enddo
          
        
        
		
          
          
         

         
        
          !sigrr=10.d0
          
		!do while (sigrr>sigtol)
		
		
			
		
		
          sigrr=0.d0

          do i=1,N
		do j=1,M
 zeta(i,j)=(sig(i+1,j)-2.d0*sig(i,j)+sig(i-1,j))/(dx**2.d0) +  &
(sig(i,j+1)-2.d0*sig(i,j)+sig(i,j-1))/(dy**2.d0) 

		enddo
		enddo

          
 		
  		  do i=1,N
          do j=1,M
          sigo=sig(i,j)
          sig(i,j)=(1.d0-ur)*sig(i,j)+ur*(Ssi(i,j)-asiE(i,j)*sig(i+1,j)-      &
          asiW(i,j)*sig(i-1,j)-asiN(i,j)*sig(i,j+1)-asiS(i,j)*sig(i,j-1))/asiP(i,j)
          sigr(i,j)=abs(sig(i,j)-sigo)
          
          if(sig(i,j)>0.999d0)then
       sig(i,j)=0.999d0
      elseif (sig(i,j)<0.001d0)then
      sig(i,j)=0.001d0
      endif
          
          
          
        !!!=====================Implementation_of_periodic_bc====================  
         ! if (sigr(i,j)>sigrr)then
          !sigrr=sigr(i,j)
          !endif
          enddo
          enddo
          
          
           do j=1,M
          dsigw=0.d0
          siggw(j)=sig(2,j)-2.d0*dx*dsigw
          !sig(0,j)=siggw(j)
          !sigv(j)=sig(0,j)
          sig(0,j)=sig(N,j)
          enddo



          do i=1,N
          dsigs=0.d0
          siggs(i)=sig(i,2)-2.d0*dy*dsigs
          !sig(i,0)=siggs(i)
          !sigh(i)=sig(i,0)
          sig(i,0)=sig(i,M)
          enddo





          do j=1,M
          dsige=0.d0
          sigge(j)=sig(N-1,j)+2.d0*dx*dsige
          !sig(N+1,j)=sigge(j)
          sig(N+1,j)=sig(1,j)
          enddo


          do i=1,N
          dsign=0.d0
          siggn(i)=sig(i,M-1)+2.d0*dy*dsign
          !sig(i,M+1)=siggn(i)
          sig(i,M+1)=sig(i,1)
          enddo
          
          
          do j=1,M
		!zetav(j)=zeta(0,j)
		zeta(0,j)=zeta(N,j)
		zeta(N+1,j)=zeta(1,j)
		enddo
		
		do i=1,N
		!zetah(i)=zeta(i,0)
		zeta(i,0)=zeta(i,M)
		zeta(i,M+1)=zeta(i,1)
		enddo
          
          
          
      !if (zetarr>sigrr)then
      !sigrr=zetarr
     ! endif
      
!print*, "time step",k, "sig iteration error is",sigrr

       
       !enddo ! end of sigtol 
       
endif ! end of  checking for solving for diffusion equation

		 do i=1,N
    	 do j=1,M
    	 if(abs(oldsig(i,j)-sig(i,j))>err)then
    	 err=abs(oldsig(i,j)-sig(i,j))
    	 endif
    	 enddo
    	 enddo  
        do i=1,N
    	 do j=1,M
    	 if(abs(oldzeta(i,j)-zeta(i,j))>err)then
    	 !err=abs(oldzeta(i,j)-zeta(i,j))
    	 endif
    	 enddo
    	 enddo  


         pvti=1+(N-1)/2
         pvtj=1+(M-1)/2
!print*,"time step",k, "itn no",itn,"error is",err,"sig",sig(pvti,pvtj),"totalsig",totsig,"sigd",dsig
          !itn=itn+1
 

enddo ! for do while ending

 
 totsig=0.d0
 dsig=0.d0
          do i=1,N
          do j=1,M
          
          totsig=totsig+dx*dy*sig(i,j)
          
          enddo
          enddo
          
          do i=1,N
          do j=1,M
          dsig=dsig+((tsig(1)-sig(i,j)))**2.d0*dx*dy
          enddo
          enddo
          

 
   do i=1,N
  do j=1,M
sig(i,j)=sig(i,j)+(tsig(1)-totsig)/1.d0
if(sig(i,j)>0.999d0)then
sig(i,j)=0.999d0
elseif (sig(i,j)<0.001d0)then
sig(i,j)=0.001d0
endif
enddo
enddo
 




endif  ! close 007


   
if(isig==iisig)then

sigmax=0.d0
do i=1,N
do j=1,M
if(sig(i,j)>sigmax)then
sigmax=sig(i,j)
endif
enddo
enddo


print*,"time step",Pw, "itn no",itn,"error is",err,"sig",sig(pvti,pvtj),"totalsig",totsig,"sigd",dsig,"sigmax",sigmax
      open(unit=28,file="sig.txt")
         do i=1,N
         write(28,33)(sig(i,j),j=1,M)
         enddo
  33     format(1x,5000f16.8)

isig=0

Pw=Pw+1
tsig(Pw)=totsig
sigd(Pw)=dsig
max_sig(Pw)=sigmax
 endif
 
 isig=isig+1




enddo  ! time loop ends

open(unit=33,file="totalsig.txt")
         do i=1,Pw
         write(33,38)(tsig(i))
         enddo
  38     format(1x,5000f16.8)
  open(unit=37,file="t.txt")
         do i=1,Pw
         write(37,42)((i-1)*iisig*dt)
         enddo
  42     format(1x,5000f16.8)
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

open(unit=40,file="sigd.txt")
         do i=1,Pw
         write(40,45)sigd(i),max_sig(i) 
         enddo
  45     format(1x,5000f18.12)

    

!print*,"lamP=",laP,"z=",z(10,10),"divm",divm

end program sumingup
