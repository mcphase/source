!  Program chibfk(2404).f95 calculates dynamical single ion susceptibility for RE ions
!  Uses as input the spin-eigenstates of the ion,   
 
module CommonData
! This module contains data used by all other subroutines 
implicit none
save 
integer :: Ns,Ms !number of states, number of dynamical variables
integer, parameter :: Np=20 !max number of states
integer, parameter :: Mp=100 !max number of transitions
integer, parameter :: Nc=250 !max number of scatt. vect. 
!integer :: nlg,nlf
real :: gl !  Lande factor
real :: diff(3,3), sf(3,3), kv(3), ksv(3), kapv(3), kapnv(3)
complex :: chi(3,3)
!real :: kapg(Nc),kapf(Nc),strg(Nc),strf(Nc)
real :: E0, Es, k1,k2,k3,ks1,ks2,ks3,kap1,kap2,kap3,kapv1,kapv2,kapv3
real :: enloss, kappa, Emax, Emin, epsilon
real :: jav(3)
integer :: mode,mst !scattering mode, file type
integer :: Npoints=200
real :: beta, g, gam, temp ! 1/k_B T, coupl. const 
real, parameter :: Pi=3.14159, ek2=2.072, meVkT=11.6, cutoff=100., r0=-0.54 
! exp. cut-off
complex, parameter :: Iunit=(0.0,1.0)
character(len=20) :: outfilename
character(len=20) :: cefname, enfilename
end module CommonData
!--------------------------------------------------------------------------

module MatrixElements
use CommonData
implicit none
save
integer :: v1(Mp), v2(Mp) !states belonging to one transition
real :: En(Np), p(Np), Pp(Mp) !energy eigen values, Boltzmann factors, static 
! susceptibilities

complex :: Ev(Np,Np),jjj(3,Np,Np) !eigenstates, spin matrixelements, magnetisation
complex :: jjx(Np,Np),jjy(Np,Np), jjz(Np,Np), jjp(Np,Np), jjm(Np,Np)
end module MatrixElements
!---------------------------------------------------------------
subroutine DataRead
!  Reads input-data from command line and input-files  
!  It asks for the coupling with conduction electrons g=J_ex N(0), 
!  temperature T and name of input-file containing 
!  Eigen-energies En(n) and eigen-states Ev(n,k) of RE ion in 
!  a crystalline electric field.
!  The fourth entry should be the name of the file containing information 
!  about the energy range. If this information is missing, 
!  it can be put in from the screen 

use CommonData
use MatrixElements
implicit none
!save
integer :: i,j, k,l, ii, status_read,n,nn,m,cline, vv1(Mp),vv2(Mp) 
integer :: zz,ww, type_of_contence,linetype(50), linenumber, nlines
real :: s,dk,kn,knn,ksn,ks,k0, Evr(Np,Np),Evi(Np,Np), v(Np),f, x, y 
real ::  x1,x2,x3,y1,y2,y3,strfak
complex ::  cs
character(len=20)::  arg 
character(len=20) :: name(5)
character(len=10) :: couplc, tempp, mod,filet
character(len=150) :: line(150),tline(150),bb, bbd ,zeile
character(len=20) :: redata(12),re(12), number


!  read information from command-line   

do i=1, iargc() ! number of entries in commandline
    call getarg(i,arg) ! entries in command line
    name(i)= arg
end do
n=iargc() !number of entries
!check if command line contains information 
if (n .le. 3) then
  write(*,*) 'command line entries incomplete'
  write(*,*) 'should contain coupling constant, temperature,' 
  write(*,*) 'name of file with crystal field levels and file' 
  write(*,*) 'containing energy range of the susceptibility '
  stop
end if
couplc=name(1)
read(couplc,*) x
g=x
tempp=name(2)
read(tempp,*) y
temp=y
cefname=name(3)
enfilename=name(4)



! Analyse file with cef-data
open(15,file=cefname, action='read')
redata(1)= ' 6   0.857   Ce '
redata(2)= ' 9   0.800   Pr '
redata(3)= ' 10  0.727   Nd '
redata(4)= ' 9   0.600   Pm '
redata(5)= ' 6   0.286   Sm '
redata(6)= ' 8   2.000   Gd '
redata(7)= ' 15  1.500   Tb '
redata(8)= ' 16  1.333   Dy '
redata(9)= ' 17  1.250   Ho '
redata(10)= ' 16  1.200   Er '
redata(11)= ' 13  1.169   Tb '
redata(12)= ' 8   1.143   Yb '
re(1)='Ce'
re(2)='Pr'
re(3)='Nd'
re(4)='Pm'
re(5)='Sm'
re(6)='Gd'
re(7)='Tb'
re(8)='Dy'
re(9)='Ho'
re(10)='Er'
re(11)='Tb'
re(12)='Yb'
n=1
do 
  read(15,'(A)',iostat=status_read) line(n)
  if(status_read /= 0) exit
  n=n+1
end do 
nlines=n-1
m=0
do n=1,nlines
  if (line(n)(1:1)=='#') then
    if (line(n)(1:2)=='#!') then
      bb=line(n)(3:150)
      nn=0
      do while (bb(1:1)==' ')
        bbd=bb(2:150)
        bb=bbd
      end do      
      if (bb(1:1)== 'd') then
        m=m+1        
        bbd=bb(2:150)
        nn=index(bbd,'=')
        tline(m)=bbd(nn+1:nn+3)
      end if
      if (bb(1:11)=='Eigenvalues') then
        m=m+1
        bbd=bb(12:150)
        nn=index(bbd,'=')
        tline(m)=bbd(nn+1:150)
      end if
    end if
  else 
    m=m+1
    tline(m)=line(n)
  end if
end do
open(16,file='cefworkfile.dat')  
do n=1,m
  write(16,'(A)') tline(n)
end do
close(16)
close(15)

!Analyse file with energy range (parfile.par)
open(17,file=enfilename,action='read')
n=1 
do 
  read(17,'(A)',iostat=status_read) line(n)
  if(status_read /= 0) exit
  n=n+1
end do 
nlines=n-1
m=0
do n=1,nlines
  if (line(n)(1:7)=='#!emin=') then
    number=line(n)(8:30)
    read(number,*) emin
  end if
  if (line(n)(1:7)=='#!emax=') then
    number=line(n)(8:30)
    read(number,*) emax
  end if
  if (line(n)(1:10)=='#!epsilon=') then
    number=line(n)(11:30)
    read(number,*) epsilon
  end if
end do
close(17)

beta = 11.6/temp
gam= 2*Pi*g**2

!open  workfile 
  open(16,file='cefworkfile.dat',action='read')
  read(16,*) Ns
  read(16,*) (En(n),n=1,Ns)
   do j=1,Ns
     read(16,*) (Evr(j,i), i=1,Ns)
   end do  
  do j=1,Ns
    read(16,*) (Evi(j,i), i=1,Ns)
  end do  
  close(16)
  do i=1,Ns
    do j=1,Ns
      Ev(j,i)=cmplx(Evr(j,i),Evi(j,i)) 
    end do
  end do

end subroutine DataRead

!--------------------------------------------------------

subroutine MatrixElementCalculation

use CommonData
use MatrixElements
implicit none
!save
integer :: i,j,k,n, nn,m,l, mm,n1,n2,it,vv1(Mp),vv2(Mp),M0
real :: Jj, mj, q, s,z, ss
real :: jp(Np,Np),jm(Np,Np),jz(Np,Np), qq(Np,Np), Ps(Mp)
complex :: cs, csz,csp,csm
! quantities used, which are defined in CommonData or MatrixElements
! Ns(number of states), Np(max number of states, 
! Ms(number of transitions considered), Mp max number of transitions
! En(Np)(energy levels), Ev(Np,Np)(complex eigenvectors) 
! beta(inv temperature), cutoff(exp cutoff)
! check of eigenvectors for orthogonality:
open(19,file='jmatrix.dat')

do i = 1,Ns
do j = 1,Ns
cs=(0.,0.)
do k = 1,Ns
cs=cs+conjg(Ev(j,k))*Ev(i,k)
end do
write(19,*) i,j,cs
if ((cabs((cs-1)*cs)>= 0.001) )then
write(*,*) 'low accuracy of eigen vectors' 
end if
end do
end do

!the file jmatrix is generated only for checking the correctnes of  the 
!matrix elements
! Calculation of spin matrix elements
Jj=(Ns-1.)/2. !j-value
do i = 1,Ns
do k = 1,Ns
jz(i,k)=0.
jp(i,k)=0.
jm(i,k)=0.
end do 
end do
do i = 1,Ns
mj=i-1-Jj
jz(i,i)=mj
end do
do i  =1,Ns-1
mj=i-1-Jj
jp(i+1,i)=sqrt((Jj-mj)*(Jj+mj+1))
end do
do i = 2,Ns
mj=i-1-Jj
jm(i-1,i)=sqrt((Jj+mj)*(Jj-mj+1))
end do
write(19,*) 'real(jjx(n,m)), aimag(jjy(n,m)), real(jjz(n,m))'
write(19,*) 'in colums'
write(19,*) '------------------------'

! in the following the matrix-elements of the total 
! angular momentum J are calculated 
do m = 1,Ns
do n = 1,Ns
csz=(0.,0.)
csp=(0.,0.)
csm=(0.,0.)
do i=1,Ns
do k=1,Ns
csz=csz+conjg(Ev(i,n))*jz(i,k)*Ev(k,m)
csp=csp+conjg(Ev(i,n))*jp(i,k)*Ev(k,m)
csm=csm+conjg(Ev(i,n))*jm(i,k)*Ev(k,m)
end do
end do
jjz(n,m)=csz !stored in MatrixElements
jjp(n,m)=csp !stored in MatrixElements
jjm(n,m)=csm !stored in MatrixElements
jjx(n,m)=(csp+csm)/2
jjy(n,m)=(csp-csm)/(2*Iunit)
write(19,*) real(jjx(n,m)), aimag(jjy(n,m)), real(jjz(n,m))
end do
write(19,*) '------------------------'
end do

do n=1,Ns
do m=1,Ns
jjj(1,n,m)=jjx(n,m) !stored in CommonData
jjj(2,n,m)=jjy(n,m) !stored in CommonData
jjj(3,n,m)=jjz(n,m) !stored in CommonData
end do
end do 

do i=1,Ns
do j=1,Ns

ss=0.
do l=1,Ns
ss=ss+ jjx(i,l)*jjx(l,j) + jjy(i,l)*jjy(l,j)+ jjz(i,l)*jjz(l,j)
end do
qq(i,j)=ss
end do
end do
write(19,*) 'Check matrix elements of JJ'
do i=1,Ns
write(19,'(10E12.3)' ) (qq(i,j), j=1,Ns)
end do

! Boltzmann factors p(n)
do n=1, Ns
s=beta*En(n)
! the Boltzmann factors are set to 0, if exponents become
! too large 
! with negative sign 
if (s > cutoff) then
p(n)=0 
else  
p(n)=exp(-beta*En(n))
end if
end do
Z=0
do n=1,Ns 
Z=Z+p(n)
end do
do n=1,Ns
p(n)=p(n)/Z
end do

!Calculate magnetisation
cs=0
do i=1,3
cs=0.
do n=1,Ns
cs=cs+jjj(i,n,n)*p(n)
end do
jav(i)=real(cs)
end do
do n=1,Ns
do i=1,3
jjj(i,n,n)=jjj(i,n,n)-jav(i)
end do
end do

write(19,'(6E12.2)') (p(n), n=1,Ns)

!  Static susceptibilities

write(19,*) 'Magnetisation'
write(19,*)  (jav(i), i=1,3)
! Selection of dynamical variables
! for each transition mm between states n1 and n2 the two states 
! are stored in the fields vv1 and vv2. There are M0=Ns x Ns transitions
mm=0
do n=1,Ns
do m =1,Ns
mm=mm+1
vv1(mm)=n
vv2(mm)=m
end do
end do
M0=mm
mm=0
do nn=1,M0
  ss=0
  n1=vv1(nn)
  n2=vv2(nn)
  if (En(n1) == En(n2)) then
      s=p(n1)
  else
      s=(p(n2)-p(n1))/(beta*(En(n1)-En(n2)))
  end if
  ss=ss+s
  Ps(nn)=ss
  write(19,*) 'Ps',  nn, Ps(nn), vv1(nn), vv2(nn)
! if Boltzmann factors are zero, the transtions are eliminated
  if (p(n1)< 0.000001.and. p(n2)<0.000001)  then
    vv1(nn)=0
    vv2(nn)=0
  end if
end  do


mm=0
do nn=1,M0
if (vv1(nn)/= 0 .and. vv2(nn)/=0) then
mm=mm+1
v1(mm)=vv1(nn)
v2(mm)=vv2(nn)
Pp(mm)=Ps(nn)
end if
end do
Ms=mm
write(*,*) 'number of dynamical variables', Ms
do mm=1,Ms
write(19,*) mm, pp(mm), v1(mm),v2(mm)
end do
close(19)
end subroutine MatrixElementCalculation

!----------------------------------------------------
subroutine MatInv(Min,Mout,N,Np)
!  Min originale Matrix der Dimension NxN, input
!  Mout invertierte Matrix der Dimension NxN, output
!  N aktuelle Dimension der zur invertierende Matrix
!  Np physikalische  Dimension der Matrizen
implicit none
integer, parameter :: Nmax=100
complex :: A(Np,Np), B(Np)
complex, intent( in) :: Min(Np,Np)
complex, intent( out ) :: Mout(Np,Np)
integer, intent( in ) ::  N, Np
integer :: indx(N)
real :: Vv(Nmax), Y(Np,Np)
real :: Aamax, c, dum, d
complex ::  sum, dsum
integer :: i,ii, j, k,l, ll, m, Imax
do i=1,N 
   do j=1,N
      A(i,j)=Min(i,j)
   end do
end do
d=1
do i=1,N
   Aamax=0
   do j=1,N 
      c=cabs(A(i,j))
      if (c > Aamax) then
         Aamax=c
      end if
   end do
   if (Aamax == 0.) then
        write(*,*) 'Matrix singulaer'
   end if

   Vv(i)=1./Aamax
end do
do j = 1,N
   do i =1,j-1
      sum=A(i,j)
      do k=1,i-1     
         sum=sum-A(i,k)*A(k,j)
      end do
      A(i,j)=sum
   end do
   Aamax=0.
   do i=j,N
      sum=A(i,j)
      do k=1,j-1
         sum=sum - A(i,k)*A(k,j)
      end do
      A(i,j)=sum
      dum= Vv(i)*cabs(sum)
      if (dum >= Aamax) then
         Imax=i
         Aamax=dum
      end if
   end do
   if (j /= Imax) then
      do k=1,N
         dsum=A(Imax,k)
         A(Imax,k)=A(j,k)
         A(j,k)=dsum
      end do
      d=-d
      Vv(Imax)=Vv(j)
   end if
   Indx(j)=Imax
   if (j /= N) then
      dsum= 1./A(j,j)
      do i=j+1,N
         A(i,j)=A(i,j)*dsum
      end do
   end if
end do
do i=1,N  
   do j =1,N
      Y(i,j)=0.
   end do
   Y(i,i)=1;
end do
do l=1,N
   do m =1,N
      B(m)=Y(l,m)
   end do
   ii=0
   do i=1,N  
      ll=Indx(i)
      sum=B(ll)
      B(ll)=B(i)
      if(ii/=0) then
         do j = ii,i-1
            sum=sum-A(i,j)*B(j)
         end do
      else if (sum /= 0 ) then
         ii=i
      end if
      B(i)=sum
   end do
   do i=N,1,-1
      sum=B(i) 
      do j=i+1,N
         sum=sum-A(i,j)*B(j)
      end do
      B(i)=sum/A(i,i)
   end do
   do m=1,N
       Mout(m,l)=B(m)
   end do
end do
return
end subroutine MatInv

!-------------------------------------------------------------------------- 

subroutine SuscepComponents(x)  
!calculates the memory function matrix for the spin-operators 
!as function of frequency x and the different spin components
!of the dynamic susceptibility/(1-exp(-beta*omega))
use CommonData
use MatrixElements
implicit none
!save
real :: x,y,zx,s,F(Np,Np),cut
complex :: mat, cs, css
integer :: t,n,m,i,j,k,l,nn,mm,n1,n2,m1,m2
complex :: Phi(Mp,Mp), Mem(Mp,Mp),Om(Mp,Mp),Ominv(Mp,Mp)
real :: Delta,s0,s1,s2
! the following quantities used in this procedure are defined in CommonData
! or MatrixElements:
! Ns(number of states), Np(max number od states),
! Ms(number of transitions), Mp(max. number of transitions),
! En(n)(energy levels), p(n) (Boltzmann factors), Pp(static susceptibilities)
! beta(inverse temperature) ,gam(coupling constant, gam=2*Pi*g^2)
! v1(nu), v2(nu)(states belonging to transition nu)   
! jjz(n1,n2), jjp(n1,n2), jjm(n1,n2) (spin matrix elements)
! enloss 
! dynamical part of the memory function
interface
   subroutine Matinv(Min,Mout,Ms,Mp)
   complex, intent(in) :: Min(Mp,Mp)
   complex, intent(out) :: Mout(Mp,Mp)
   integer :: Ms,Mp
   end subroutine Matinv 
end interface

cut=cutoff/2
do n=1,Ns
  do m=1,Ns
    Delta = En(n)-En(m)
    s0=sqrt(p(n)*p(m))/beta
    y=beta*Delta 
    zx=beta*x
    if (Delta == 0) then
      F(n,m) = p(m)/beta 
    else if (abs(zx) <= 0.0001) then
      if (y>cut) then
        F(n,m)=p(n)*Delta
      else if (-y > cut) then
        F(n,m)= -p(m)*Delta
      else 
        F(n,m)= s0*y/(exp(y/2) - exp(-y/2)) 
      end if
    else if(x == Delta) then
      if (y-zx>cut) then
        F(n,m)= p(m)/beta/y
      else if (-(y-zx)>cut) then        
        F(n,m)= -p(n)/beta/y
      else
        F(n,m)= s0/y*(p(m)-p(n))
      end if
    else 
      if ((zx> cut) .and. (zx-y > cut) ) then
        F(n,m)=p(m)/beta*(x-Delta)/x
      else if ((zx > cut).and. (abs(zx-y)<cut)) then 
        F(n,m)=p(m)/beta*(x-Delta)/x/(1-exp(y-zx))
      else if ((abs(zx) < cut) .and. (zx-y > cut)) then 
        F(n,m)=p(m)/beta*(x-Delta)/x*(1-exp(-zx))
      else if ((abs(zx)<cut) .and. (-(zx-y) > cut)) then
        F(n,m)= p(n)/beta*(x-Delta)/x*(1-exp(zx))
      else if ((-zx > cut) .and. (abs(zx-y)< cut)) then
        F(n,m)= p(n)/beta*(x-Delta)/x/(1-exp(zx-y))
      else if ((-zx > cut) .and. (-(zx-y)>cut)) then
        F(n,m)= p(n)/beta* (x-Delta)/x/(1-exp(zx-y))
      else if  ((-zx > cut) .and. (-(zx-y)>cut)) then
        F(n,m)= p(n)*(beta)*(x-Delta)/x 
      else if ((abs(zx)<cut) .and. (abs(zx-y)<cut)) then
        s1=(exp(zx/2) - exp(-zx/2))
        s2=(exp((zx-y)/2) - exp(-(zx-y)/2))
        F(n,m)= s0*(x-Delta)/x*s1/s2
      else
        F(n,m)=0
      end if
    end if
if (F(n,m)<0) then
write(*,*) n,m,F(n,m)
end if
  end do
end do
!memory function matrix
do nn=1,Ms
  do mm=1,Ms
    n1=v1(nn)
    n2=v2(nn)
    m1=v1(mm)
    m2=v2(mm)
    cs=0.
    if (n1 == m1) then
      do t=1,Ns
        mat=jjz(m2,t)*jjz(t,n2)+0.5*jjp(m2,t)*jjm(t,n2)+0.5*jjm(m2,t)*jjp(t,n2)
        cs=cs+mat*F(n1,t)
      end do
    end if
    if (n2 == m2) then
      do t=1,Ns
        mat=jjz(n1,t)*jjz(t,m1)+0.5*jjp(n1,t)*jjm(t,m1)+0.5*jjm(n1,t)*jjp(t,m1)
        cs=cs+mat*F(t,n2)
      end do
    end if
    mat=jjz(n1,m1)*jjz(m2,n2)+0.5*jjp(n1,m1)*jjm(m2,n2)+0.5*jjm(n1,m1)*jjp(m2,n2)
    cs=cs-mat*(F(n1,m2)+F(m1,n2))
    Mem(nn,mm)=-iunit*gam*cs
  end do
end do
!  Omega-Matrix
do nn=1,Ms
do mm=1,Ms
cs=0.
if (nn == mm) then
n1=v1(nn)
n2=v2(nn)
cs=cs+(x-En(n1)+En(n2))*Pp(nn)
end if
cs=cs-Mem(nn,mm)
OM(nn,mm)=cs
end do
end do
! Inversion and final result
call matinv(OM,OMinv,Ms,Mp)

do nn=1,Ms
do mm=1,Ms
Phi(nn,mm)=Pp(nn)*OMinv(nn,mm)*Pp(mm)
end do
end do

do k=1,3
do l=1,3
   cs=0
   css=0
   do n=1,Ms
   do m=1,Ms
     n1=v1(n)
     n2=v2(n)
     m1=v1(m)
     m2=v2(m)
     cs = cs+jjj(k,n1,n2)*Phi(n,m)*jjj(l,m2,m1)
     if( n==m )then
        css=css + jjj(k,n1,n2)*Pp(n)*jjj(l,m2,m1)
     end if  
    end do
    end do 
    chi(k,l)=beta*(css-x*cs)
!     chi(k,l)=beta*(-x*cs)
      
    if (abs(x) <= 0.0001) then
       sf(k,l)=-aimag(cs)
     else  
       sf(k,l)=-aimag(cs)*beta*x /(1-exp(-beta*x))                    
    end if
end do
end do



return
end subroutine SuscepComponents   
!--------------------------------------------


subroutine OutputResults
use CommonData
use MatrixElements
implicit none
integer :: i,k,l,n,nn

integer :: type, status_read,m
real :: deltax,x,y,sum, s,st, sff(3,3), ss(3,3)
real :: dksv(3),dkv(3),dkapv(3),xksv(3),xkv(3),xkapv(3)
real :: k0,kn,ks,ksn,kapn,kq,ksq,dcs
! the following quantities used here are defined in CommonData
! outfilename (name of file with results), cefname (name of file 
! with cef-data, strfilename (name of original file with formfactor)
! 
character(len=20) :: f, erg(6)

interface
subroutine SuscepComponents(z)
real :: z
end subroutine
end interface
!write(*,*) 'Hier angekommen'
write(*,'(A)') 'results written into Results/chibfk.res'
open(21,file='./Results/chibfk.res',action='write')

write(21,*) '# results for the dynamical susceptibility of RE ions'
write(21,*) '# with bfk theory  using level-scheme ', cefname
write(21,'(A15,F6.2,A18,F6.2)') '# temperature ', temp, ' coupling constant ', g
write(21,'(A33,3F10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
write(21,*) '# energy loss, real and imaginary part of the dynamical susceptibility'
write(21,*) '# chi(1,1)    chi(1,2)   chi(1,3)' 
write(21,*) '# chi(2,1)    chi(2,2)   chi(2,3)'
write(21,*) '# chi(3,1)    chi(3,2 )  chi(3,3)'  
   

deltax=epsilon
Npoints=(emax-emin)/epsilon
  do n=0,Npoints
    x=emin+n*deltax
    if (abs(x) < 0.00001) then
      x=x+0.00001 
     end if
    call SuscepComponents(x)
    write(21,*) x
    do i=1,3
      write(21,'(5X,6F12.6)') (chi(i,k),k=1,3)
    end do
  end do

end subroutine OutputResults
!--------------------------------------------


program bft
 
call DataRead
call MatrixElementCalculation
call OutputResults

end program
!-------------------------------------------------------------------

