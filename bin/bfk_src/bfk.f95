
!  Program bfk(2707).f95 calculates dynamical  susceptibility and 
!  q-dependent inelastic neutron scattering cross-section for single  
!  RE ions.
!  Input-files: a) file with eigenvalues and eigenstates of the RE ions in a 
!  crystalline electric field, b) file with formfactor, c) parameter file 
!  containing coupling constants, energy range etc. d) table of q, omega values 
!  for which the calculation shall be performed. 

Module CommonData
! This module contains data used by all other subroutines
implicit none
save 
integer :: Ns,Ms !number of states, number of dynamical variables
integer, parameter :: Np=20 !max number of states
integer, parameter :: Mp=100 !max number of transitions
integer, parameter :: Nc=250 !max number of scatt. vect. 
integer :: w1,w2,w3,w4,w5,w6,w7,w8,w9,w10
real :: gl !  Lande factor
real :: diff(3,3), sf(3,3), kv(3), ksv(3), kapv(3), kapnv(3)
complex :: chi(3,3)
real :: E, E0, Es, k1,k2,k3,ks1,ks2,ks3,kap1,kap2,kap3,kapv1,kapv2,kapv3
real :: enloss, kappa 
real :: jav(3)
complex :: Phi(Mp,Mp), Om(Mp,Mp), Ominv(Mp,Mp) !Relaxation function matrix
real :: Mem(Mp,Mp)
real :: emin,emax, k11,k12,k13,k21,k22,k23
real :: a,b,c,qqq(3)
integer ::  Npoints
integer :: mode,mst !scattering mode, type of output
real :: beta,gg, g,ex, gam, temp, cutoff ! 1/k_B T, coupl. const 
real, parameter :: Pi=3.14159, ek2=2.072, meVkT=11.6, r0=-0.54 
complex, parameter :: Iunit=(0.0,1.0)
character(len=20) :: outfilename
character(len=30) :: cefname, parfilename, formfactorname,scatfilename
end module CommonData
!--------------------------------------------------------------------------

module MatrixElements
use CommonData
implicit none
save
integer :: v1(Mp), v2(Mp) !states belonging to one transition
real :: En(Np), p(Np), Pp(Mp) !energy eigen values, Boltzmann factors, static 
!susceptibilities
complex :: Pq(Mp,Mp),Pqinv(Mp,Mp)
real :: Zzz(Np,Np,Np,Np), Yyy(Np,Np,Np,Np)
real ::  jex(3,3),exc(3)
complex :: Ev(Np,Np),jjj(3,Np,Np) !eigenstates, spin matrixelements, magnetisation
complex :: jjx(Np,Np),jjy(Np,Np), jjz(Np,Np), jjp(Np,Np), jjm(Np,Np)
end module MatrixElements
!---------------------------------------------------------------
module FormfactorPreparation
use CommonData
! transformation of table of structure factor kapg(n), strg(n)
! to table kapf(n), strf(n) with equidistant kappa steps
implicit none
save
real :: kapg(Nc),kapf(Nc),strg(Nc),strf(Nc)
integer :: nlg, nlf
contains 
subroutine FormfactorTransformation
implicit none
integer :: m, n, nn
real :: dk, y1,y2,x1,x2,s,x
  dk=0.1;
  nlf=int(kapg(nlg)/dk+0.001)-1
  strf(1)=1.0
  do m=1,nlf
    kapf(m)=dk*(m-1)
  end do
  n=1
  do m=2,nlf
    do while(kapg(n)<=kapf(m)) 
      nn=n
      n=n+1
    end do
    y1=strg(nn)
    y2=strg(nn+1)
    x1=kapg(nn)
    x2=kapg(nn+1)
    x=kapf(m)
    s=y1*(x-x2)/(x1-x2)
    s=s+y2*(x-x1)/(x2-x1)
    strf(m)=s
  end do
end subroutine FormFactorTransformation
real function Strfak(kapa)
! Berechnet Strukturfaktor für beliebige Werte von kappa 
!durch Interpolation
!kapa Betrag des Streuvektors,
!strf(n) Tafel der Werte des Struktur-Faktors
!max kapa-value
implicit none
save
integer :: n,nn
real, intent(in) :: kapa
real :: s1, s2, kappa1,kappa2,dkap
dkap=kapf(2)-kapf(1)
nn=int(kapa/dkap)
if (nn .gt. nlf-2) then
write(*,*) 'Mistake in number of lines in form factor table' 
end if 
kappa1=dkap*nn
kappa2=dkap*(nn+1)
s1=strf(nn+1)
s2=strf(nn+2)
strfak=s1*(kapa-kappa2)/(kappa1-kappa2) + s2*(kapa-kappa1)/(kappa2-kappa1)
end function



end module FormfactorPreparation
!----------------------------------------------------------------------------
subroutine DataRead
!  Reads input-data from command line and input-files  
!  The entries in the command line are:
!  temperature T, coupling constant for the interaction 
!  with conduction electrons,
!  mode of scattering, type of output-file
!  The 5th entry is the name of a file containing information 
!  about the RE ion: number of levels,
!  Eigen-energies En(n) and eigen-states Ev(n,k) of RE ion in 
!  a crystalline electric field.
!  The 6th entry should be the name of a parameter-file containing 
!  information about the energy and momentum range of the calculation,  
!  the names of a file with q and omega values for which the scattering
!  cross section should be calculated, and a file with  the atomic
!  formfactor of the RE ion.  

use CommonData
use MatrixElements
use FormfactorPreparation
implicit none
!save
integer :: i,j, k,kl,n,m, l, ii ,nn, cline, vv1(Mp),vv2(Mp), w, status_read 
integer :: ww, type_of_contence,linetype(50), linenumber, nlines
real :: s,dk,kn,knn,ksn,ks,k0, Evr(Np,Np),Evi(Np,Np), v(Np),f, x, y,gJ,dd
integer ::  io_error,ddi
real ::  x1,x2,x3,y1,y2,y3,strfak
! kv(i) direction of incoming beam, ksv(i) of scattered beam, 
! and kap(i) of vector kappa
complex ::  cs
character(len=30)::  arg 
character(len=30) :: name(7)
character(len=10) :: couplc, tempp, mod
character(len=100) :: line(50), tline(50), bb, zeile
character(len=20) :: redata(14),re(14)
character(len=100) :: number
character(len=100) :: vector
!  read information from command-line
do i=1, iargc() ! number of entries in commandline
    call getarg(i,arg) ! entries in command line
    name(i)= arg
end do
nn=iargc() !number of entries
!check if command line contains information 
if (nn < 6) then
  write(*,*) 'command line emptyor uncomplete.' 
  write(*,*) 'It should contain temperature, electron-coupling constant,'
  write(*,*) 'calculation-mode, output-mode,'
  write(*,*) 'name of file with cef-data, name of parameterfile'
  write(*,*) 'start again'
  stop
end if
w1=0
w2=0
w3=0
w4=0
w5=0
w6=0
w7=0
w8=0
w9=0
w10=0
if (nn>0) then
  tempp=name(1)
  read(tempp,*) x
  temp=x
  w1=1
end if
  if (nn>1) then
  couplc=name(2)
  read(couplc,*) g
  w1=w1+1
end if
if (nn>2) then
  mod=name(3)
  read(mod,*) mode
  w1=w1+1
end if
if (nn>3) then
  mod=name(4)
  read(mod,*) mst
  w1=w1+1
end if
if (nn > 4 ) then
  cefname=name(5)
  w2=1
end if
if (nn > 5) then
  parfilename=name(6)
  w3=1
end if
if (w2==1) then
! Analyse file with cef-data
  open(11,file=cefname,action='read')
  n=1
  do 
    read(11,'(A)',iostat=status_read) line(n)
    if(status_read /= 0) exit
    n=n+1
  end do 
  nlines=n-1
  m=0
  do n=1,nlines
    if (line(n)(1:1)=='#') then
      if ((line(n)(1:4)=='#!d=').or.(line(n)(1:4)=='#!J=')) then
         number=line(n)(5:100)
         read(number,*) dd
         if(line(n)(1:4)=='#!J=') then
          dd=2*dd+1
         end if
          ddi=dd
         m=m+1
         bb=line(n)(5:100)
           do l=1,35
             if (bb(l:l+8) == 'sipffile=') then
                bb=bb(l+9:100)
                do i=1,100
                 if(bb(i:i) == ' ') then
                  bb=bb(1:i)
                 end if
                end do                
                open(12,file=bb,action='read')
                do 
                 read(12,'(A)',iostat=status_read) bb
                 do i=1,100
                  if (bb(i:i+2)=='GJ=') then
                   number=bb(i+3:100)
                   read(number,*) gJ
                  endif
                end do          
                if(status_read /= 0) exit
               end do 
               close (12)
               write(tline(m),'(I5,F12.8)') ddi,gJ 
               write(*,*) 'dimension    gJ'
               write(*,*) tline(m)
             end if
           end do
      end if
      if (line(n)(1:14)=='#! Eigenvalues') then
        m=m+1
        bb=line(n)(15:100)
        k=index(bb,'=')
        tline(m)=bb(k+1:100)
      end if
    else 
      m=m+1
      tline(m)=line(n)
    end if
  end do
  open(12,file='cefworkfile.dat')  
  do n=1,m
    write(12,'(A)') tline(n)
  end do
  close(11)
  close(12)
end if !w2
if (w3==1) then
!Analyse file with parameters
  open(13,file=parfilename,action='read')
  n=1 
  do 
    read(13,'(A)',iostat=status_read) line(n)
    if(status_read /= 0) exit
    n=n+1
  end do
  m=0
  nlines=n-1
  do n=1,nlines
    if (line(n)(1:9)=='#!cutoff=') then
      number=line(n)(10:100)
      read(number,*) cutoff
      w4=1
    end if
    if (line(n)(1:7)=='#!emin=') then
      number=line(n)(8:100)
      read(number,*) emin
      w5=1
    end if
    if (line(n)(1:7)=='#!emax=') then
      number=line(n)(8:100)
      read(number,*) emax
      w5=w5+1
    end if
    if (line(n)(1:10)=='#!Npoints=') then
      number=line(n)(11:100)
      read(number,*) Npoints
      w5=w5+1    
    end if
    if (line(n)(1:4)=='#!E=') then
      number=line(n)(5:100)
      read(number,*) E
      w6=1
    end if
    if (line(n)(1:5)=='#!k1=') then
      vector=line(n)(6:100)
      read(vector,*) k11, k12, k13
      w6=w6+1
    end if
    if (line(n)(1:5)=='#!k2=') then
      vector=line(n)(6:100)
      read(vector,*) k21, k22, k23
      w6=w6+1
    end if
    if (line(n)(1:17)=='#!formfactorname=') then
       formfactorname=line(n)(18:100)
       w7=1
    end if
    if (line(n)(1:15)=='#!scatfilename=') then
       scatfilename=line(n)(16:100)
       w8=1
    end if
  end do
end if ! w3=1
close(13)
!write data into bfkdata-file  
open(14,file='bfkdata.dat',action='write',iostat=io_error)
if (io_error /=0) then 
 write(*,*) 'io_error opening 14'
end if
write(14,*)  w1, w2, w3, w4, w5, w6, w7,w8,w9,w10 
if (w1 == 5) then 
  write(14,*) temp,g,mode, mst
end if
if (w4==1) then
  write(14,*) cutoff
end if
if (w5==3) then
  write(14,*) emax,emin, Npoints
end if
if (w6==3) then
  write(14,*) E
  write(14,*) k11,k12,k13
  write(14,*) k21,k22,k23
end if
if (w7==1) then
  write(14,'(A)') formfactorname
end if  
if( w8==1) then
  write(14,'(A)') scatfilename
end if
    
close(14)
if (w7==1) then !if formfactor data are availablle, they are read-in
  open(17,file=formfactorname,action='read',iostat=io_error)
  if (io_error /= 0) then
     write(*,*) 'io_error while opening file with formfactor'
     stop
  end if
  n=2
  kapg(1)=0.
  strg(1)=1.
  m=0
  einlesen: do
      read(17,*, iostat=status_read) kapg(n), strg(n)
      if (status_read ==0) then
        n=n+1
        m=1
        cycle
      end if
      if (status_read /=0 .and. m/=0) exit
  end do einlesen
  nlg=n-1
  close(17)
    
!transfer of structure factor data to field strworkfile.dat with 
!equidistant steps kap(n)
  call FormfactorTransformation
  open(18,file='strworkfile.dat')
  do n=1,nlf
    write(18,*) kapf(n),strf(n)
  end do
  close(18) 
end if !w7 
!prepare necessary data for the calculations  

beta = 11.6/temp
gam= 2*Pi*g**2
!open different workfiles 
if (w2==1) then
  open(12,file='cefworkfile.dat',action='read')
  read(12,*) Ns, gl 
  read(12,*) (En(n),n=1,Ns)
  do j=1,Ns
    read(12,*) (Evr(j,i), i=1,Ns)
  end do  
  do j=1,Ns
    read(12,*) (Evi(j,i), i=1,Ns)
  end do  
  close(12)
end if
do i=1,Ns
  do j=1,Ns
    Ev(j,i)=cmplx(Evr(j,i),Evi(j,i)) 
  end do
end do
open(14,file='bfkdata.dat')
read(14,*) w1,w2,w3,w4,w5,w6,w7,w8,w9,w10
if (w1==5) then
read(14,*) temp,g,ex,mode, mst
end if
if (w4 ==1) then
  read(14,*) cutoff
end if
if (w5==3) then 
    read(14,*) emax,emin, Npoints
end if
if (w6== 3) then
    read(14,*) E
    read(14,*) k11,k12,k13
    read(14,*) k21,k22,k23
end if
if (w7==1) then
    read(14,'(A)') formfactorname
end if  
if( w8==1) then
    read(14,'(A)') scatfilename
end if
close(14)
if (w7==1) then
  open(18,file='strworkfile.dat')
  n=1
  m=0
  read: do
    read(18,*,iostat=status_read) kapf(n),strf(n)
    if (status_read==0) then
      n=n+1
      m=1   
      cycle
    end if
  if (status_read /= 0 .and. m/=0) exit
  end do read
  nlf=n
close(18)
end if
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
do i = 1,Ns
do j = 1,Ns
cs=(0.,0.)
do k = 1,Ns
cs=cs+conjg(Ev(j,k))*Ev(i,k)
end do
if ((cabs((cs-1)*cs)>= 0.001) )then
write( *,*) 'Mistake in Eigenvektor, accuracy low' 
end if
end do
end do

!the file jmatrix is generated only for checking the correctnes of  the 
!matrix elements
open(19,file='jmatrix.dat')

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
!  if (n1==n2) then
!    vv1(nn)=0
!    vv2(nn)=0
!  end if 
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

subroutine Relmatrix(x)  
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
    s=0.
    if (n1 == m1) then
      do t=1,Ns
        mat=jjz(m2,t)*jjz(t,n2)+0.5*jjp(m2,t)*jjm(t,n2)+0.5*jjm(m2,t)*jjp(t,n2)
        s=s+mat*F(n1,t)
      end do
    end if
    if (n2 == m2) then
      do t=1,Ns
        mat=jjz(n1,t)*jjz(t,m1)+0.5*jjp(n1,t)*jjm(t,m1)+0.5*jjm(n1,t)*jjp(t,m1)
        s=s+mat*F(t,n2)
      end do
    end if
    mat=jjz(n1,m1)*jjz(m2,n2)+0.5*jjp(n1,m1)*jjm(m2,n2)+0.5*jjm(n1,m1)*jjp(m2,n2)
    s=s-mat*(F(n1,m2)+F(m1,n2))
    Mem(nn,mm)=gam*s
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
cs=cs+iunit*Mem(nn,mm)
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
end subroutine Relmatrix

!--------------------------------------------------------
subroutine SuscepComponents(x)
use CommonData
use MatrixElements
implicit none
integer :: k,l,n,m,n1,n2,m1,m2
complex :: cs,css
real :: x
interface
  subroutine Relmatrix(x)
  real :: x
  end subroutine
end interface
call Relmatrix(x)
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
    if (abs(x) <= 0.0001) then
       sf(k,l)=-aimag(cs)
     else  
       sf(k,l)=-aimag(cs)*beta*x /(1-exp(-beta*x))                    
    end if
end do
end do

return
end subroutine SuscepComponents

!---------------------------------------------- 





!-------------------------------------------------
real function ScatFunction(qv,omega)
! calculates the essential part of the scattering cross section  
! as function of the scattering wave vector kappav and energy loss omega
use Commondata
use MatrixElements
implicit none
save
integer :: n,m,k,l,n1,n2,m1,m2
real :: omega, s,q, qv(3),qnv(3)
complex :: cs
interface
subroutine SuscepComponents(z)
real :: z
end subroutine !calculates scattering components J(i,n) Phi(n,m)J(k,m)
end interface
call SuscepComponents(omega)
! calculates different spin components of the dynamical suscept. sf(i,k)
q= sqrt(qv(1)**2+qv(2)**2 + qv(3)**2)
do k=1,3
  qnv(k)=qv(k)/(q + 0.000001)
end do
do k=1,3
  do l=1,3
     diff(k,l)=-qnv(l)*qnv(k)
     if (l==k) then
       diff(k,l)=diff(k,l)+1
     end if
  end do
end do
cs=0
do k=1,3
do l=1,3
  cs=cs+diff(k,l)*sf(k,l)
end do
end do
Scatfunction=cs
return
end function    
!--------------------------------------------

subroutine OutputResults
use CommonData
use MatrixElements
use FormfactorPreparation
implicit none
integer :: i,k,l,n,nn
integer :: type,m, status_read
real :: deltax,x,y,sum, s,st, sff(3,3), ss(3,3),qq(3)
real :: dksv(3),dkv(3),dkapv(3),xksv(3),xkv(3),xkapv(3)
real :: k0,kn,ks,ksn,kapn,kq,ksq,dcs,kmin,kmax,q1,q2,q3
real :: k0k0, ksks,kv0(3), kvs(3), k0v(3)
! the following quantities used here are defined in CommonData
! outfilename (name of file with results), cefname (name of file 
! with cef-data, strfilename (name of original file with formfactor)
! 
character(len=20) :: f(7),erg(7)

interface
real function ScatFunction(kapv,omega)
real :: kapv(3), omega
end function
subroutine SuscepComponents(z)
real :: z
end subroutine
end interface



if (mode==0) then
if (mst == 1) then  
open(21,file='./results/bfk0.res',action='write')
open(22,file='./results/plot0.res',action='write')
end if
if (mst == 2) then
open(21,file='./results/bfk0.res',action='write',position='append')
end if

write(*,*) 'Results are written into ./results/bfk0.res and plot0.res'
write(21,*) '#results of program bfk using level-scheme ', cefname
write(21, *)'# calculation mode', mode
write(21,*) '# results for the dynamical susceptibility of RE ions'
write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
write(21,'(A33,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
write(21,'(A14)') '# energy loss'
write(21,'(A65)') '# real and imaginary part of the dynamical susceptibility'
write(21,'(A65)') '# chi(1,1)                chi(1,2)               chi(1,3)' 
write(21,'(A65)') '# chi(2,1)                chi(2,2)               chi(2,3)'
write(21,'(A65)') '# chi(3,1)                chi(3,2 )              chi(3,3)'  

deltax =(emax-emin)/Npoints
  do n=0,Npoints
    x=emin+n*deltax
    if (abs(x) < 0.00001) then
      x=x+0.00001 
     end if
    call SuscepComponents(x)
    write(21,*) x
    do i=1,3
      write(21,'(5X,6F12.6)') (chi(i,k),k=1,3)
      write(22,'(6F10.3)') x, real(chi(1,1)),aimag(chi(1,1)),real(chi(2,2)), aimag(chi(2,2))
!      write(*,'(6F10.3)') x, real(chi(1,1)),aimag(chi(1,1)),real(chi(2,2)), aimag(chi(2,2))
    end do
  end do
close(21)
close(22)
end if


!-----------------------------------------------------------
if (mode==1) then
if (mst == 1) then  
open(21,file='./results/bfk1.res',action='write')
end if
if (mst == 2) then
open(21,file='./results/bfk1.res',action='write',position='append')
end if
open(22,file='./results/plot1.res',action='write')

write(*,*) 'Results are written into ./results/bfk1.res and plot1/res'
write(21,*) '#results of program bfk using level-scheme ', cefname
write(21, *)'# calculation mode', mode
write(21,*) '# results for the dynamical susceptibility of RE ions'
write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
write(21,'(A33,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
write(21,'(A12,3F6.2)') ' # q-vector ', qqq(1), qqq(2), qqq(3)
write(21,*) '#Calculates S(i,i)=Im Chi(i,i)*coth(beta*omega/2)'
write(21,*) '#and checks sum-rule (S(1,1)+S(2,2)+S(3,3))/Pi = J*(J+1)' 
write(21,*) '#energy loss,    S(1,1)         S(2,2)       S(3,3)'
 
  deltax=(Emax-emin)/Npoints
  s=0
  m=4
  write(22,'(6F10.3)') x, real(chi(1,1)),aimag(chi(1,1)),real(chi(2,2)), aimag(chi(2,2))
!  write(*,'(6E10.3)') x, real(chi(1,1)),aimag(chi(1,1)),real(chi(2,2)), aimag(chi(2,2))

  do n=1,10

    y=n/10.
    x=y**m*deltax 
    call SuscepComponents(x)
    do i=1,3
      do k=1,3
        ss(i,k)=sf(i,k)*(1 + exp(-x*beta))
      end do
    end do
    write(21,'(4G14.5)') x,(ss(i,i),i=1,3)
    write(22,'(6F10.3)') x, ss(1,1), ss(2,2), ss(3,3)

!    write(*,'(6F10.3)') x, ss(1,1),ss(2,2),ss(3,3)
  end do
  do n=1,10

    s=s+(ss(1,1)+ss(2,2)+ss(3,3))*m*x/(y*10.)
  end do

  do n=2,Npoints
    x=(n+0.00001)*deltax
    call SuscepComponents(x)
    do i=1,3
      do k=1,3
        ss(i,k)=sf(i,k)*(1 + exp(-x*beta))
      end do
    end do
    write(21,'(4G14.5)') x,(ss(i,i),i=1,3)
    write(22,'(6F10.3)') x, ss(1,1), ss(2,2), ss(3,3)
!    write(*,'(6F10.3)') x, ss(1,1),ss(2,2),ss(3,3)

    s=s+(ss(1,1)+ss(2,2)+ss(3,3))*deltax
  end do
  sum=s/Pi
  write(21,*)'#Integral (ss(1,1)+ss(2,2)+ss(3,3)/Pi: ', sum, ( Ns**2-1.)/4 
close(21)
close(22)
end if 

!--------------------------------------------------------------------------
if (mode==2) then
if (mst == 1) then  
open(21,file='./results/bfk2.res',action='write')
end if
if (mst == 2) then
open(21,file='./results/bfk2.res',action='write',position='append')
end if
open(22, file='./results/plot2.res',action='write')
write(*,*) 'Results are written into ./results/bfk2.res and plot2.res'

write(21,*) '#results of program bfk using level-scheme', cefname
write(21,*)'# calculation mode', mode
write(21,*) '# Calculation of scattering cross section for given set of '
write(21,*) '# scattering wave vector q and energy loss'
write(21,*) '# provided by input-file ', scatfilename
write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
write(21,'(A33,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
write(21,*) '# formfactor   ', formfactorname 
 if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if
write(21,*) '# q (A^-1) energy loss (meV) cross section (1/meV,steradian)'
    open(23,file=scatfilename, action='read')

    calculation: do
      read(23,*,iostat=status_read) kapv(1), kapv(2), kapv(3), enloss
      if (status_read /= 0 ) exit
    kappa=sqrt(kapv(1)**2+kapv(2)**2+kapv(3)**2)
      do i=1,3
        kapnv(i)=kapv(i)/kappa
      end do
      if (w7==1) then
        st=strfak(kappa)
      else
        st=1
      end if
      dcs=(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)

      write(21,'(3G14.5)') kappa,enloss,dcs
!      write(*,'(3G14.5)') kappa,enloss,dcs
      write(22,'(3G14.5)') kappa,enloss,dcs
      
    end do calculation
    close(22)
    close(21)
    close(23)
  end if      



!-------------------------------------------
if (mode==3) then
if (mst == 1) then  
open(21,file='./results/bfk3.res',action='write')
end if
if (mst == 2) then
open(21,file='./results/bfk3.res',action='write',position='append')
end if
open(22,file='./results/plot3.res',action='write')
write(*,*) 'Results are written into ./results/bfk3.res and ./results/plot3.res'
write(21,*) '#results of program bfk using level-scheme ', cefname
write(21,*)'# calculation mode', mode
write(21,*) '# calculation of different components of the scattering functions'
write(21,*) '# as function of q and omega provided by input-file ',scatfilename 
write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
write(21,'(A33,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
 if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if
write(21,*) '#  Q,  energy-loss '
write(21,*) '#                      S(1,1),S(1,2),S(1,3)'
write(21,*) '#                      S(2,1),S(2,2),S(2,3)'
write(21,*) '#                      S(3,1),S(3,2),S(3,3)'

open(23,file=scatfilename,action='read')
status_read=0
calc: do 
    read(23,*,iostat=status_read) kapv(1),kapv(2),kapv(3),enloss
    if (status_read /= 0 ) exit
    kappa=sqrt(kapv(1)**2+kapv(2)**2+kapv(3)**2)
    do i=1,3
      kapnv(i)=kapv(i)/kappa
    end do
    if (w7==1) then
      st=strfak(kappa)
    else
      st=1
    end if  
    call SuscepComponents(enloss)
    write(21,'(3G14.5)') kappa,enloss
    do k=1,3
    do l=1,3
     sff(k,l)= (r0*gl/2*st)**2/Pi*sf(k,l)
     end do
     end do
     do k=1,3
      write(21,'(5X,3G14.5)') (Sff(k,l), l=1,3)
     end do  
     write(22,'(11G12.3)') kappa, enloss, Sff(1,1),Sff(1,2),Sff(1,3),Sff(2,1),&
     Sff(2,2),Sff(2,3),Sff(3,1),Sff(3,2),Sff(3,3)
    
!write(*,*) kappa, enloss, Sff(1,1)
end do calc
close(22)
close(21)
close(23)
end if      

!-----------------------------------------------
if (mode==4) then

  if (mst == 1) then  
    open(21,file='./results/bfk4.res',action='write')
  end if
  if (mst == 2) then
    open(21,file='./results/bfk4.res',action='write',position='append')
  end if
write(*,*) 'Results are written into ./results/bfk4.res'

  write(21,'(A43,A10)')'#results of program bfk using level-scheme', cefname
  write(21,'(A32,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
  write(21,*) '# calculation mode m=  ', mode
  write(21,'(A)') ' # calculation of the cross section for fixed direction k and'
  write(21,'(A)') ' # energy of the  incident beam, fixed direction but variable' 
  write(21,'(A)') ' # length of the scattering wave vector q '
  write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
 if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if
  E0=E
  k1=k11
  k2=k12
  k3=k13
  kap1=k21
  kap2=k22
  kap3=k23
  write(21,'(A)') ' # E0,    k1,   k2,   k3;   q1,   q2,   q3 '
  write(21,'(7F6.2)') E,k1,k2,k3,kap1,kap2,kap3
  write(21,*) '#q (in A^(-1)), energy-loss (meV), cross section(barns/meV,stradian)'

  k0=sqrt(E0/ek2)
  kn=sqrt(k1**2 + k2**2 +k3**2)
  kv(1)=k1*k0/kn
  kv(2)=k2*k0/kn
  kv(3)=k3*k0/kn
  kapn=sqrt(kap1**2 + kap2**2 + kap3**2)
  kapv(1)=kap1*k0/kapn 
  kapv(2)=kap2*k0/kapn
  kapv(3)=kap3*k0/kapn
  do i=1,3
    dkapv(i)=kapv(i)/Npoints
  end  do
  do n=-Npoints,Npoints
     do i=1,3
       xkapv(i)=n *dkapv(i) +0.0001
       ksv(i)= kv(i)-xkapv(i)
     end do
     kappa=sqrt(xkapv(1)**2 + xkapv(2)**2 + xkapv(3)**2)
     do i=1,3     
       kapnv(i)=xkapv(i)/kappa
     end do     
     ksq=ksv(1)*ksv(1)+ksv(2)*ksv(2)+ksv(3)*ksv(3)
     ks=sqrt(ksq)
     if (w7==1) then
        st=strfak(kappa)
     else
       st=1
    end if
    enloss=ek2*(k0**2-ksq)
     dcs=ScatFunction(xkapv,enloss)
     dcs=ks/k0*(r0*gl/2*st)**2/Pi*dcs
      write(21,'(3G14.5)') kappa,enloss,dcs
!      write(*,'(3G14.5)') kappa,enloss,dcs 
!      write(*,*) qqq(1),qqq(2),qqq(3)  
   end do
end if

!----------------
if (mode==5) then
if (mst == 1) then  
open(21,file='./results/bfk5.res',action='write')
end if
if (mst == 2) then
open(21,file='./results/bfk5.res',action='write',position='append')
end if

write(*,*) 'Results are written into ./results/bfk5.res'

  write(21,'(A43,A10)') '#results of program bfk using level-scheme ',cefname
  write(21,*) '# calculation mode m=  ', mode
  write(21,'(A)') ' # calculation of cross section for fixed direction k and energy E'
  write(21,'(A)') ' # of incident beam, fixed direction ks but variable energy of '
  write(21,'(A)') ' # scattered beam. vector q = vector k0 - vector ks '
  write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
  write(21,'(A29,F6.2)') ' # spin exchange coupling ex= ', ex
  write(21,'(A33,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
 if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if
    E0=E
  k1=k11
  k2=k12
  k3=k13
  ks1=k21
  ks2=k22
  ks3=k23
  write(21,'(A)') ' # E,    k1,   k2,   k3;   ks1,  ks2,  ks3 '
  write(21,'(7F6.2)') E0,k1,k2,k3,ks1,ks2,ks3
  write(21,*) '#q (in A^(-1)), energy-loss (meV), cross section(barns/meV,stradian)'

  k0=sqrt(E0/ek2)
  kn=sqrt(k1**2+k2**2+k3**2)
! normalised wave vector of incoming particles
  kv(1)=k1*k0/kn
  kv(2)=k2*k0/kn
  kv(3)=k3*k0/kn
! normalised wave vector of scattered particles  
  ksn=sqrt(ks1**2+ks2**2+ks3**2)
  ksv(1)=ks1*k0/ksn
  ksv(2)=ks2*k0/ksn
  ksv(3)=ks3*k0/ksn
  do i=1,3
    dksv(i)=ksv(i)/Npoints
  end  do
    do n=-Npoints,Npoints
!      if (n/=0) then 
      do i=1,3
        xksv(i)=float(n)*dksv(i)
        xkapv(i)=-xksv(i)+kv(i) ! scattering vector

      end do


      ksq=xksv(1)**2+xksv(2)**2+xksv(3)**2
      ks=sqrt(ksq)
! length of scattering vector:
      kappa=sqrt(xkapv(1)**2+xkapv(2)**2 + xkapv(3)**2+0.00001)
      do i=1,3
         kapnv(i)=xkapv(i)/kappa ! normalised length of scattering wave vector
      end do
      if (w7==1) then
        st=strfak(kappa)
      else
        st=1
      end if
! energy loss
      enloss=ek2*(k0**2-ks**2)
!differential crosssection 
      dcs=ScatFunction(xkapv,enloss)
!      write(*,'(3G14.5)') kappa,enloss,dcs
      write(21,'(3G14.5)') kappa,enloss,dcs
      dcs=ks/k0*(r0*gl/2*st)**2/Pi*dcs
!   end if  
  end do
  
end if

!-----------------------------------------------
if (mode==6) then
  if (mst == 1) then  
    open(21,file='./results/bfk6.res',action='write')
  end if
  if (mst == 2) then
    open(21,file='./results/bfk6.res',action='write',position='append')
  end if
write(*,*) 'Results are written into ./results/bfk6.res'

  write(21,'(A43,A10)') '#results of program bfk using level-scheme ', cefname
  write(21,*) '# calculation  mode m=  ', mode
  write(21,'(A)') ' # calculation of cross section for fixed direction ks and '
  write(21,'(A)') ' # energy Es of scattered beam, fixed direction k but variable' 
  write(21,'(A)') ' # energy of incident beam'
  write(21,'(A17,F6.2,A42,F6.2)') ' # temperature T=', temp, 'K, coupling &
with conduction electrons g= ', g 
  write(21,'(A33,3G10.3)') ' # magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
 if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if

  Es=E
  k1=k11
  k2=k12
  k3=k13
  ks1=k21
  ks2=k22
  ks3=k23

  write(21,'(A)') ' # Es,     k1,    k2,   k3;   ks1,  ks2,  ks3 '
  write(21,'(7F6.2)') Es,k1,k2,k3,ks1,ks2,ks3
  write(21,*) '#q (in A^(-1)), energy-loss (meV), cross section(barns/meV,stradian)'


  ks=sqrt(Es/ek2)
  kn=sqrt(k1**2+k2**2+k3**2)
!normalised wave vector of incoming particles 
  kv(1)=k1*ks/kn
  kv(2)=k2*ks/kn
  kv(3)=k3*ks/kn
  ksn=sqrt(ks1**2+ks2**2+ks3**2)
! wave vector scattered particle  
  ksv(1)=ks1*ks/ksn
  ksv(2)=ks2*ks/ksn
  ksv(3)=ks3*ks/ksn

  do i=1,3
    dkv(i)=kv(i)/Npoints
  end  do

  do n=-Npoints,Npoints
      do i=1,3
        xkv(i)=kv(i)+n*dkv(i)
        kapv(i)=-ksv(i)+xkv(i)
      end do
      kq=xkv(1)*xkv(1)+xkv(2)*xkv(2)+xkv(3)*xkv(3)
      k0=sqrt(kq)
      kappa=sqrt(kapv(1)*kapv(1)+kapv(2)*kapv(2) + kapv(3)*kapv(3))
      do i=1,3
        kapnv(i)=kapv(i)/kappa
      end do
      if (w7==1) then
        st=strfak(kappa)
      else
       st=1
      end if
      do k=1,3
        do l=1,3
          diff(k,l)=-kapnv(l)*kapnv(k)
          if (l==k) then
             diff(k,l)=diff(k,l)+1
          end if
        end do
      end do
      ksq=ks**2
      kq=k0**2
      enloss=ek2*(kq-ksq)
!differential crosssection
      dcs=ks/k0*(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)
      write(21,'(3G14.5)') kappa,enloss,dcs
!      write(*,'(3G14.5)') kappa,enloss,dcs
  end do
end if


if (mode==6) then
  if (mst == 1) then  
    open(21,file='./results/bfk6.res',action='write')
  end if
  if (mst == 2) then
    open(21,file='./results/bfk6.res',action='write',position='append')
  end if
write(*,*) 'Results are written into ./results/bfk6.res'

  write(21,'(A43,A10)') '#results of program bfk using level-scheme ', cefname
  write(21,'(A12,F6.2,A18,F6.2)') '#temperature ', temp, ' coupling constant ', g
  write(21,'(A30,3G10.3)') '#magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
  write(21,*) '#calculation  mode m=  ', mode
 if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if
  write(21,'(A)') '#calculation of cross section for fixed direction ks and energy Es'
  write(21,'(A)') '#of scattered beam, fixed direction k but variable energy of '
  write(21,'(A)') '#incident beam'

  Es=E
  k1=k11
  k2=k12
  k3=k13
  ks1=k21
  ks2=k22
  ks3=k23

  write(21,'(A)') '# Es,     k1,    k2,   k3;   ks1,  ks2,  ks3 '
  write(21,'(7F6.2)') Es,k1,k2,k3,ks1,ks2,ks3
  write(21,*) '#q (in A^(-1)), energy-loss (meV), cross section(barns/meV,stradian)'


  ks=sqrt(Es/ek2)
  kn=sqrt(k1**2+k2**2+k3**2)
!normalised wave vector of incoming particles 
  kv(1)=k1*ks/kn
  kv(2)=k2*ks/kn
  kv(3)=k3*ks/kn
  ksn=sqrt(ks1**2+ks2**2+ks3**2)
! wave vector scattered particle  
  ksv(1)=ks1*ks/ksn
  ksv(2)=ks2*ks/ksn
  ksv(3)=ks3*ks/ksn
  
  do i=1,3
    dkv(i)=kv(i)/Npoints
  end  do

  do n=-Npoints,Npoints
!    if (n /= 0 ) then
      do i=1,3
        xkv(i)=kv(i)+n*dkv(i)
        kapv(i)=-ksv(i)+xkv(i)
      end do
      kq=xkv(1)*xkv(1)+xkv(2)*xkv(2)+xkv(3)*xkv(3)
      k0=sqrt(kq)
      kappa=sqrt(kapv(1)*kapv(1)+kapv(2)*kapv(2) + kapv(3)*kapv(3))
      do i=1,3
        kapnv(i)=kapv(i)/kappa
      end do
      if (w7==1) then
        st=strfak(kappa)
      else
        st=1
      end if      
      do k=1,3
        do l=1,3
          diff(k,l)=-kapnv(l)*kapnv(k)
          if (l==k) then
             diff(k,l)=diff(k,l)+1
          end if
        end do
      end do
      ksq=ks**2
      kq=k0**2
      enloss=ek2*(kq-ksq)
!differential crosssection
      dcs=ks/k0*(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)
      write(21,'(3G14.5)') kappa,enloss,dcs
 !   end if
  end do
end if

!-------------------------------------------------------------------
if (mode==7) then
  if (mst == 1) then  
    open(21,file='./results/bfk7.res',action='write')
  end if
  if (mst == 2) then
    open(21,file='./results/bfk7.res',action='write',position='append')
  end if
  write(*,*) 'Results are written into ./results/bfk7.res'
  write(21,'(A43,A10)') '#results of program bfk using level-scheme ', cefname
  write(21,'(A12,F6.2,A18,F6.2)') '#temperature ', temp, ' coupling constant ', g
  write(21,'(A30,3G10.3)') '#magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
  write(21,*) '#calculation  mode m=  ', mode
  if (w7==0) then 
  write(21,*) '# formfactor F(Q)=1'
  end if 
 write(21,'(A)') '#calculation of cross section for fixed scattering vector q' 
  write(21,'(A)') '#fixed direction k but variable energy of '
  write(21,'(A)') '#incident beam'

!direction of incoming particle 
  k1=k11
  k2=k12
  k3=k13
!scattering wavevector
  q1=k21
  q2=k22
  q3=k23
  qq(1)=q1
  qq(2)=q2
  qq(3)=q3
  write(21,'(A)') '#  q1,    q2,   q3;   k1,  k2,  k3 '

  write(21,'(7F6.2)') q1,q2,q3,k11,k12,k13
  write(*,'(7F6.2)') q1,q2,q3,k11,k12,k13
  write(21,*) '#q (in A^(-1)), energy-loss (meV), cross section(barns/meV,stradian)'

  kmin = Emin/ek2
  kmax = Emax/ek2


  kn=sqrt(k1**2+k2**2+k3**2)
!normalised wave vector of incoming particles 
  kv0(1)=k1*kmin/kn
  kv0(2)=k2*kmin/kn
  kv0(3)=k3*kmin/kn
  dkv(1)=k1*(kmax-kmin)/kn/Npoints
  dkv(2)=k2*(kmax-kmin)/kn/Npoints
  dkv(3)=k3*(kmax-kmin)/kn/Npoints


  

  do n=1,Npoints
      do i=1,3
        kv(i)=kv0(i)+n*dkv(i)
        ksv(i)=kv(i)+ qq(i)
      end do
      k0k0=kv(1)*kv(1) +kv(2)*kv(2) +kv(3)*kv(3) 
      ksks=ksv(1)*ksv(1) +ksv(2)*ksv(2) +ksv(3)*ksv(3) 
      k0=sqrt(k0k0)
      ks=sqrt(ksks) 
      kappa=sqrt(qq(1)*qq(1)+qq(2)*qq(2) +qq(3)*qq(3))
      do i=1,3
        kapnv(i)=kapv(i)/kappa
       end do
       if (w7==1) then
         st=strfak(kappa)
       else
         st=1
      end if
      do k=1,3
        do l=1,3
          diff(k,l)=-kapnv(l)*kapnv(k)
          if (l==k) then
             diff(k,l)=diff(k,l)+1
          end if
        end do
      end do
      enloss=ek2*(k0k0-ksks)
!differential crosssection
      dcs=ks/k0*(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)
      write(21,'(3G14.5)') kappa,enloss,dcs
      write(*,'(3G14.5)') kappa,enloss,dcs
 !   end if
  end do
end if
close(21)
end subroutine OutputResults
!--------------------------------------------

program bft
 
call DataRead
call MatrixElementCalculation
call OutputResults

end program
!-------------------------------------------------------------------
