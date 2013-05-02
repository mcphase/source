!  Program bfk2103.f95 calculates unpolarised neutron scattering from RE ions
!  Uses as input the spin-eigenstates of the ion,  
!  Output is the scattering cross section for different 
!  scattering geometries.
 
module CommonData
! This module contains data and variables used by all other subroutines 
implicit none
save 
integer :: Ns,Ms !number of states, number of dynamical variables
integer, parameter :: Np=20 !max number of states
integer, parameter :: Mp=100 !max number of transitions
integer, parameter :: Nc=250 !max number of scatt. vect. 
integer :: Npoints=200
real :: gl !  Lande factor
real ::  sf(3,3), kv(3), ksv(3), kapv(3), kapnv(3)
real :: E0, Es, k1,k2,k3,ks1,ks2,ks3,kap1,kap2,kap3,kapv1,kapv2,kapv3
real :: enloss, kappa, Emax
real :: jav(3) !calculated magnetisation 
integer :: mode,mst !scattering mode, output-file type
real :: beta, temp !temperature (K), rec. temperatur (meV)
real :: g, gam ! dimensionless coupl. const JexN(0), 2Pi g^2 
real, parameter :: Pi=3.14159, ek2=2.072, meVkT=11.6, cutoff=100., r0=-0.54 
complex, parameter :: Iunit=(0.0,1.0)
character(len=20) :: outfilename
character(len=20) :: cefname, strfilename, tablename !names of input files
end module CommonData
!--------------------------------------------------------------------------

module MatrixElements
use CommonData
implicit none
save
integer :: v1(Mp), v2(Mp) !states belonging to one transition
real :: En(Np), p(Np), Pp(Mp) !energy eigen values, Boltzmann factors, static 
! susceptibilities
complex :: Ev(Np,Np),jjj(3,Np,Np) !eigenstates, spin matrixelements 
complex :: jjx(Np,Np),jjy(Np,Np), jjz(Np,Np), jjp(Np,Np), jjm(Np,Np)
end module MatrixElements
!---------------------------------------------------------------

module Ffp
use CommonData
! transformation of table of structure factor kapg(n), strg(n)
! to table kapf(n), strf(n) with equidistant kappa steps
implicit none
save
real :: kapg(Nc),kapf(Nc),strg(Nc),strf(Nc)
integer :: nlg, nlf
contains 
subroutine Fft
!transforms any formfactor table into table with equidistant
!wave vector steps 
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
end subroutine Fft
end module Ffp

!----------------------------------------------------------------------------
subroutine DataRead
!  Reads input-data. 
!  Eigen-energies En(n) and eigen-states Ev(n,k) of RE ion in 
!  a crystalline electric field.
!  It asks for the coupling with conduction electrons g=J_ex N(0), 
!  temperature T and name of input-file. In following runs this information can  
!  be provided by command-line input.
!  Controll parameters to run the program are asked for in the first run and 
!  stored in the files bfkdata0.dat, bfkdata12.dat or bfkdata345.dat
 

use CommonData
use MatrixElements
use Ffp
implicit none
save
integer :: i,j, k,l, status_read,n,nn,m,cline, vv1(Mp),vv2(Mp) 
integer :: zz,ww, nlines
real :: s,dk,kn,knn,ksn,ks,k0, Evr(Np,Np),Evi(Np,Np), v(Np),f, x, y 
real ::  x1,x2,x3,y1,y2,y3,strfak
complex ::  cs
character(len=10)::  arg 
character(len=10) :: name(5)
character(len=10) :: couplc, tempp, mod,filet
character(len=150) :: line(150),tline(150), bb, zeile
character(len=20) :: redata(12),re(12)


!  read information from command-line   

k=0
do i=1, iargc() ! number of entries in commandline
    call getarg(i,arg) ! entries in command line
    name(i)= arg
end do
n=iargc() !number of entries
!check if command line contains information 
if (n==0) then 
  k=0
end if
if (n==1) then
  k=1
  couplc=name(1)
  read(couplc,*) x
  g=x
end if
if (n == 2) then
  k=2
  couplc=name(1)
  read(couplc,*) x
  g=x
  tempp=name(2)
  read(tempp,*) y
  temp=y
end if
if (n == 3  ) then
  k=3
  couplc=name(1)
  read(couplc,*) x
  g=x
  tempp=name(2)
  read(tempp,*) y
  temp=y
  mod=name(3)
  read(mod,*) zz
  mode=zz
end if
if (n == 4  ) then
  k=4
  couplc=name(1)
  read(couplc,*) x
  g=x
  tempp=name(2)
  read(tempp,*) y
  temp=y
  mod=name(3)
  read(mod,*) zz
  mode=zz
  filet=name(4)
  read(filet,*) ww
  mst=ww
end if
if (n == 5) then
  k=5 
  write(*,*) 'more than 4 entries in commandline, input ignored'
end if
if (k==0 .or. k==5) then
  cline=0
end if
if (k==1) then 
  write(*,*) 'Input for temperature is missing'
  write(*,*) 'temperature (in K)'
  read(*,*) temp
  cline=1
end if
if (k==2 .or. k==3 .or. k==4) then
  cline=1
end if 
! cline=0 commandline empty, generate new datafile
! cline=1, no new datafile required
if (cline == 0) then
  write(*,*) 'There were no data transferred via the command-line'
  write(*,*) 'You may use the command-line in the next run'
  write(*,*) 'to specify the dimensionless coupling constant g=Jex N(0),'
  write(*,*) 'temperature, mode of calculation m=0,1,2,3,4,5 '
  write(*,*) 'and type of output m=1 (overwrite), m=2 append'
  write(*,*) 'all other quantities are assumed to be unchanged.'   

  write(*,*) 'type dimensionless coupling constant g'
  read(*,*) g
  write(*,*) 'temperature (in K)'
  read(*,*) temp

  write(*,*) 'give name of input file with CEF data (normally levels.cef' 
  write(*,*) 'generated by so1ion'
  read(*,*) cefname 
  outfilename='bfkresults.dat'
  write(*,*) 'give name of file with formfactor (output file formfactor'
  write(*,*) 'of program formfactor)'
  read(*,*) strfilename

  write(*,*) 'Determine  type of calculation and scattering geometry'
  write(*,*) 'mode=0: calculates Im chi(i,i)*coth(beta*omega)'
  write(*,*) 'for 0<omega<Emax at Npoints'
  write(*,*) 'mode=1: calculates scattering function S(Q,omega) for' 
  write(*,*) 'given set of scattering wave vectors (in 1/A)' 
  write(*,*) 'and energy loss omega (in meV) provided by a data file'
  write(*,*) 'mode=2: calculates different components of the scattering'
  write(*,*) 'function for given wave vector and energy loss'
  write(*,*) 'for modes 3-5 different scattering geometries are used'
  write(*,*) 'mode=3: direction and energy of incoming beam fixed,'
  write(*,*) 'scattering wave vector variable ' 
  write(*,*) 'mode=4: direction and energy of incoming beam fixed,'
  write(*,*) 'direction of scattered beam fixed, energy variable' 
  write(*,*) 'mode=5: direction of incoming beam fixed, energy variable'
  write(*,*) 'direction and energy of scattered beam fixed' 
  write(*,*) 'give number of mode m=0,1,2,3,4,5'
  read(*,*) mode
  write(*,*) 'give storage mode of results: mst=1: results of previous run'
  write(*,*) 'are overwritten, mst=2: all results' 
  write(*,*) 'are stored consecutively, new results are  appended.' 
  read(*,*)  mst   

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
  open(15,file=cefname)
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
      if (line(n)(1:4)=='#!J=') then
         m=m+1
         bb=line(n)
         do k=1,12
           do l=1,30
             i=5+l 
             if (bb(i:i+1) == re(k)) then
               tline(m)=redata(k) 
             end if
           end do
         end do
      end if
      if (line(n)(1:15)=='#! Eigenvalues=') then
        m=m+1
        tline(m)=line(n)(16:150)
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


  open(12,file='bfkdata0.dat',action='write')
  open(13,file='bfkdata12.dat',action='write')
  open(14,file='bfkdata345.dat',action='write')
  write(12,'(A)') cefname
  write(12,'(A)') outfilename
  write(13,'(A)') cefname
  write(13,'(A)') outfilename
  write(13,'(A)') strfilename
  write(14,'(A)') cefname
  write(14,'(A)') outfilename
  write(14,'(A)') strfilename

  if ( mode==0) then
    write(*,*) 'give maximum value of energy loss and number of points'
    read(*,*)   Emax, Npoints
    write(12,*) Emax, Npoints
  end if

  if ( mode==1 .or. mode == 2) then
    write(*,*) 'give name of file with table of scattering vectors' 
    write(*,*) 'and values of energy loss'
    read(*,*) tablename
    write(13,'(A)') tablename
  end if  

  if (mode == 3) then
    write(*,*) 'E0 (energy of incoming beam)'
    read(*,*) E0
    write(*,*) 'k1,k2,k3 (direction of incoming beam)'    
    read(*,*) k1,k2,k3
    write(*,*) 'kap1,kap2,kap3 (direction of scattering wave-vector kappa)'
    read(*,*) kap1, kap2,kap3
    write(14,*)  E0
    write(14,*) k1,k2,k3
    write(14,*) kap1,kap2,kap3
  end if

  if (mode == 4) then
    write(*,*) 'E0 (energy of incoming beam in mev)'
    read(*,*) E0
    write(*,*) 'k1,k2,k3 (direction of incoming beam)'
    read(*,*) k1,k2,k3
    write(*,*) 'ks1,ks2,ks3 (direction of scattered beam)'
    read(*,*) ks1,ks2,ks3
    write(14,*) E0
    write(14,*) k1,k2,k3
    write(14,*) ks1,ks2,ks3
  end if

  if (mode == 5) then
    write(*,*) 'Es (energy of scattered beam)'
    read(*,*) Es
    write(*,*) 'k1,k2,k3 (direction of incoming beam)'  
    read(*,*) k1,k2,k3 
    write(*,*) ' ks1, ks2, ks3 (direction  of scattered beam)'
    read(*,*) ks1,ks2,ks3
    write(14,*) Es   
    write(14,*) k1,k2,k3
    write(14,*) ks1,ks2,ks3
  end if

  close(12)
  close(13)
  close(14)

  open(17,file=strfilename,action='read')

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
  nlg=n

  close(17)
!transfer of structure factor data to field strworkfile.dat with 
!equidistant steps kap(n)

!-------------------------------------

  call Fft
  open(18,file='strworkfile.dat')
  do n=1,nlf
     write(18,*) kapf(n),strf(n)
  end do
  close(18) 
end if  !end cline=0

beta = 11.6/temp
gam= 2*Pi*g**2

!open different workfiles 
if (mode == 0) then
  open(12,file='bfkdata0.dat')
  read(12,'(A)') cefname
  read(12,'(A)') outfilename
  read(12,*,iostat=status_read) Emax, Npoints
  if (status_read /= 0) then
    write(*,*) 'data missing in input-file bfkdata0.dat'
    write(*,*) 'give values for Emax and Npoints'
    read (*,*) Emax, Npoints
    write(12,*) Emax, NPoints
  end if  
end if

if((mode == 1) .or. (mode == 2))then
  open(13,file='bfkdata12.dat')
  read(13,'(A)') cefname
  read(13,'(A)') outfilename
  read(13,'(A)') strfilename
  read(13,'(A)',iostat=status_read) tablename  
  if (status_read /= 0) then
    write(*,*) 'data missing in input-file bfkdata12.dat'
    write(*,*) 'give name of file with data for kappa and energy loss'
    read (*,*) tablename
    write(13,'(A)') tablename 
  end if   
end if
  
if((mode==3).or.(mode==4).or.(mode==5)) then 
  open(14,file='bfkdata345.dat')
  read(14,'(A)') cefname
  read(14,'(A)') outfilename
  read(14,'(A)') strfilename
end if   

if (mode == 3)  then
  read(14,*,iostat=status_read) E0
  if (status_read /= 0) then
    write(*,*) 'input data missing'
    write(*,*) 'give energy value and direction of wave vectors'
    write(*,*) 'in the following lines' 
    write(*,*) 'E0 (energy of incoming beam in mev)'
    read(*,*) E0     
    write(*,*) 'k1,k2,k3 (direction of incoming beam)'
    read(*,*)  k1,k2,k3
    write(*,*) 'kap1,kap2,kap3 (direction of scattering vector)'
    read(*,*)  kap1,kap2,kap3
    write(14,*) E0
    write(14,*)  k1,k2,k3
    write(14,*) kap1,kap2,kap3
  else
    read(14,*)  k1,k2,k3
    read(14,*) kap1,kap2,kap3
  end if
end if


if (mode == 4) then
  read(14,*,iostat=status_read) E0
  if (status_read /= 0) then
    write(*,*) 'input data missing'
    write(*,*) 'give energy value and direction of wave vectors'
    write(*,*) 'in the following lines' 
    write(*,*) 'E0 (energy of incoming beam in mev)'
    read(*,*) E0
    write(*,*) 'k1,k2,k3 (direction of incoming beam)'
    read(*,*) k1,k2,k3
    write(*,*) 'ks1,ks2,ks3 (direction of scattered beam)'
    read(*,*) ks1,ks2,ks3
    write(14,*) E0
    write(14,*) k1,k2,k3
    write(14,*) ks1,ks2,ks3
  else
    read(14,*) k1,k2,k3
    read(14,*) ks1,ks2,ks3
  end if  
end if 

if (mode == 5) then   
  read(14,*,iostat=status_read) Es
  if (status_read /= 0) then
    write(*,*) 'input data missing'
    write(*,*) 'give energy value and direction of wave vectors'
    write(*,*) 'in the following lines' 
    write(*,*) 'Es (energy of scattered  beam in mev)'
    read(*,*) Es     
    write(*,*) 'k1,k2,k3 (direction of incoming beam)'
    read(*,*)  k1,k2,k3
    write(*,*) 'ks1,ks2,ks3 (direction of scattered beam)'
    read(*,*)  ks1,ks2,ks3
    write(14,*) Es
    write(14,*)  k1,k2,k3
    write(14,*) ks1,ks2,ks3
  else
    read(14,*)  k1,k2,k3
    read(14,*) ks1,ks2,ks3
  end if
end if

  close(12)
  close(14)
  close(13)
  open(16,file='cefworkfile.dat',action='read')
  read(16,*) Ns, gl 
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
  open(18,file='strworkfile.dat')
  n=1 
  do
    read(18,*,iostat=status_read) kapf(n), strf(n)
    if (status_read /= 0 ) exit
    n=n+1
  end do
  nlf =n

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
write( *,*) 'Mistake in Eigenvektor' 
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
cs=cs-(x-En(n1)+En(n2))*Pp(nn)
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
   do n=1,Ms
   do m=1,Ms
     n1=v1(n)
     n2=v2(n)
     m1=v1(m)
     m2=v2(m)
     cs = cs+jjj(k,n1,n2)*Phi(n,m)*jjj(l,m2,m1)
    end do
    end do     
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

real function Strfak(kapa)
! Berechnet Strukturfaktor für beliebige Werte von kappa 
!durch Interpolation
!kapa Betrag des Streuvektors,
!strf(n) Tafel der Werte des Struktur-Faktors
!max kapa-value
use CommonData
use Ffp
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
!-------------------------------------------------------------------

real function ScatFunction(qv,omega)
! calculates the essential part of the scattering cross section  
! as function of the scattering wave vector kappav and energy loss omega
use Commondata
use MatrixElements
implicit none
save
integer :: n,m,k,l,n1,n2,m1,m2
real :: diff(3,3), omega, s,q, qv(3),qnv(3)
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
implicit none
integer :: i,k,l,n,nn
real :: diff(3,3)
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
real function strfak(kapa) !calculates formfactor by interpolation  
real :: kapa
end function
real function ScatFunction(kapv,omega)
real :: kapv(3), omega
integer :: n
end function
subroutine SuscepComponents(z)
real :: z
end subroutine
end interface

erg(1)='0.dat'
erg(2)='1.dat'
erg(3)='2.dat'
erg(4)='3.dat'
erg(5)='4.dat'
erg(6)='5.dat'


do n=0,5
if (mode==n) then
  i=index(outfilename,'.')
  if (i==0) then 
    f=trim(outfilename)//erg(n+1)
  else
    f=outfilename(1:i-1)//erg(n+1)
  end if
end if  
end do

write(*,'(A)') 'results written into',  f
if (mst == 1) then  
open(21,file=f,action='write')
end if
if (mst == 2) then
open(21,file=f,action='write',position='append')
end if

write(21,*) '#results of program bfk using leve-scheme ', cefname
write(21,'(A14,F6.2,A18,F6.2)') '#temperature ', temp, ' coupling constant ', g
write(21,'(A33,3G10.3)') '#magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)



if (mode==0) then
  write(21,*) '#mode =', mode , 'Calculates S(i,i)=Im Chi(i,i)*coth(beta*omega/2)'
  write(21,*) '#and checks sum-rule (S(1,1)+S(2,2)+S(3,3))/Pi = J*(J+1)' 
  write(21,*) '#energy loss,    S(1,1)         S(2,2)       S(3,3)'

  deltax=Emax/Npoints
  s=0
  m=4
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
    s=s+(ss(1,1)+ss(2,2)+ss(3,3))*deltax
  end do
  sum=s/Pi
  write(21,*)'#Integral (ss(1,1)+ss(2,2)+ss(3,3)/Pi: ', sum, ( Ns**2-1.)/4 

else

  if ((mode==1) .or. (mode==2)) then

    write(21,*) '#strfactor ', strfilename 
    write(21,*) '#scattering mode m=  ', mode   
  end if  
  if (mode==1) then
    open(22,file=tablename, action='read')
    write(21,*) '# Calculation of scattering cross section for given set of '
    write(21,*) '# scattering wave vector q and energy loss'
    write(21,*) '# provided by input-file ', tablename
    write(21,*) '# q (A^-1) energy loss (meV) cross section (1/meV,steradian)'

    calculation: do
      read(22,*,iostat=status_read) kapv(1), kapv(2), kapv(3), enloss
      if (status_read /= 0 ) exit
      kappa=sqrt(kapv(1)**2+kapv(2)**2+kapv(3)**2)
      do i=1,3
        kapnv(i)=kapv(i)/kappa
      end do
      st=strfak(kappa)
      dcs=(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)

      write(21,'(3G14.5)') kappa,enloss,dcs
    end do calculation
    close(22)
  end if      

!-------------------------------------------

    if (mode==2) then
       write(21,*) '#calculation of different components of the scattering functions'
       write(21,*) '#as function of q and omega provided by input-file ',tablename 

       write(21,*) '#Q,  energy-loss '
write(21,*) '#                      S(1,1),S(1,2),S(1,3)'
write(21,*) '#                      S(2,1),S(2,2),S(2,3)'
write(21,*) '#                      S(3,1),S(3,2),S(3,3)'

open(22,file=tablename,action='read')
status_read=0
calc: do 
    read(22,*,iostat=status_read) kapv(1),kapv(2),kapv(3),enloss
    if (status_read /= 0 ) exit
    kappa=sqrt(kapv(1)**2+kapv(2)**2+kapv(3)**2)
    do i=1,3
      kapnv(i)=kapv(i)/kappa
    end do
    st=strfak(kappa)
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
end do calc
close(22)
end if      
!-----------------------------------------------------

if (mode==3) then
write(21,*) '# calculation of cross section for fixed direction k and energy E0'
write(21,*) '# of incoming beam, fixed direction but variable length of scattering'
write(21,*) '# wave vector q '
write(21,'(A)') ' #E0,    k1,   k2,   k3;   q1,   q2,   q3 '
write(21,'(7F6.2)') E0,k1,k2,k3,kap1,kap2,kap3
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
    dkapv(i)=kapv(i)/200
  end  do

  do n=1,200
     do i=1,3
       xkapv(i)=n *dkapv(i)
       ksv(i)= kv(i)-xkapv(i)
     end do
     kappa=sqrt(xkapv(1)**2 + xkapv(2)**2 + xkapv(3)**2)
     do i=1,3     
       kapnv(i)=xkapv(i)/kappa
     end do     
     ksq=ksv(1)*ksv(1)+ksv(2)*ksv(2)+ksv(3)*ksv(3)
     ks=sqrt(ksq)
     st=strfak(kappa)
     enloss=ek2*(k0**2-ksq)
     dcs=ScatFunction(xkapv,enloss)
     dcs=ks/k0*(r0*gl/2*st)**2/Pi*dcs
      write(21,'(3G14.5)') kappa,enloss,dcs
   end do
end if

!-----------------------------------------------------------------------
if (mode==4) then
write(21,*) '# calculation of cross section for fixed direction k and energy E0'
write(21,*) '# of incoming beam, fixed direction ks but variable energy of '
write(21,*) '# scattered beam. vector q = vector k0 - vector ks '

write(21,'(A)') ' #E0,    k1,   k2,   k3;   ks1,  ks2,  ks3 '
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
    dksv(i)=ksv(i)/200
  end  do
    do n=1,200
      do i=1,3
        xksv(i)=float(n)*dksv(i)
        kapv(i)=-xksv(i)+kv(i) ! scattering vector
      end do
      ksq=xksv(1)**2+xksv(2)**2+xksv(3)**2
      ks=sqrt(ksq)
! length of scattering vector:
      kappa=sqrt(kapv(1)**2+kapv(2)**2 + kapv(3)**2+0.00001)
      do i=1,3
         kapnv(i)=kapv(i)/kappa ! normalised length of scattering wave vector
      end do
      st=strfak(kappa)
      k0=sqrt(E0/ek2)
! energy loss
      enloss=ek2*(k0**2-ksq)
!differential crosssection 
      dcs=ks/k0*(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)
      write(21,'(3G14.5)') kappa,enloss,dcs
  end do
end if
!-----------------------------------------------
if (mode==5) then
!write(21,'(A32,3G10.3)') '#magnetisation <Jx>, <Jy>, <Jz> ', jav(1), jav(2), jav(3)
write(21,*) '# calculation of cross section for fixed direction ks and energy Es'
write(21,*) '# of scattered beam, fixed direction k but variable energy of '
write(21,*) '# incoming beam. vector q = vector k - vector ks '


write(21,'(A)') ' #Es,     k1,    k2,   k3;   ks1,  ks2,  ks3 '
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
    dkv(i)=kv(i)/200
  end  do

  do n=1,100
    if (n /= 0 ) then
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

      st=strfak(kappa)
      ksq=ks**2
      kq=k0**2
      enloss=ek2*(kq-ksq)
!differential crosssection
      dcs=ks/k0*(r0*gl/2*st)**2/Pi*ScatFunction(kapv,enloss)
      write(21,'(3G14.5)') kappa,enloss,dcs
    end if
  end do
end if

close(21)
end if
end subroutine OutputResults
!--------------------------------------------


program bft
 
call DataRead
call MatrixElementCalculation
call OutputResults

end program
!-------------------------------------------------------------------


