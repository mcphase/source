      program rcn36
c     program rcn36(in36,tape10=in36,out36,tape9=out36,tape11,
c    1  indiel,tape12=indiel,tape4,tape2n,tape2=tape2n,tape7)
c
c          This program copyright by los alamos national laboratory
c               and the u.s. government.  the program may be freely
c               distributed gratis anywhere in the world, but
c               may not be sold.
c
c
c     hartree-fock-slater self-consistent atomic field program(kut=1,0)
c          with options for hartree (kut=-1, exfm1=0.0) or for
c               hartree-plus-statistical-exchange (kut=-1, exfm1.gt.0)
c               or for hartree-slater (kut=-2)
c               or for hartree-fock (ihf=2)
c
c      originally  written by sherwood skillman,
c           rca laboratories, princeton, new jersey, spring 1961
c      modified by frank herman, summer 1961
c      further modified by richard kortum,  lockheed research
c           laboratories, palo alto, california,  summer 1962
c      further mod by r. d. cowan, lasl, jan-mar1964, oct 65-mar 66,
c           may-june 1966, nov-dec 1966.
c      further mod by k. l. andrew, lasl, may 1968, to simplify
c           configuration input
c      method of calculating correlation-energy correction
c           modified january, 1971.
c      calculation of continuum functions modified june, 1972
c           to handle kinetic energies up to 500 rydbergs
c           (no mesh-size doubling beyond point idb at which
c                deltar=0.25/sqrt(emx)   )
c      finite-boundary and relativistic options added spring,1973
c      hartree-fock options (ihf=2) added april and august 1981
c      more flexible configuration-card input introduced (mod 33)
c           at purdue, july-august 1987
c      Vinti-integral option added at Royal Holloway Bedford New College,
c           (England) october 1987
c      option to delete writing of output on tape2 when the configuration
c           card contains an asterisk anywhere in columns 21 to 26
c           added at the Zeeman Laboratory (Amsterdam), december 1987.
c
c      rcn mod 34, adapted to the vax 8200 at Lund University (Sweden),
c           made march 1988.  further modified to use generic functions,
c           etc., at LANL, july-december, 1988
c      rcn mod 35, incorporating subroutines
c           diel and analyz1 to set up input for dielectronic-
c           recombination runs, made at LANL january 1989
c      This version (mod 36) with commons of the same name having the
c           same length in all routines, and with timing routine
c           SECONDS coded for CRAY, VAX and MACINTOSH, SUN, or
c           IBM System/6000 RISC, made at LANL February-April 1992
c               (uncomment appropriate section of this routine).
c      Modifications for Macintosh Centris (using Language Systems
c           Fortran), and feature to write some output to the monitor
c           screen if IW6<0, added May 1993.
c      Vers. 36.0.1 (August 1996): Minor change at rcn3s 930+ and 
c           a(5,5) -> a(ko,5) in blank common, to remove unimportant
c           problems when compiled with bounds checking.
c      Vers. 36.0.2 (September 28, 1996): Additions before statement
c           200 of subroutine scheq to cure a bug for irel=1, by
c           smoothing out a singularity in the relativistic potential
c           from the term  (dP/dr)/P -1/r in Eq. (7.60) of my book.
c      Vers. 36.0.3 (March 1997): For neutral and singly ionized 
c           lanthanide and actinide systems, change default values of 
c           alfmin and alfmax from 0.2 and 1.0 to 0.15 and 0.75;  
c           at rcn 540, change way in which alfmin and alfmax are 
c           modified.  Also, modify code in region 710-730.
c      Vers. 36.1.0 (June 1998).  All equivalence statements removed.
c      Vers. 36.1.1 (Dec 1998).  Blank commons changed to labeled.
c      Vers. 36.1.2 (May 1999).  Removed numerous unused variables that
c           give warning diagnostics with some compilers.
c      Vers. 36.2.0 (Aug 1999).  Added the possibility of using fractional
c           occupation numbers (maximum of three decimal places; i.e.,
c           maximum of 4 to 6 characters, according as the integral part
c           of the occupation number is 0, 1, or 2 digits, respectively).
c         NOTE: If any occupation number is non-integral, nothing will
c           be written on the file tape2n for that "configuration", 
c           as fractional-occupation-number results are not appropriate 
c           for either rcn2 or rcg.
c                Also changed the action of the diagnostic variable npr,
c           improved the algorithm for normalizing continuum functions,
c           and added warning prints to the monitor if the value of emx
c           is inappropriate for highly excited or continuum functions.
c      Vers. 36.2.1 (Feb 9, 2000)  Corrected errors in calculation of
c           Breit energies when wavefunctions are printed (norbpt.gt.0).
c           If irel.gt.2, print qnl at mesh points 1, 11... instead of
c           pnl at points 6,16..., and print diagnostic output for irel
c           gt.3 instead of gt.2.
c      Vers. 36.2.2 (Feb 28, 2000).  Added line of print to monitor output
c           giving Z, no. of electrons, ionization stage, IREL, etc.
c      Vers. 36.2.3 (Feb 13, 2002).  Corrected a bug in subroutine scheq
c           (3 statements after 248) that prevented convergence on
c           orbitals with large l (greater than 6 or so).
c      Vers. 36.2.4 (Apr 5, 2002).  Revised subroutine diel, which had 
c           not been made compatible with the 1998 changes.
c      Vers. 36.2.5 (Mar 19,29, 2004). Modified statement before 708 to keep
c           divisor from being very small; also increased test from
c           gt.5.0e-06 to gt.5.0e-04.  After 700 of scheq, added pnlo(k)=0.0.
c
c               -----------------------------------------------
c
c
c     Binary tape 2  contains self consistent atomic potential,
c      energy eigenvalues, and radial wave functions
c      successive solutions separated by end of logical record.
c     z = atomic number.  ion = ionicity.  zzz = ion+1.
c     rho=total electron density
c     rhom=electron density of min(wwnl(m),2.0) electrons in orbital m
c     ruee=total classical electron-electron potential
c     rueem=same for one electron of orbital m
c     ruexch=-2.0*(3.0*rho/pi)**(1.0/3.0)
c     ru=-2.0*z+ruee+exf10*ruexch
c     kut= 1, rv=min1f(ru,-2.0*ion) ---hfs (no tail cutoff)
c     kut= 0, rv=min1f(ru,-2.0*zzz) ---hfs (with latter tail cutoff)
c     kut=-1, ---hx or hartree
c       rv=-2*z+ruee-rueem+exfm1*ruexch*(1.0-rhom/rho)*f1*f2
c          f1=(rho-rhom)/(rho-rhom+ca0/(n-l))
c          f2=1.0 except for l(m) greater than 1 and either
c            l(m-1)=l(m), or l(m-2)=l(m) and wwnl(m-1)=1,
c                 in which cases
c          f2=1.0+ca1*(r(j)-r(i))/r(j)  for i less than j,
c            where r(j) is the position of the kth node of pnl(m)
c               (k=no. of orbitals with l=l(m) and n less than n(m).)
c     ihf=1, write output on tape7 as input to program hf8
c     ihf=2, make self-contained hartree-fock calculation within rcn
c
c     maxit = maximum no. of iterations
c     exf=factor multiplying slater exchange term
c     corrf=factor multiplying correlation term
c
c
      character*20 inf
c
c          ----------------------------------------------------------
c          Comment out IMPLICIT REAL*8 cards for 64-bit machines,
c               uncomment them for VAX, SUN, etc.
c          ----------------------------------------------------------
c          Open files, if not done via program card
c               (Use iii>1 and program card at start
c                of program for cray and cdc computers;
c                otherwise, comment out program card that defines files and
c                set iii=1 to use standard file names, or
c                set iii=0 to name files interactively, e.g. for vax)
c           ---------------------------------------------------------
c
      iii=1
      if (iii.eq.1) go to 150
      if (iii.gt.1) go to 160
c
      write(6,*) ' input file name  '
      read(5,11) inf
   11 format(a20)
      open(10,file=inf,status='old')
      write(6,*) ' output file name  '
      read(5,12) inf
   12 format (a)
      open(9,file=inf,status='unknown')
      write(6,*) ' give name of data file for the rcn2-program '
      read(5,12) inf
      if(inf.ne.'   ') go to 120
      open(2,status='scratch',form='unformatted')
      go to 155
  120 open(2,file=inf,status='unknown',form='unformatted')
      go to 155
c
c          use these open statements for standard file names
c               (set iii=1 in statement above)
c
  150 open(10,file='in36',status='old')
      open(9,file='out36',status='unknown')
      open(2,file='tape2n',status='unknown',form='unformatted')
  155 continue
      open(7,status='scratch')
      open(4,status='scratch',form='unformatted')
c
c
  160 rewind 10
      read (10,10,end=200) inf(1:5)
   10 format (a5)
      if (inf(1:5).ne.'   -1') go to 300
c
c          if tape10 is empty or the first record is an exit card, 
c               call diel to set up rcn input for dielectronic recombination runs
c
  200 if (iii.lt.2) open(12,file='indiel',status='old')
      if (iii.lt.2) open(11,status='scratch')
      call diel
c
c          call rcn for (final) calculation of radial wavefunctions
c
  300 rewind 10
      call rcn
      end
      subroutine rcn
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      common/lc2/ra(kmsh),rb(kmsh)
c
      common/c9/ qq(kmsh),delp(ko),alfm(ko)
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      character*6 wwwnl(25),wwnl8(8)
      common/charww/wwwnl,wwnl8
      integer fg
c
      dimension ru0(91),deltn(200)
      character*8 chxid
      data hxid1,hxid2,hfid1,hfid2 /2hxb,2hxr,2hhf,2hhr/
      data ru0/1.0,0.99094,0.98164,0.97218,0.96266,0.95311,0.94359,
     1  0.93413,0.92475,0.91549,0.90634,0.88844,0.87109,0.85429,
     2  0.83802,0.82227,0.80700,0.79220,0.77785,0.76393,0.75044,
     3  0.72473,0.70066,0.67816,0.65710,0.63732,0.61866,0.60097,
     4  0.58416,0.56814,0.55285,0.52427,0.49813,0.47416,0.45214,
     5  0.43184,0.41304,0.39553,0.37916,0.36377,0.34927,0.32262,
     6  0.29887,0.27789,0.25953,0.24348,0.22926,0.21645,0.20468,
     7  0.19376,0.18354,0.16493,0.14852,0.13412,0.12152,0.11053,
     8  0.10094,0.09255,0.08517,0.07866,0.07289,0.06311,0.05556,
     9  0.048,0.042,0.037,0.033,0.030,0.025,0.021,0.017,0.014,0.011,
     a  0.009,0.007,0.006,0.005,0.004,0.003,0.002,0.001,0.0,0.0,0.0,
     b  0.0,0.0,0.0,0.0,0.0,0.0,0.0/
      data jhlf/102/
      data jhuF/70/
c
c      kmsh=1801
c      ko=20
c
c
c     position tapes 2, 4, and 7 for next run
c
  160 rewind 2
      rewind 4
      rewind 7
      meshx=kmsh
      mesh=meshx
      izo=-7
      iono=-7
      do 165 m=1,20
      eeo(m)=-7.0
      wwnlo(m)=0.0
  165 nlbcdo(m)=' '
      ncspvo=0
      ns=1
      nconf=0
      nconft=0
      nconft2=0
      rsint1=0.0
      do 170 i=1,mesh
      ruexch(i)=0.0
  170 ruee(i)=0.0
c
c          expand atomic potential from data block
c
  195 m=0
      do 196 i=1,361,4
      m=m+1
  196 ru(i)=ru0(m)
      call seconds(t0)
      t0=t0
      t1=t0
c
c          read control card.
c
  200 read (10,20) itpow,iptvu,ipteb,norbpt,izhxbw,iphfwf,ihf,ibb,
     1  tolstb,tolkm2,tolend,thresh,
     2  kutd,kut1,ivinti,irel,maxit,npr,exf10,exfm1,emx,corrf,iw6
   20 format (3i1,i2,i1,2i2,i3,f2.1,e5.1,2e10.1,i2,i2,2i1,2i2,4f5.5,i5)
      ca0=0.50
      ca1=0.70
      ihf1=ihf
      if (ihf.gt.1) ihf=0
      thres1=thresh
      itemp1=kmshq-kmsh
      itemp2=koq-ko
      if ((itemp1.lt.0.or.itemp2.lt.0).and.irel.ge.2) irel=1
      hxid=hxid1
      if (irel.ne.0) hxid=hxid2
      if (ihf1.gt.1) hxid=hfid1
      if (ihf1.gt.1.and.irel.ne.0) hxid=hfid2
      chxid='        '
      if (ihf1.gt.1) chxid='hf      '
      if (ihf1.gt.1.and.irel.ne.0) chxid='hfr     '
      if (kutd.eq.0) kutd=-2
      maxit0=maxit
      if (emx.gt.0.0) go to 204
c          set emx = largest continuum energy (if any)
      iread=0
c 201 read (10,19) n,iz,ee8
c  19 format (i1,i4,69x,0pf6.5)
  201 read (10,21) n,iz,ibb1,nspect,inpcard
      call analyz
      iread=iread+1
      if (iz.le.0) go to 202
      emx=max(emx,ee8)
      if (n.eq.0) go to 201
      do 1201 i=1,n
      read (10,18) nfg
   18 format (10x,i5)
      if (nfg.gt.3) read (10,18)
      if (nfg.gt.7) read (10,18)
 1201 iread=iread+1+nfg/4
      go to 201
c
  202 if (iread.eq.0) go to 204
      do 203 i=1,iread
  203 backspace 10
c
  204 if (emx.le.0.0) emx=1.0e-20
      mesh=meshx
      if (emx.le.1.0e-20) mesh=min0(mesh,641)
      ib=ibb
      if (ibb.eq.0.and.mesh.lt.meshx) ib=mesh
      if (ibb.eq.0.and.ib.eq.0) ib=mesh-1
      icut=mesh
      if (irel.eq.0) go to 210
      do 205 i=1,mesh
  205 vr(i)=0.0
c
c          read z, stage of ionization, and configuration.
c
c 210 read (10,21) nhftrm,iz,ibb1,nspect,conf,ialfmin,ialfmax,
c    1  (nlbcd8(m),wwnl8(m), m=1,8), kut,ee8
c  21 format (i1,i4,i3,i2,3a6,2i2,8(a3,a2),i2,0pf6.5)
c     alfmin=0.01*ialfmin
c     alfmax=0.01*ialfmax
  210 read (10,21) nhftrm,iz,ibb1,nspect,inpcard
   21 format (i1,i4,i3,i2,a70)
      n1scf=1
      if (ibb1.le.ko.and.nconft.gt.0) n1scf=max0(1,ibb1+1)
      if (ibb1.le.ko) ibb1=0
      if (iz.eq.-999) return
      if (iz) 212,200,214
  212 end file 2
      rewind 2
      end file 9
      rewind 9
      if (ihf.ne.0) go to 213
      stop '(normal exit)'
  213 ncspvs=0
      if (rhohf.eq.0.0) rhohf=-3.0
      if (h.eq.0.0) h=0.065
      if (no.eq.0) no=200
      if (nufsh.eq.0) nufsh=0
      write (7) conf,nconf,iz,z,ncspvs,rhohf,h,no,nhftrm,irel,etotrl,
     1  nufsh,izhxbw,ns,iphfwf
      end file 7
      rewind 7
      stop '(cff hf exit)'
  214 if (inpcard(1:6).ne.'      ') go to 216
c          skip over iz completed runs on tape 2
      do 215 n=1,iz
      read (2) nnn,l,wwnl,ncspvs,vprm,conf,nconf,iz1,z,ion,irel,hxid,
     1  mesh,c,idb,exf,corrf,ca1,ca0,kut,npar,nlbcd,nnlz,nkkk,ee,norb,
     2  r,ru
  215 write (9,22) nconf, (conf(i), i=1,3)
   22 format (36h0on tape 2, skip over configuration ,i3,1h=,3x,3a6)
      nconft=nconf
      nconft2=nconf
      go to 210
  216 if (ibb1.gt.0) ib=ibb1
      if (nspect.eq.0)  nspect = 1
      noelec = iz + 1 - nspect
c          analyze inpcard
      call analyz
c
      nf(1)=0
      ng(1)=0
      if (ihf.eq.0) go to 217
      if (irel.gt.0.and.noelec.gt.2) go to 217
      if (tolend.gt.5.e-5) tolend=5.e-8
      if (thresh.gt.1.e-8) thresh=1.e-11
      thres1=thresh
  217 if (nhftrm.eq.0) go to 230
      do 220 i=1,nhftrm
      read (10,23) t1hf(i),t2hf(i),nfg,
     1  (cfg(j,i),fg(j,i),kfg(j,i),ifg(j,i),jfg(j,i), j=1,3)
   23 format (a7,a3,i5,5x,3(f9.8,a1,i1,1x,a3,1x,a3,1x))
      if (nfg.gt.3) read (10,24)
     1  (cfg(j,i),fg(j,i),kfg(j,i),ifg(j,i),jfg(j,i), j=4,nfg)
   24 format (4(f9.8,a1,i1,1x,a3,1x,a3,1x))
      nf(i)=0
      do 218 j=1,nfg
      if (fg(j,i).eq.jhlf.or.fg(j,i).eq.jhuF) nf(i)=nf(i)+1
  218 continue
  220 ng(i)=nfg-nf(i)
  230 if (alfmin.eq.0.)  alfmin=0.20
      if (alfmax.eq.0.)  alfmax=1.0
      alfmax0=alfmax
      alffix=0.4*alfmax
      if (kut.eq.-1) alffix=alfmax
      if (kut.eq.0)  kut = kut1
      if (iz.ne.izo.and.kut.lt.0) kut=0
      maxit=maxit0
      z=iz
      twoz=z+z
      nvales=1
      c=0.88534138/z**(1.0/3.0)
      nblock=(mesh)/40
c
      if (itpow.lt.8) go to 240
      write (9,26) iz,nspect,conf,alfmin,alfmax, (nlbcd8(m),wwnl8(m),
     1  m=1,8), ee8,alffix
   26 format (29h1card image except for alphas,2(2x,i3),3a6, f5.3,1x,
     1  f5.3, 8(a3,a6),f8.5, 8halffix =,f5.3 //  )
c
c           set up orbitals and estimate eigenvalues
c
  240 call setcfg
      if (ierror.ne.0) go to 210
      if (ion.gt.1) go to 250
      if (iz.lt.55.or.iz.gt.102) go to 250
      if (iz.gt.70.and.iz.lt.87) go to 250
      temp=0.75
      alfmax=min(alfmax,temp)
      alfmin=0.2*alfmax
      alffix=0.4*alfmax
      alfmax0=alfmax
      if (kut.eq.-1) alffix=alfmax
  250 if (ihf.eq.0.or.nnn(ncspvs).le.10) go to 260
      if (tolend.gt.5.e-5) tolend=5.e-8
      if (thresh.gt.1.e-8) thresh=1.e-11
      thres1=thresh
  260 if (n1sc.gt.1.or.n1scf.gt.1) go to 300
c
c          construct atomic potential
c
  270 if (nconft.gt.0) go to 280
      do 272 i=1,357,4
  272 ru(i)=-ru(i)*twoz
      ru(361)=ru(357)
      ru(365)=ru(357)
      m=9
      do 275 i=1,357,4
      m=m-1
      if (m.ge.0) go to 274
      ru(i+1)=(22.0*ru(i)+11.0*ru(i+4)-ru(i+8))/32.0
      ru(i+2)=(10.0*ru(i)+15.0*ru(i+4)-ru(i+8))/24.0
      ru(i+3)=( 6.0*ru(i)+27.0*ru(i+4)-ru(i+8))/32.0
      m=9
      go to 275
  274 ru(i+1)=(21.0*ru(i)+14.0*ru(i+4)-3.0*ru(i+8))/32.0
      ru(i+2)=( 3.0*ru(i)+ 6.0*ru(i+4)-    ru(i+8))/ 8.0
      ru(i+3)=( 5.0*ru(i)+30.0*ru(i+4)-3.0*ru(i+8))/32.0
  275 continue
      do 276 i=362,mesh
  276 ru(i)=ru(361)
      go to 285
c
  280 if (abs(twoz+ru(1)).le.0.001) go to 285
      eras=(-twoz-ru(mesh))/(ru(1)-ru(mesh))
      do 282 i=1,mesh
  282 ru(i)=ru(mesh)+eras*(ru(i)-ru(mesh))
  285 rum=-twoion
      if (ru(mesh).le.rum) go to 290
      eras=(rum-ru(1))/(ru(mesh)-ru(1))
      do 289 i=2,mesh
  289 ru(i)=ru(1)+eras*(ru(i)-ru(1))
c
c          construct r mesh
c
  290 i=1
      r(1)=0.0
      idb=1
      rdb=0.0
      deltar=0.0025*c
      drx=0.40/sqrt(emx)
      do 295 j=1,nblock
      do 294 k=1,40
      i=i+1
      r(i)=r(i-1)+deltar
  294 continue
      if (deltar.ge.drx) go to 295
      idb=i
      rdb=r(i)
      deltar=deltar+deltar
  295 continue
      r(1)=1.0e-10
      v(1)=-twoz/r(1)
c
c
  300 ihftrm=max0(1,nhftrm)
      if (ihf.eq.1) ihftrm=1
      write (nla,29) conf(1)
   29 format (a6)
      ndel=1
      do 305 k=1,6
      if (nla(k:k).eq.'*') ndel=0
  305 continue
      do 850 ii=1,ihftrm
      nfg=nf(ii)+ng(ii)
      nconft=nconft+1
      nconft2=nconft2+ndel
      nconf=0
      if (ndel.eq.1) nconf=nconft2
      ion=iz-noelec
      delta=2000.0
      delt=2000.0
      exf=exf10
      if (kut.eq.-1) exf=exfm1
  310 niter=0
      end=0.0
c
  330 if (tolend.lt.1.0e-04) write (9,30)
   30 format (1h1)
      write (9,31) conf,nconf,iz,ncores,nvales,ion,c,ca1
   31 format (3h   ,3a6,2x,10h    nconf=,i3,6h    z=,i3,11h    ncores=,
     1  i2,11h    nvales=,i2, 8h    ion=,i2, 14h    r(i)/x(i)=,f10.8,
     2  6x,4hca1=,f5.3)
      n1scfm1=n1scf-1
      n1sc=n1scf
      if (itpow.eq.2) go to 390
      write (9,32)tolstb,tolkm2,tolend,thresh,kut, exf,corrf,ca0,ihf1,
     1  mesh,idb,rdb,emx,r(mesh),irel,ib
   32 format (/8h tolstb=,f5.3,3x,7htolk-2=,f7.5,3x,7htolend=,e8.1,
     1  3x,7hthresh=,e8.1,3x,4hkut=,i2,4x,4hexf=,f5.3,9h   corrf=,f5.3,
     2  6x,4hca0=,f6.3//18h0rcn mod 36  ihf1=,i2,
     3  8x,5hmesh=,i4,3x,4hidb=,i4,3x,4hrdb=,f8.3,3x,4hemx=,f7.2,
     4  3x,8hr(mesh)=,f8.3,5x,5hirel=,i2,5x,3hib=,i4)
      if (ihf1.ne.2) go to 360
      if (nhftrm.gt.0) go to 345
      write (9,33) chxid
   33 format (1h0,a8,15h energy = e(av))
      go to 360
  345 write (9,34) t1hf(ii),t2hf(ii),(cfg(j,ii),fg(j,ii),kfg(j,ii),
     1  ifg(j,ii),jfg(j,ii), j=1,nfg)
   34 format (12h0energy = e(,a7,a3,9h) = e(av),4(3h + ,f9.6,1h*,a1,i1,
     1  1h(,a3,1h,,a3,1h))
     2  /31x,4(3h + ,f9.6,1h*,a1,i1,1h(,a3,1h,,a3,1h)))
  360 if (n1sc.gt.1) write (9,36)
   36 format (50x,45h calculate only outermost (continuum) orbital)
      if (n1scf.le.1) go to 370
      write (9,27) n1scfm1
   27 format ('0*****first',i3,' orbital(s) held fixed')
c
  370 write (9,37) (nlbcd(m),m=1,ncspvs)
   37 format (/31h it  time   delta  idel kut itp,1x,15(a3,3x)/
     1  34x,15(a3,3x))
      if (ifrac.eq.0) then
        write (9,38) (wwnl(m),m=1,ncspvs)
      else
        write (9,39) (wwnl(m),m=1,ncspvs)
      end if
  390 if (iw6.ge.0) go to 400
      write (6,1039) (conf(i),i=1,3),iz,ion,noelec,irel,emx,idb,mesh,
     1  r(mesh)
 1039 format (//1x,3a6,' Z=',i3,'   Ion stage=',i3,'   Noelec=',i3,
     1  '   IREL=',i1,'   EMX=',f6.2,'   IDB=',i4,
     2  '   R(',i4,')=',f6.1)
      write (6,37) (nlbcd(m),m=1,ncspvs)
      if (ifrac.eq.0) then
        write (6,38) (wwnl(m),m=1,ncspvs)
      else
        write (6,39) (wwnl(m),m=1,ncspvs)
      end if
   38 format (30x,15f6.0/32x,15f6.0)
   39 format (31x,15f6.3/33x,15f6.3)
c
c          start scf iteration
c
  400 niterp=0
      if (kut.eq.-1) niterp=1
  401 niterp=niterp+1
  402 niter=niter+1
      thresh=thres1
      if (alffix.ge.0.7) go to 406
c     if (delta/tolend.lt.100.0) go to 406
c     thresh=sqrt(thres1*delta)
c 406 if (npr.le.0) go to 409
  406 go to 409
c             This diagnostic removed 22 Aug 1999
c     write (9,40) niter, (ru(i), i=1,400,20), (ruee(i), i=1,400,20),
c    1  (ruexch(i), i=1,400,20)
c  40 format (1h1,i2,5x,14hru,ruee,ruexch/(10f12.6/10f12.6//))
  409 do 410 i=1,mesh
      rscore(i)=0.0
  410 rsvale(i)=0.0
      nod=0
      nof=0
c
c          calc v for kut greater than -1
c
  411 do 413 i=1,mesh
  413 rc(i)=rsatom(i)
      eras=-twozzz
      if (kut) 430,418,417
  417 eras=-twoion
  418 do 425 i=2,mesh
      a0=min(ru(i),eras)
      v(i)=a0/r(i)
  425 x2(i)=v(i)
c
c     solve schroedinger equation for each orbital in turn.
c            also, calculate core and valence charge densities.
c
  430 do 499 m=1,ncspvs
      if (end.gt.0.0.and.m.gt.ncspvs) go to 800
      kkk=nkkk(m)
      if (end.eq.0.0.and.n1sc.gt.1.and.m.lt.n1sc) go to 491
      e=ee(m)
      nn=nnn(m)
      lam=l(m)
      nodf=0
      if (lam-2) 433,431,432
  431 nod=nod+1
      nodf=nod
      go to 433
  432 if (lam.ne.3) go to 433
      nof=nof+1
      nodf=nof
  433 if (e.le.0.0) go to 434
      nn=99999
  434 if (kut.lt.0) go to 435
      if (end) 450,450,449
  435 if (m.gt.ncspvo) go to 450
      do 437 i=1,mesh
      if (e.le.0.0) go to 436
      rsatom(i)=0.0
      go to 437
  436 rsatom(i)=pnl(i,m)**2
  437 continue
c
c     calc xi(i)=coulomb energy of one electron in orbital m
c
      call quad2(1)
c
c          if (ihf1.gt.1) calc hf potential
      if (ihf1.gt.1) call hfpot(m)
      two=2.0
      a4=min(wwnl(m),two)
      j=jjj(m)
      do 447 i=2,mesh
      if (ihf1.gt.1.and.irel.eq.0) go to 446
      if (kut.eq.-2) go to 442
c          calc hx potential
      eras=rc(i)-a4*rsatom(i)
      a0=exfm1*eras/(eras+a0m(m))
      if (i.ge.j) go to 441
      a0=a0*(1.0+ca1*(r(j)-r(i))/r(j))
  441 continue
      a2=0.0
      if (rc(i).ne.0.0) a2=eras/rc(i)
      x1(i)=ru(i)-xi(i)+(a0*a2-exf10)*ruexch(i)
      go to 444
c          calc hs potential
  442 eras1=r(i)
      eras=a4*rsatom(i)
      eras=2.0*rsatom(i)
      call subcor
      x1(i)=ru(i)-xi(i)-exf10*b3
      if (wwnl(m).gt.1.0) go to 444
      eras=rc(i)+rsatom(i)
      call subcor
      x1(i)=x1(i)+exf10*(b3-ruexch(i))
  444 continue
      x2(i)=x1(i)/r(i)
      if (ihf1.le.1) xj(i)=x1(i)
  446 v(i)=xj(i)/r(i)
c          add correlation contribution to potential
      if (corrf.le.0.0) go to 447
      eras1=r(i)
      eras=rc(i)-rsatom(i)
      call subcor
      v(i)=v(i)+corrf*b2/eras1
  447 continue
      if (end.le.0.0) go to 450
      kkk=nkkk(m)
      if (ihf.eq.1.and.irel.eq.0.and.ncspvs.gt.1) go to 449
      if (iptvu.lt.5) go to 449
      eras=delp(m)/pnl(kkk-47,m)
      eras1=l(m)*(l(m)+1)
      do 448 i=1,400,20
  448 rsatom(i)=xj(i)+eras1/r(i)
      write (9,49) niter,nlbcd(m),rsint1,eras, (xi(i), i=1,400,20),
     1  (xj(i), i=1,400,20), (rsatom(i),i=1,400,20)
      if (iptvu.lt.6.or.m.lt.ncspvs) go to 449
      write (9,1570) conf,nlbcd(m),ee(m)
 1570 format (1h1,3a6,1h(,a3,10h function,,3x,2he=,f11.6,1h)///5x,1hi,
     1  10x,1hr,13x,3hr*v,10x,6hr*veff,12x,1hp,9x,10h(veff-e)*p,3x,
     2  14hd2p/dr2 (5 pt),1x,14hd2p/dr2 (3 pt),5x,6hveff-e/)
      ix=min0(nkkk(m),500)
      do 1571 i=2,ix
      rsatom(i)=xj(i)+eras1/r(i)
      a4=rsatom(i)/r(i)-ee(m)
      b0=a4*pnl(i,m)
      if (i.eq.2.or.mod(i,40).eq.0) go to 1569
      j=i+2
      jj=40*(i/40)+1
      a0=1.0/(r(j)-r(j-1))**2
      do 1568 i1=1,5
      pnlo(553-i1)=pnl(j,m)
      j=j-1
      if (j.lt.jj) j=j-1
 1568 continue
      b2=a0*(16.0*(pnlo(551)+pnlo(549))-30.0*pnlo(550)-pnlo(552)
     1  -pnlo(548))/12.0
      b1=a0*(pnlo(551)-2.0*pnlo(550)+pnlo(549))
      go to 1571
 1569 b2=0.0
      b1=4.0*a0*(pnl(i+1,m)-2.0*pnl(i,m)+pnl(i-1,m))
 1571 write (9,1572) i,r(i),xj(i),rsatom(i),pnl(i,m),b0,b2,b1,a4
 1572 format (i6,7f15.6,1p,e19.6)
  449 write (4) (v(i), i=1,mesh)
      go to 499
c
c     integrate differential equation to obtain radial wavefunction pnl
c
  450 mm=m
      if (m.lt.n1scf) go to 491
      call scheq
      jjj(m)=j
      if (m.lt.3) go to 451
      if (l(m).eq.l(m-1)) go to 452
      if (l(m).ne.l(m-2)) go to 451
      if (wwnl(m-1).eq.1.0) go to 452
  451 jjj(m)=1
c
  452 ee(m)=e
      ehf(m)=-e
      nkkk(m)=kkk
      kkk=nkkk(m)
      imat(m)=imatch
      if (niter.gt.10) go to 454
      nsch(niter,m)=nprint
c
c          stabilize pnl if e neg and niterp gr than 1
c
  454 eras=0.0
      eras=pnlo(kkk-47)-pnl(kkk-47,m)
      if (niterp.eq.1) go to 460
      if (e.gt.0.0) go to 459
      if (eras.eq.0.0) go to 461
      if (niterp-2) 460,470,480
  459 alfm(m)=1.0
      go to 461
  460 alfm(1)=alffix
      alfm(m)=0.0
  461 do 462 i=1,kkk
  462 pnl(i,m)=pnlo(i)
      go to 487
  470 alfm(m)=alffix
      go to 485
  480 if (delp(m).eq.0.0) go to 484
      b0=eras/delp(m)
      if (b0.gt.-0.20) go to 482
      alfm(m)=max(alfmin,alfm(m)*0.66667)
      go to 485
  482 if (b0.lt.0.20) go to 484
      alfm(m)=min(alfmax,alfm(m)*1.5)
  484 if (alfm(m).ge.1.0) go to 461
  485 a0=1.0-alfm(m)
      do 486 i=1,kkk
  486 xj(i)=alfm(m)*pnlo(i)+a0*pnl(i,m)
      if (ee(m).gt.0.0) go to 487
      call quad5(xj,2,1,kkk,pnorm)
      pnorm=1.0/sqrt(pnorm)
      do 1486 i=1,kkk
 1486 pnl(i,m)=pnorm*xj(i)
  487 if (kkk.ge.mesh) go to 490
      k1=kkk+1
      do 489 i=k1,mesh
  489 pnl(i,m)=0.0
  490 delp(m)=eras
c
  491 if (m.gt.ncores) go to 493
      if (ee(m).gt.0.0) go to 495
      do 492 i=1,kkk
  492 rscore(i)=rscore(i)+wwnl(m)*pnl(i,m)**2
      go to 495
  493 do 494 i=1,kkk
  494 rsvale(i)=rsvale(i)+wwnl(m)*pnl(i,m)**2
c 495 if (npr.le.0) go to 499
  495 go to 499
c              This diagnostic removed 22 Aug 1999
c     do 497 i=1,400
c 497 xj(i)=v(i)*r(i)
c     write (9,49) niter,nlbcd(m),rsint1,e,(xi(i),i=1,400,20), (xj(i),
c    1  i=1,400,20), (pnlo(i), i=1,400,20),(pnl(i,m),i=1,400,20)
   49 format (/1h0,i2,3x,a3,f15.7,f15.6,3x,8hr*vself,,2x,9hr*vtotal,,2x,
     1  7hr*veff,,3x,59h(if diagnostic print out-r*vs, r*vt, calc. wvfn,
     2 mod. wvfn)  / (10f12.6/10f12.6 // ))
  499 continue
      ncspvo=ncspvs
      if (end.gt.0.0) go to 800
c
c     calc total charge density and exchange and correlation potentials
c
  500 do 509 i=1,mesh
      if (ee(ncspvs).le.0.0) go to 502
      rsatom(i)=rscore(i)
      go to 503
  502 rsatom(i)=rscore(i)+rsvale(i)
  503 eras1=r(i)
      eras=rsatom(i)
      call subcor
      rc(i)=b0
      recorr(i)=b1
      rucorr(i)=b2
      ruexch(i)=b3
  509 continue
      rc(1)=rc(2)+rc(2)-rc(3)
c
c       calc ruee=r*(electron-electron energy), xj=r*(total energy)
c
      call quad2(2)
c
c
c          calculate delta
c
  530 delta=0.0
      delto=delt
      idelto=idelta
      do 540 i=2,mesh
      a0=xj(i)-ru(i)
      eras=abs(a0)
      if (eras.lt.delta) go to 540
  535 delta=eras
      delt=a0
      idelta=i
  540 continue
      deltn(niter)=delt
      if (ion.gt.1) go to 545
      if (iz.lt.55.or.iz.gt.102) go to 545
      if (iz.gt.70.and.iz.lt.87) go to 545
      if (niter.gt.2) then
           temp=delt/delto
           if (temp.lt.-0.3) alfmax=0.8*alfmax
      end if
  541 if (niter.lt.8) go to 543
      do 542 n=niter-2,niter-1
  542 if (deltn(n)*deltn(niter).lt.0.0) go to 543
      if (abs(deltn(niter)).gt.0.3*abs(deltn(niter-1)))
     1  alfmax=min(alfmax0,1.2*alfmax)
  543 alfmin=0.2*alfmax
      go to 548
  545 if(niterp.gt.2.and.abs(delt).gt.abs(delto).and.alfmax.gt.0.5) then
           half=0.5
           alfmax=max(0.8*alfmax,half)
           alfmin=min(0.2*alfmax,alfmin)
      end if
  548 call seconds(t2)
      time=t2-t1
      if (iw6.lt.0) write (6,55) niter,time,delt,idelta,kut,niterp,
     1  (alfm(m), m=1,ncspvs)
c            New npr diagnostic, 22 Aug 1999
      if (npr.ne.0.and.niter.ge.npr) call power
  550 if (itpow.ne.1.and.itpow.lt.3) go to 600
      write (9,55) niter,time,delt,idelta,kut,niterp,
     1  (alfm(m), m=1,ncspvs)
   55 format (1x,i2,f6.2,f10.6,i4,i3,i4,1x,15f6.3/f39.3,15f6.3)
c
c          calculate next trial potential.
c             (stabilize ru if niterp=1)
c
  600 a1=1.0
      if (niterp.gt.1) go to 650
      eras=delt/delto
      if (eras.gt.-0.3) go to 603
      alffix=alffix*0.8
      go to 610
  603 if (eras.lt.0.3) go to 610
      if (iabs(idelta-idelto).gt.50) go to 610
      alffix=alffix*1.25
  610 a1=alffix
  650 a0=1.0-a1
      do 660 i=1,mesh
  660 ru(i)=a1*xj(i)+a0*ru(i)
c
c     test self-consistency of atomic potential,
c          and (if thresh=thres1) of the wavefunctions.
c
  700 if (niter.le.1) go to 402
      if (delta.ge.tolend) go to 710
      if (thresh.gt.thres1) go to 710
      if (ihf.eq.1.and.ncspvs.gt.1) go to 780
      if (tolkm2.ge.0.1.and.kut.gt.-1) go to 710
      do 708 m=n1sc,ncspvs
      kkk=nkkk(m)
      if (abs(delp(m))/(abs(pnl(kkk-47,m))+1.0e-07).gt.5.0e-04) goto 710
  708 continue
      if (ee(ncspvs).gt.0.0.and.abs(delp(ncspvs)).gt.1.0e-7) go to 710
      go to 780
  710 if (niter.lt.maxit) go to 730
      if (maxit.eq.maxit0) go to 725
      write (9,71) maxit0,(i,deltn(i),i=niter-2,niter),tolend
   71 format (//' ***************************************************'/
     1          ' *     POOR SCF CONVERGENCE:                       *'/
     2   ' *      maxit0=',i3,'   deltn(',i3,')=',f11.8,'        *'/
     3   ' *                   deltn(',i3,')=',f11.8,'        *'/
     4   ' *                   deltn(',i3,')=',f11.8,'        *'/
     5   ' *                       tolend=',f11.8,'        *'/
     6          ' ***************************************************'/)
      if (iw6.lt.0) write (6,71) maxit0,(i,deltn(i),i=niter-2,niter),
     1  tolend
      go to 780
  725 if (delta.le.10.0*tolend) go to 780
      maxit=maxit+10
      if (delta.lt.100.0*tolend) maxit=maxit+20
  730 if (kut.lt.0) go to 740
      if (kut.ge.0.and.tolkm2.gt.tolend.and.delta.lt.tolkm2) go to 735
      if (delta.lt.tolend) go to 740
      if (delta.ge.tolkm2) go to 740
  735 kut=-1
      exf=exfm1
      if (kutd.ne.-1) go to 400
      kut=-2
      exf=exf10
      go to 400
  740 if (niterp.ge.2) go to 401
      if (delta-tolstb) 401,402,402
c
  780 end=1.0
c     if (npr.eq.9999) npr=npr0
      write (4) ru,ruee,ruexch,rsatom,recorr,rucorr
      if (ihf.eq.0) go to 411
      if (tolend.gt.1.0e-4) go to 800
      if (ihf.eq.1.and.irel.ne.0.and.ncspvs.gt.1) go to 800
      go to 411
c
  800 if (itpow.ne.2) go to 810
      write (9,32)tolstb,tolkm2,tolend,thresh,kut, exf,corrf,ca0,ihf1,
     1  mesh,idb,rdb,emx,r(mesh),irel,ib
      if (ihf1.ne.2) go to 806
      if (nhftrm.gt.0) go to 803
      write (9,33) chxid
      go to 806
  803 write (9,34) t1hf(ii),t2hf(ii),(cfg(j,ii),fg(j,ii),kfg(j,ii),
     1  ifg(j,ii),jfg(j,ii), j=1,nfg)
  806 if (n1sc.gt.1) write (9,36)
      if (n1scf.le.1) go to 810
      write (9,27) n1scfm1
c
  810 if (ihf.eq.0) call clock
      call outpt
      if (nnn(ncspvs).gt.98) call phshift
      if (ihf.ne.0.or.nhftrm.eq.0) go to 840
      if (ii.eq.1) write (7,55)
      do 820 i=2,npar
  820 vprm(i)=vprm(i)/109.73731
      write (7,82) conf,t1hf(ii),t2hf(ii),(vprm(i),kpar(i),i=1,npar)
   82 format (2a6,a2,a7,a3,f12.4,i1,8(f9.5,i1)/28x,9(f9.5,i1))
c
  840 if (ivinti.gt.0) call vinti(ivinti)
      call clock
      if (iw6.ge.0) go to 845
      tr=t2-t1
      tj=t2-t0
      write (6,85) conf,tr,tj
   85 format (' finished  ',2a6,a2,'  at   t(run)=',f6.2,',',5x,
     1  't(job)=',f7.2,' seconds')
  845 t1=t2
  850 continue
c
  900 if (tolkm2.ne.tolend.and.kut.gt.-1) go to 210
  991 kut=kut-1
      if (kut.eq.kutd) go to 991
      if (kut.lt.-2) go to 210
      go to 300
c
      end
      subroutine analyz
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70,chtemp*70,lbcd(21)*1
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c5/izo,iono,idum(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      character*6 wwwnl(25),wwnl8(8)
      common/charww/wwwnl,wwnl8
c
      data lbcd/'s','p','d','f','g','h','i','k','l','m',
     1          'n','o','q','r','t','u','v','w','x','y','z'/
c
      i=6
  110 nbl=0
  112 i=i+1
      if (inpcard(i:i).ne.' ') go to 110
      nbl=nbl+1
      if (nbl.ge.3.or.i.ge.18) go to 120
      go to 112
c          configuration identification
  120 chtemp=' '
      i=min0(i,18)
      chtemp(1:i)=inpcard(1:i)
      read (chtemp,12) conf
   12 format (3a6)
c
      i3=i
      do 125 j=1,8
      nlbcd8(j)='   '
  125 wwnl8(j)='  ' 
      kut=0
      ee8=0.0
      nalf=0
      alfmin=0.0
      alfmax=0.0
      norb=0
c
c
  130 i=i3
      i1=0
      i2=0
      i3=0
      idec=0
      ineg=0
      chtemp=' '
c
  132 i=i+1
      if (i.gt.70) go to 900
      if (inpcard(i:i).ne.' ') go to 134
      if (i3.ne.0) go to 150
      go to 132
  134 if (i1.eq.0) i1=i
      if (i2.gt.0) go to 138
      do 136 j=1,20
      if (inpcard(i:i).eq.lbcd(j)) i2=i
  136 continue
      if (inpcard(i:i).eq.'.') idec=i
      if (inpcard(i:i).eq.'-') ineg=i
  138 i3=i
      if (i.lt.70.or.i1.eq.0) go to 132
c
  150 chtemp='     '
      if (i2.ne.0.and.i1.ge.i2-2) go to 300
      if (norb.eq.0.and.nalf.eq.0.and.i3.gt.i1.and.ineg.eq.0) go to 200
      if (idec.gt.0) go to 500
      if (ineg.gt.0.or.i1.eq.i3) go to 400
      go to 500
c          alfmax and (if i3-i1.gt.1) alfmin
  200 if (i2.gt.0) i3=i2-3
      chtemp(i1-i3+4:4)=inpcard(i1:i3)
      read (chtemp,20) ialfmin,ialfmax
   20 format (2i2)
      alfmin=0.01*ialfmin
      alfmax=0.01*ialfmax
      nalf=1
      go to 130
c          orbital
  300 i3=min0(i3+1,i2+6)
      chtemp(3+i1-i2:9)=inpcard(i1:i3)
      norb=norb+1
      read (chtemp,30) nlbcd8(norb),wwnl8(norb)
   30 format (a3,a6)
      go to 130
c          kut
  400 chtemp(1:2)=inpcard(i3-1:i3)
      read (chtemp,40) kut
   40 format (i2)
      go to 130
c          ee8
  500 chtemp=inpcard(i1:i3)
      read (chtemp,50) ee8
   50 format (f70.3)
      go to 130
c
  900 continue
c     write (6,90) inpcard
c  90 format (/3x,a70)
c     write (6,91) conf,alfmin,alfmax,(nlbcd8(m),wwnl8(m),m=1,8),
c    1  kut,ee8
c  91 format (1x,3a6,2f3.1,8(a3,a6),i2,f6.2)
      return
      end
      block data
 
      implicit real*8 (a-h,o-z)
      character*6 wwwnl(25),wwnl8(8)
      common/charww/wwwnl,wwnl8
      parameter (ko=20)
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
c
      data iblank /2h  /
      data (morb(i),i=1,7)   /1,3,5,8,11,15,19/
      data (lbcd(i),i=1,21)/1hs,1hp,1hd,1hf,1hg,1hh,1hi,1hk,1hl,1hm,
     1                      1hn,1ho,1hq,1hr,1ht,1hu,1hv,1hw,1hx,1hy,1hz/
      data (wwwnl(m),m=1,25) /'0 ','1 ','2 ','3 ','4 ','5 ','6 ','7 ',
     1           '8 ','9 ','10','11','12','13','14','00','01','02','03',
     2           '04','05','06','07','08','09'/
c
      end
      subroutine setcfg
c
c     set up heaviest noble gas configuration with no. of electrons
c          .le.iz-nspect+1, and modify as indicated on config. input car
c     estimate eigenvalue for each orbital
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
      dimension ne(7)
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
      character*6 wwwnl(25),wwnl8(8)
      common/charww/wwwnl,wwnl8
c
      data (ne(i),i=1,7)     /2,10,18,36,54,86,118/
c
c     set up noble gas core with atomic number less than or equal to noe
c
c     first find number of orbitals in this core from table
c
      ierror=0
  218 mcore=0
      necore=0
      ifrac=0
      do 220 n=1,7
      if (ne(n).gt.noelec) go to 221
      necore=ne(n)
  220 mcore=morb(n)
  221 m=0
      do 225 n=1,7
      do 225 lp1=1,n
      if (mcore.eq.m) go to 230
      m=m+1
      nnn(m)=n
      l(m)=lp1-1
      i=l(m)+1
c     encode (3,22,nlbcd(m))  nnn(m),lbcd(i)
      write (nlbcd(m),22)  nnn(m),lbcd(i)
   22 format (i2,a1)
      wwnl(m)=4*lp1-2
      if (necore.eq.54.and.n.eq.4.and.l(m).eq.3)  m=m-1
      if (necore.eq.86.and.n.eq.5.and.l(m).eq.3)  m=m-1
      if (necore.eq.86.and.n.eq.5.and.l(m).eq.4)  m=m-1
      if (necore.eq.118.and.n.eq.5.and.l(m).eq.4)  m=m-1
      if (necore.eq.118.and.n.eq.6.and.l(m).eq.3)  m=m-1
      if (necore.eq.118.and.n.eq.6.and.l(m).eq.4)  m=m-1
      if (necore.eq.118.and.n.eq.6.and.l(m).eq.5)  m=m-1
  225 continue
c
c     decode and store n, l, and occupation nos. for subshells
c          outside of or which modify the core
c
  230 do 241 i=1,8
      if (nlbcd8(i)(1:3).eq.'   ') go to 241
      do 233 j=1,25
      if (wwnl8(i).ne.wwwnl(j)) go to 233
      ww8=j-1
      if (j.lt.16) go to 234
      ww8=j-16
      go to 234
  233 continue
      ww8=1.0
      if (ww8.eq.1.0.and.(wwnl8(i)(3:3).ne.' '.or.wwnl8(i)(2:2).eq.'.'
     1  .or.wwnl8(i)(1:1).eq.'.')) read (wwnl8(i),1023) ww8
 1023 format (f6.3)
      if (mod(ww8,1.0).ne.0.0) ifrac=1
  234 if (m.eq.0) go to 237
      m1=m
      do 236 j=1,m1
      if (nlbcd8(i).ne.nlbcd(j)) go to 236
      if (ww8.ne.0.0.or.j.lt.n1scf) go to 240
      m=m-1
      if (j.eq.m+1) go to 241
      do 235 mp=j,m
      nnn(mp)=nnn(mp+1)
      l(mp)=l(mp+1)
      nlbcd(mp)=nlbcd(mp+1)
  235 wwnl(mp)=wwnl(mp+1)
      go to 241
  236 continue
  237 if (ww8.eq.0.0.and.m+1.ge.n1scf) go to 241
      m=m+1
      j=m
      nlbcd(m)=nlbcd8(i)
c     decode (3,22,nlbcd(m))  nnn(m),li
      read (nlbcd(m),22)  nnn(m),li
      do 238 iii=1,21
      if (li.eq.lbcd(iii)) go to 239
  238 continue
      write (9,23)  m,nlbcd(m)
   23 format (7h0nlbcd(,i2, 4h) = ,a3,
     1  43h was not found in table.  next case please.  ///)
      go to 310
  239 l(m)=iii-1
      if (iii.gt.nnn(m)) go to 1241
  240 wwnl(j)=ww8
  241 continue
      ncspvs=m
      ncores=m-1
      if (m.gt.0) go to 242
 1241 write (9,24) iz,noelec,conf,alfmin,alfmax, (nlbcd8(m),wwnl8(m),
     1  m=1,8), ee8
   24 format (55h0ncspvs = 0  you have troubles with the following case-
     1  //2(2x,i3),3a6,2f5.2,8(a3,a6),0pf8.5//22h try the next element
     2  //1h1  )
      go to 310
c
c     accumulate no. of electrons, estimate eigenvalue for each orbital
c
  242 www=0.
      do 243 m=1,ncspvs
      www=www+wwnl(m)
      if (m.eq.ncores) ncelec=www
      eras=nnn(m)-l(m)
  243 a0m(m)=ca0/eras
      if (ee8.gt.0.0) www=ncelec
      noelec=www
      nvelec=noelec-ncelec
      ion=iz-noelec
      twoion=ion+ion
      zzz=ion+1
      twozzz=zzz+zzz
      n1sc=1
      if (eeo(ncspvs).lt.0.0.or.ee(ncspvs).lt.0.0.or.ee8.le.0.0)
     1  go to 265
      if (iz.ne.izo.or.ion.ne.iono+1) go to 265
      do 262 m=1,ncspvs
      if (m.eq.ncspvs) go to 262
      if (nlbcd(m).ne.nlbcdo(m).or.wwnl(m).ne.wwnlo(m)) go to 265
  262 continue
      n1sc=ncspvs
      kut=-1
  265 do 250 m=1,ncspvs
      if (n1sc.gt.1) go to 246
      smsig=0.
      if (iz.ne.izo) go to 244
      if (ion.ne.iono.and.m.ge.ncspvs-2) go to 244
      if (m.ge.ncspvs) go to 244
      if (nlbcd(m).eq.nlbcdo(m)) go to 246
  244 eras=l(m)
      sig1=0.6+0.05*eras
      do 245  mp=1,ncspvs
      eras=nnn(m)-nnn(mp)
      sig=sig1+0.2*eras
      if (sig.lt.0.0) go to 245
      if (sig.gt.1.0) sig=1.0
      wwnlx=wwnl(mp)
      if (nnn(m).eq.nnn(mp).and.l(m).eq.l(mp)) wwnlx=wwnlx-1.0
      smsig=smsig+sig*wwnlx
  245 continue
      if (nnn(m)+l(m).eq.1.and.wwnl(m).eq.2.0) smsig=0.5*smsig+0.3
      ee(m)=-((z-smsig)/float(nnn(m)))**2
  246 if (n1sc+n1scf.gt.2.and.m.lt.ncspvs) ee(m)=eeo(m)
      if (m.eq.ncspvs.and.ee8.gt.0.0) ee(m)=ee8
      if (itpow.lt.8.or.n1sc.gt.1) go to 250
      write (9,25) m,nnn(m),l(m),wwnl(m),nlbcd(m),wwnl(m),ee(m), smsig
   25 format (4h m =,i2, i3,i2,f7.3,3x,a3,f6.3,f14.5, 9h rydbergs, 5x,
     1 38hscreening of electron in mth orbital =  ,f9.5  )
  250 continue
      if (itpow.lt.8) go to 260
      write (9,26) mcore, necore, www
   26 format (8h0mcore =,i3, 3x,8hnecore =,i3, 8h   www =,f3.0  //  )
  260 do 268 m=1,ncspvs
      wwnlo(m)=wwnl(m)
  268 nlbcdo(m)=nlbcd(m)
      m=ncspvs
      if (ee(m).ne.0.0) go to 300
      fn=nnn(m)
      ee(m)=-(zzz/fn)**2
  300 return
  310 ierror=7
      return
      end
      subroutine diel
c
c     subroutine to construct continuum configuration cards for
c          rcn/rcg autoionization-problem runs
c
c          to use this option, make tape10=in36 an empty file, 
c               or a file containing a single record with "-1" in columns 4-5;
c
c          input (in file called tape12=indiel) consists of:
c     (1) rcn36 control card (with appropriate value of emx, or zero)
c     (2) ion config. card(s) to which a continuum electron is to be added
c          (1 to 5 in number)
c     (3) first-parity (radiative-decay) configurations
c     (4) second-parity (doubly excited) configuration(s)
c          (1 to 10 in number)
c     (5) card with -1 in columns 4-5 (rcn exit card)
c          further sets of cards (2)-(5) may be added as desired
c
c          action of subroutine:
c     (a) first portion of program writes control card (1) on
c          tape10 (input file for preliminary rcn run) and on
c          tape11=in34 (input file for final rcn/rcg run)
c     (b) statements 100-300 copy cards (2) and (4)-(5) on tape 10
c     (c) rcn is called to provide data for calculation of
c          kinetic energy of free electron
c     (d) statements 410-799 copy cards (3)-(4) onto tape11 and
c          add continuum configuration cards
c     (e) parts (b)-(d) are repeated for each additional
c          set of input data cards (2)-(5)
c     (f) statement 800 adds the rcn exit card to tape11
c
c          final tape11 (copied to tape10) is ready to be used directly as
c            input file for final rcn/rcn2/rcg run, except for the following
c          (i) it may be necessary to revise the value of emx on the
c            control card if it was not adequately estimated on input card (1)
c          (ii) if the parities of cards (3)--or the parities of cards (4)--
c               were not the same in contiguous sets but the values of z
c               were the same, then it is necessary to break tape11
c               into two pieces, one for each parity, because of the requirements
c               of rcn2
c               (this can be done automatically by including at the appropriate
c               place on tape12 a single ion card of different z,
c               plus an rcn exit card (5). )
c
c
      implicit real*8 (a-h,o-z)
      character*8 ia(10),ib(3),inpcard*70,iba,icard*80,form*42,hibo*5
      character*2 ioccch,iwwnl(26),iochion(8,10)
      character*1 lbcd(21),lh,lhion(8,10)
      common/ch/ia,ib,inpcard,iba
      parameter (ko=20)
      common/d0/nnn(ko),lll(ko),wwnl(ko),ncspvs,vprm
      dimension nion(8,10),
     1  lion(8,10),ioccion(8,10),mionk(10),eavion(10),alf(10)
      common/d3/natom(8,10),latom(8,10),ioccatm(8,10),matmk(10),la(2)
c
      data lbcd/'s','p','d','f','g','h','i','k','l','m',
     1          'n','o','q','r','t','u','v','w','x','y','z'/
      data (iwwnl(m),m=1,26) /'0 ','1 ','2 ','3 ','4 ','5 ','6 ','7 ',
     1         '8 ','9 ','10','11','12','13','14','00','01','02','03',
     2         '04','05','06','07','08','09','  '/
      data alf/1ha,1hb,1hc,1hd,1he,1hf,1hg,1hh,1hi,1hj/
      data fff1,fff2,fff3/5hf6.4,,5hf6.2,,5hf6.1,/
c
      rewind 12
      rewind 11
      rewind 10
c          read rcn control card
      read (12,10) ia
   10 format (10a8)
      write (11,10) ia
      write (10,10) ia
      write (9,11) ia
   11 format (//1h ,10a8)
c          read ion configuration card(s)
  100 rewind 10
      read (10,10) ia
      ird=0
      iparty0=-7
      i=0
  150 i=i+1
      read (12,22,end=800) iz,ispect,inpcard
      ird=ird+1
      if (i.eq.1) isp1=ispect
      if (i.eq.1) izo=iz
      if (i.gt.2.or.iz.gt.0) go to 160
      ispect=izo-1
      hibo='dummy'
      write (11,22) izo,ispect,hibo
      go to 100
  160 if (isp1.ne.ispect) go to 230
      write (10,22) iz,ispect,inpcard
      write (9,30) iz,ispect,inpcard
      go to 150
c          select second-parity configuration cards
  220 read (12,22) iz,ispect,inpcard
   22 format (2i5,a70)
      ird=ird+1
  230 if (iz.le.0) go to 300
  240 iparity=0
      call analyz1
      do 280 i=1,8
c     decode (5,24,ia(i)) n,lh,ioccch
      read (ia(i),24) n,lh,ioccch
   24 format (i2,a1,a2)
      if (n.le.0) go to 290
      do 245 j=1,21
      if (lh.eq.lbcd(j)) go to 250
  245 continue
  250 l=j-1
      do 255 m=1,26
      if (ioccch.eq.iwwnl(m)) go to 260
  255 continue
      m=1
  260 iocc=m-1
      if (m.gt.15) iocc=m-16
      if (m.eq.26) iocc=1
      iparity=iparity+l*iocc
  280 continue
  290 iparity=mod(iparity,2)
      if (iparty0.eq.-7) iparty0=iparity
      if (iparity.eq.iparty0) go to 220
c
  300 if (iz.lt.0) iz=-999
      write (10,22) iz,ispect,inpcard
      write (9,30) iz,ispect,inpcard
   30 format (1h ,2i5,a70)
      if (iz.gt.0) go to 220
c
      rewind 10
c          run rcn to calculate free-electron energies
      call rcn
      rewind 2
c
      do 410 i=1,ird
  410 backspace 12
c          read ion configuration(s)
      k=0
  420 k=k+1
      read (12,52) iz,ispion,inpcard
      if (k.eq.1) isp=ispion
      if (ispion.eq.isp) go to 430
      kion=k-1
      backspace 12
      go to 500
  430 iparity=0
      call analyz1
      do 480 i=1,8
c     decode (5,54,ia(i)) nion(i,k),lhion(i,k),iochion(i,k)
      read (ia(i),54) nion(i,k),lhion(i,k),iochion(i,k)
      if (nion(i,k).gt.0) mionk(k)=i
      do 435 j=1,21
      if (lhion(i,k).eq.lbcd(j)) go to 450
  435 continue
      j=1
  450 lion(i,k)=j-1
      do 455 m=1,26
      if (iochion(i,k).eq.iwwnl(m)) go to 460
  455 continue
      m=1
  460 iocc=m-1
      if (m.gt.15) iocc=m-16
      if (m.eq.26) iocc=1
      if (nion(i,k).eq.0) iocc=0
      iparity=iparity+lion(i,k)*iocc
      ioccion(i,k)=iocc
  480 continue
      i=1
  490 read (2) nnn,lll,wwnl,ncspvs,vprm
      eavion(k)=vprm
      go to 420
c          copy bound-state configurations
  500 iparty0=-7
      k=1
  520 read (12,52) iz,ispect,inpcard
   52 format (2i5,a70)
      if (iz.le.0) go to 600
      write (11,52) iz,ispect,inpcard
      write (9,52) iz,ispect,inpcard
  540 iparity=0
      call analyz1
      do 580 i=1,8
c     decode (5,54,ia(i)) n,lh,ioccch
      read (ia(i),54) n,lh,ioccch
   54 format (i2,a1,a2)
      if (n.le.0) go to 590
      natom(i,k)=n
      do 545 j=1,21
      if (lh.eq.lbcd(j)) go to 550
  545 continue
      j=1
  550 latom(i,k)=j-1
      do 555 m=1,26
      if (ioccch.eq.iwwnl(m)) go to 560
  555 continue
      m=1
  560 iocc=m-1
      if (m.gt.15) iocc=m-16
      if (m.eq.26) iocc=1
      iparity=iparity+latom(i,k)*iocc
      ioccatm(i,k)=iocc
  580 continue
  590 iparity=mod(iparity,2)
      if (iparty0.eq.-7) iparty0=iparity
      katm=k
      matmk(k)=i-1
      iza=iz
      ispa=ispect
      iba=ib(1)
      if (iparity.ne.iparty0) k=k+1
      go to 520
c          set up continuum cards
  600 do 799 ka=1,katm
      matom=matmk(ka)
      read (2) nnn,lll,wwnl,ncspvs,vprm
      do 798 k=1,kion
      mion=mionk(k)
      eps=vprm-eavion(k)
      write (6,60) ka,k,eps
      if (eps.gt.0.0) go to 605
      write (9,60) ka,k,eps
   60 format (/' for atom config.',i2,' and ion config.',
     1  i2,',   eps=',f10.5)
      go to 798
  605 ji=0
      do 649 mi=1,mion
      do 609 ma=1,matom
      if (nion(mi,k).ne.natom(ma,ka).or.lion(mi,k).ne.latom(ma,ka))
     1  go to 609
      if (ioccion(mi,k).le.ioccatm(ma,ka)) go to 649
  609 continue
  610 ji=ji+1
      if (ji.eq.1) go to 620
      j=1
      ma=8
      write (9,61) j,mi,ma
   61 format (11h0*****error,i1,2i5/)
      go to 798
  620 li1=lion(mi,k)
  649 continue
      ja=0
      do 699 ma=1,matom
      iocca=ioccatm(ma,ka)
      do 659 mi=1,mion
      if (nion(mi,k).ne.natom(ma,ka).or.lion(mi,k).ne.latom(ma,ka))
     1  go to 659
      if (ioccion(mi,k).ge.iocca) go to 699
      go to 680
  659 continue
      mi=8
  680 ja=ja+1
      if (ja.le.2) go to 685
      j=2
      mi=8
      write (9,61) j,mi,ma
      go to 798
  685 la(ja)=latom(ma,ka)
      iocca=iocca-1
      if (iocca.gt.ioccion(mi,k)) go to 680
  699 continue
c
  700 k1=iabs(la(1)-la(2))
      k2=la(1)+la(2)
      lcontx=li1+k2
      lcontn=mod(lcontx,2)
  702 if (iabs(li1-lcontn).le.k2.and.li1+lcontn.ge.k1) go to 703
      lcontn=lcontn+2
      go to 702
  703 continue
      fff=fff2
      if (eps.gt.999.0) fff=fff3
      if (eps.lt.9.99) fff=fff1
      mp1=mion+1
c     encode (42,71,form) fff,mp1,fff
      write (form,71) fff,mp1,fff
   71 format (11h(2i5,a5,a1,,a5,6ha1,9x,,i1,11h(i2,a1,a2),,3h3x,,a4,1h))
      kk=99
      if (lcontx.gt.9) lcontx=9
      if (lcontx.lt.lcontn) go to 798
      do 720 l=lcontn,lcontx,2
      lh=lbcd(l+1)
      write (11,form) iza,ispa,iba,alf(k),eps,lh, (nion(i,k),
     1  lhion(i,k),iochion(i,k), i=1,mion), kk,lh,iwwnl(26),eps
  720 write (9,form) iza,ispa,iba,alf(k),eps,lh, (nion(i,k),
     1  lhion(i,k),iochion(i,k), i=1,mion), kk,lh,iwwnl(26),eps
  798 continue
  799 continue
      go to 100
c
  800 write (11,80)
   80 format (3x,2h-1,10x)
      rewind 11
      rewind 10
  850 read (11,81,end=900) icard
   81 format (a)
      write (10,82) icard
   82 format (a80)
      go to 850
  900 rewind 10
      rewind 9
      return
      end
      subroutine analyz1
c
      implicit real*8 (a-h,o-z)
      character*8 ia(10),ib(3),inpcard*70,iba
      common/ch/ia,ib,inpcard,iba
      character chtemp*70,lbcd(21)*1
c
      data lbcd/'s','p','d','f','g','h','i','k','l','m',
     1          'n','o','q','r','t','u','v','w','x','y','z'/
c
      i=6
  110 nbl=0
  112 i=i+1
      if (inpcard(i:i).ne.' ') go to 110
      nbl=nbl+1
      if (nbl.ge.3.or.i.ge.18) go to 120
      go to 112
c          configuration identification
  120 chtemp=' '
      i=min0(i,18)
      chtemp(1:i)=inpcard(1:i)
      read (chtemp,12) ib
   12 format (3a8)
c
      i3=i
      kut=0
      ee8=0.0
      nalf=0
      alfmin=0.0
      alfmax=0.0
      norb=0
      do 125 j=1,10
  125 ia(j)='        '
c
c
  130 i=i3
      i1=0
      i2=0
      i3=0
      idec=0
      ineg=0
      chtemp=' '
c
  132 i=i+1
      if (i.gt.70) go to 900
      if (inpcard(i:i).ne.' ') go to 134
      if (i3.ne.0) go to 150
      go to 132
  134 if (i1.eq.0) i1=i
      if (i2.gt.0) go to 138
      do 136 j=1,20
      if (inpcard(i:i).eq.lbcd(j)) i2=i
  136 continue
      if (inpcard(i:i).eq.'.') idec=i
      if (inpcard(i:i).eq.'-') ineg=i
  138 i3=i
      if (i.lt.70.or.i1.eq.0) go to 132
c
  150 chtemp='     '
      if (i2.ne.0.and.i1.ge.i2-2) go to 300
      if (norb.eq.0.and.nalf.eq.0.and.i3.gt.i1.and.ineg.eq.0) go to 200
      if (idec.gt.0) go to 500
      if (ineg.gt.0.or.i1.eq.i3) go to 400
      go to 500
c          alfmax and (if i3-i1.gt.1) alfmin
  200 if (i2.gt.0) i3=i2-3
      chtemp(i1-i3+4:4)=inpcard(i1:i3)
      read (chtemp,20) ialfmin,ialfmax
   20 format (2i2)
      alfmin=0.01*ialfmin
      alfmax=0.01*ialfmax
      nalf=1
      go to 130
c          orbital
  300 i3=min0(i3+1,i2+6)
      chtemp(3+i1-i2:9)=inpcard(i1:i3)
      norb=norb+1
      read (chtemp,30) ia(norb)
   30 format (a8)
      go to 130
c          kut
  400 chtemp(1:2)=inpcard(i3-1:i3)
      read (chtemp,40) kut
   40 format (i2)
      go to 130
c          ee8
  500 chtemp=inpcard(i1:i3)
      read (chtemp,50) ee8
   50 format (f70.3)
      go to 130
c
  900 continue
      write (9,90) inpcard
   90 format (/3x,a70)
      write (9,91) ib,alfmin,alfmax,(ia(m),m=1,8),kut,ee8
   91 format (1x,2a8,a2,2f3.1,8a9,i2,f6.2)
      return
      end
      subroutine scheq
c
c     compute energy eigenvalue and wave function
c      originally  written by sherwood skillman
c      rca laboratories, princeton, new jersey, spring 1961
c      modified by frank herman, summer 1961
c      further modified by richard kortum,  lockheed research
c       laboratories, palo alto, california,  summer 1962
c      further modified by r. d. cowan, lasl, jan-feb, 1964
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
      common/c9/ qq(kmsh),delp(ko),alfm(ko)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      dimension p(5),d(5),q(5),t(5)
      save
c
c
c
c     set up constants and initialize
c
      erasx=0.0
  100 many=99
      emore=-1.0e-10
      eless=-1.0e+12
      eg=0.0
      de=0.0
      nprint=0
      imatch=0
      lamp=lam+1
      xlp=lamp
      ndcr=nn-lamp
      if (nn.gt.98) ndcr=9999
      b0=lam*lamp
      h=r(2)
      h1=h*h
      b3=(v(3)-v(2))/h-z/h1
      y=h+h
      flps=4*lam+6
      slpt=6*lam+12
      elpt=8*lam+20
      a1=-z/xlp
      b1=-z-z
      ab1=a1*b1
      ab3=a1*b3
c
c     raise h and y to lam+1
c
      htl=h
      ytl=y
      if (lam) 210,140,120
  120 do 130 i=1,lam
      htl=htl*h
  130 ytl=ytl*y
  140 bohs=b0/h1
      boh=b1/h
      bth=b3*h
      bq3=bohs+boh+bth
      bq4=bohs/4.0+boh/2.0+bth+bth
  150 if ((niter.eq.1.and.n1sc.eq.1).or.irel.eq.0) go to 200
      if (lam.ne.0) go to 200
c          calc darwin relativistic term
      kkk=nkkk(mm)
      nblk=kkk/40
      i=2
      k1=2
      k2=40
      do 160 j=1,nblk
      if (j.eq.nblk) k2=39
      eras=0.5/(r(i+1)-r(i))
      do 155 k=k1,k2
      i=i+1
      eras1=eras*(pnl(i+1,mm)-pnl(i-1,mm))
      eras1=eras1/pnl(i,mm)-1.0/r(i)
  155 vr(i)=eras1*eras*(x2(i+1)-x2(i-1))
      k1=1
      if (i.gt.idb.or.j.eq.nblk) go to 160
      eras=0.5*eras
      eras1=eras*(pnl(i+1,mm)-pnl(i-2,mm))
      eras1=eras1/pnl(i,mm)-1.0/r(i)
      vr(i)=eras1*eras*(x2(i+1)-x2(i-2))
  160 continue
      vr(2)=4.0*vr(3)
      vr(mesh)=vr(mesh-1)
      if (iptvu.ge.6.and.mod(niter-2,5).eq.0)
     1  write (9,7788) nn,lam,eras, (vr(k), k=2,47)
 7788 format (2i10/(1p,8e15.5))
      if (irel.le.3) go to 180
c     eras=1.0/(12.0*(r(2)-r(1)))
c     do 165 i=4,39
c     eras1=eras*(8.0*(pnl(i+1,mm)-pnl(i-1,mm))-pnl(i+2,mm)+pnl(i-2,mm))
c     eras1=eras1/pnl(i,mm)-1.0/r(i)
c 165 vr(i)=eras1*eras*(8.0*(x2(i+1)-x2(i-1))-x2(i+2)+x2(i-2))
c     eras1=log(vr(4)/vr(5))/log(r(5)/r(4))
c     vr(3)=vr(4)*((r(4)/r(3))**eras1)
c     vr(2)=vr(4)*((r(4)/r(2))**eras1)
c     vr(mesh)=vr(mesh-1)
      if (iptvu.ge.5.and.mod(niter-2,5).eq.0)
     1  write (9,7788) nn,lam,eras, (vr(k), k=2,47)
  180 do 185 i=20,mesh-2
      if (vr(i+1)*vr(i).gt.0.0) go to 185
      j=i
      if (abs(vr(i+1)).gt.abs(vr(i))) j=i+1
      vr(j)=0.5*(vr(j+1)+vr(j-1))
  185 continue
c
c     start outward integration
c
  200 nprint=nprint+1
      j=iabs(npr)
      if (j.eq.0.or.j.gt.9) go to 880
      if (mod(nprint,j).ne.0) go to 880
c     if (mm.lt.ncspvs-2) go to 880
      write (9,17) niter,nprint,nn,lam,ik,imatch,kkk,ncross,
     1  de,eg,emore,eless
   17 format (//3i3,i1,4i5,4f10.6)
      write (9,18) (qq(j),j=201,kkk,10)
   18 format (1x,10f12.4)
      write (9,19) (pnlo(j),j=201,kkk,10)
   19 format (1x,10f12.7)
  880 continue
      eps =e-eg
      eg =e
      if (nprint-many) 220,220,201
  201 write (9,20) nlbcd(mm),nprint,eless,emore,e, de,ncross,imatch,i,ik,
     1  ikmax, (pnlo(j), j=1,mesh,2)
   20 format (' *****No convergence on ',a3,i4,4f17.8,5i4/
     1  /' EMX may need to be non-zero, but small'//(3x,1p,10e11.3))
      write (6,20) nlbcd(mm),nprint,eless,emore,e, de,ncross,imatch,i,ik,
     1  ikmax
      write (6,1020) r(mesh),emx,nlbcd(mm)
 1020 format (' Possibly r(mesh)=', f10.5,' is too small because',
     1  ' EMX=0.0 (mesh=641 instead of 1801),'/
     2  '         or because EMX is too large (idb too small--increase',
     3  ' idb by 40 by decreasing EMX=',f7.2,' by a factor 4):'/
     4  '     r(mesh) needs to be 5 or 10 times as large as <r> for ',
     5  'P(',a3,').')
      if (abs(de/e)-8.0*thresh) 900,205,205
  205 end=1.0
      go to 990
  210 nstop=210
  211 write (9,21)nstop
   21 format(5h0stop,i4,9h in scheq)
      stop 210
  220 do 221 i=1,mesh
  221 pnlo(i)=0.0
c
c     calc qq.    if e neg, make certain that qq(mesh) is positive,
c         and that qq(i) is negative for some i less than ikmax
c
      n2bg=nprint+nn
      if (nprint-1) 210,230,235
  230 do 231 i=1,mesh
  231 qq(i) = v(i)+b0/(r(i)*r(i))-e
      if ((niter.eq.1.and.n1sc.eq.1).or.irel.eq.0) go to 233
c          calc mass-velocity relativistic term
      eras=(0.5/137.036)**2
      do 232 i=2,kkk
      eras1=e-v(i)
      eras2=eras*eras1
      qq0(i)=qq(i)
      eras3=eras1*eras2
      if (lam.eq.0) eras3=eras3+eras*vr(i)/(1.0+eras2)
  232 qq(i)=qq(i)-eras3
      if (niter.ne.8.or.irel.lt.4) go to 233
      write(9,22) niter,nprint,mm,b0,e
   22 format (1h1,3i5,2f17.7//)
      write(9,23)(i,r(i),v(i),vr(i),qq0(i),qq(i),pnl(i,mm), i=11,kkk,10)
   23 format (i5,6f18.7)
  233 ik=mesh+1
      if (e) 240,240,250
  235 do 236 i=1,mesh
  236 qq(i) = qq(i)-eps
      if (nprint-20) 240,240,237
  237 if (abs(eps/e)-0.00001) 250,250,240
  240 m=mesh
      ik=m+1
      do 246 i=10,mesh
      if (qq(m)) 245,246,246
  245 ik=m+20
      go to 247
  246 m=m-1
      eps=-0.3*e
      go to 248
  247 ikmax=mesh-min0( 60, 6000/(nprint+nn)**2 )
      if (ik.le.ikmax+5) go to 250
      mesh=mesh
      eps=min(qq(ikmax-10),qq(mesh)-0.000001)
  248 e=e+eps
      eg=e
      n2bg=n2bg+1
      if (n2bg.le.60) go to 235
      if (nprint.lt.many) go to 250
      write (9,24) nn,lam, (qq(i), i=201,mesh)
   24 format (//2i3//(5x,10f11.5))
c
c          calc starting values for irel=0
c     b0=lam*(lam+1)
c     b1= -2.*z
c     b3=(v(3)-v(2))/h -z/hsq
c     a1= -z/(lam+1)
c     a2=(a1*b1+b2)/(4*lam+6)
c     a3=(a2*b1+a1*b2+b3)/(6*lam+12)
c     a4=(a3*b1+a2*b2+a1*b3)/(8*lam+20)
c     p(3)=(1.0+a1*h+a2*h**2+a3*h**3+a4*h**4)*h**(xl+1.0)
c     p(4)=(1.0+a1*y+a2*y**2+a3*y**3+a4*y**4)*y**(xl+1.0)
c     q(3)=(b+b1*h+b2*h**2+b3*h**3)/h**2
c     q(4)=(b+b1*y+b2*y**2+b3*y**3)/y**2
c
  250 h=r(2)
      y=h+h
      b2=3.0*z/h-e+2.0*v(2)-v(3)
      a2=(ab1+b2)/flps
      a3=(a2*b1+a1*b2+b3)/slpt
      a4=(a3*b1+a2*b2+ab3)/elpt
      p(3)=(1.0+h*(a1+h*(a2+h*(a3+h*a4))))*htl
      p(4)=(1.0+y*(a1+y*(a2+y*(a3+y*a4))))*ytl
      p(5)=0.0
      q(3)=bq3+b2
      q(4)=bq4+b2
      i=3
      dx=h
      h1=h**2
      h2=h1/12.0
      if (irel.eq.0) go to 255
c          calc starting values for irel.ne.0 (relativistic terms includ
      q(3)=qq(2)
      q(4)=qq(3)
      eras1=(z/137.036)**2
      d2=1.0/(1.0+(274.072**2+e)*r(3)*0.5/z)
      if (lam.ne.0) d2=0.0
      eras=0.5*(1.0-d2)+sqrt(b0+0.25*((1.0+d2)**2)-eras1)
      p(3)=r(2)**eras
      p(4)=r(3)**eras
      eras2=eras+1.0+d2
      a1=-(2.0*z+eras1*e/z)/(eras2*eras-b0+eras1)
      p(3)=p(3)*(1.0+a1*r(2))
      p(4)=p(4)*(1.0+a1*r(3))
c
  255 pnlo(2)=p(3)
      pnlo(3)=p(4)
      t(3)=p(3)*(1.0-h2*qq(2))
      t(4)=p(4)*(1.0-h2*qq(3))
      d(4)=t(4)-t(3)
      ncount=2
      nint=2
      if (nprint-40) 280,260,260
  260 write (9,26) nn,lam,niter,nprint,ikmax,ik,imatch,ncross,e,eps,
     1  eless,emore,alfm(mm),delta
   26 format (1h ,i2,i1,6i6,4f19.11,f6.3,f10.6)
  280 ncross=0
      ib1=ib+1
      ptthree=0.3
c
c          integrate outward
c
  300 i=i+1
      if (i.lt.ib1) go to 310
      if (e.ge.0.0) go to 960
      if (ncross-ndcr) 370,303,360
  303 if (ik-ib1) 370,850,850
c
  310 q(5)=qq(i)
      if (abs(p(5)).gt.1.0e30) go to 312
      if (i.lt.ik) go to 314
  312 if (ncount.gt.1) go to 314
      if (nint.gt.4) go to 350
  314 d(5)=d(4)+h1*q(4)*p(4)
      t(5)=d(5)+t(4)
      eras=h2*q(5)
      if (abs(eras).gt.abs(erasx)) erasx=eras
      if (abs(eras).gt.ptthree) eras=sign(ptthree,eras)
  315 p(5)=t(5)/(1.0-eras)
      pnlo(i)=p(5)
      if (p(5).gt.0.0.and.p(4).gt.0.0) go to 331
      if (p(5).lt.0.0.and.p(4).lt.0.0) go to 331
      if (p(5).eq.0.0) go to 331
c          count changes in sign
  330 ncross=ncross+1
  331 ncount=ncount+1
      if (ncount.eq.6) ncount=1
      nint=nint+1
      if (i.gt.idb) go to 341
      if (nint.lt.40) go to 341
  340 dx=dx+dx
      h=dx
      h1=h**2
      h2=h1/12.0
      nint=0
      t(5)=p(5)*(1.0-h2*q(5))
      t(3)=p(3)*(1.0-h2*q(3))
      d(5)=t(5)-t(3)
  341 do 342 k=1,4
      p(k)=p(k+1)
      t(k)=t(k+1)
      d(k)=d(k+1)
  342 q(k)=q(k+1)
      go to 300
c
c     matching radius has been reached going out
c
  350 if (ncross-ndcr) 370,380,360
c          too many crossings, increase absf(e)
  360 emore=min(emore,e)
      e=max(2.00*e,0.5*(eless+emore))
      go to 200
c          too few crossings, decrease absf(e)
  370 eless=max(eless,e)
      de=emore-eless
      if (de.lt.0.0001) emore=emore+0.1*de
      if (emore.ge.0.0) emore=-1.0e-8
      e=min(0.50*e,0.5*(eless+emore))
      go to 200
c
c     number of crossings is correct
c     check to see that wave is in the damped region (absolute value
c       decreasing and signs alike)
c
  380 if (abs(pnlo(i-1))-abs(pnlo(i-2))) 382,381,381
c     large absolute value of p in what should be the damped region
c       indicates too few peaks, decrease absf(e)
  381 if (abs(p(5))-1.0e05) 314,370,370
c
  382 if (p(5)) 383,314,384
  383 if (pnlo(i-2)) 400,314,314
  384 if (pnlo(i-2)) 314,314,400
c
c     now ndcr = ncross and matching radius lies in damped region
c         calculate logarithmic derivative
c
  400 imatch=i-2
      xmatch=r(i-2)
      ppout=(t(4)-t(2)-0.5*(p(4)-p(2)))/h
      s6=ppout/p(3)
c
c  integrate pnlo**2
c
      call quad5(pnlo,2,1,imatch,sum1)
c
      s5=sum1/p(3)**2
      pmatch=p(3)
c
c          start inward integration
c
c       choose outermost mesh point, kkk
  500 xinw=10.0*xmatch
c     if (nn.gt.6) xinw=6.0*xmatch
c     if (nn.gt.9) xinw=5.0*xmatch
c     if (nn.gt.12) xinw=4.0*xmatch
  510 do 511 i=41,mesh,40
      if (r(i).ge.xinw) go to 512
  511 continue
      kkk=mesh
      go to 550
  512 kkk=i
c       calc starting values
  550 i=kkk
      dx=r(i)-r(i-1)
      h=-dx
      h1=h*h
      h2=h1/12.0
      q(3)=qq(i)
      eras1=sqrt(q(3))
      eras2=r(i)*eras1
      p(3)=exp(-eras2)
      i=i-1
      q(4)=qq(i)
      eras3=r(i)*sqrt(q(4))
      p(4)=exp(-eras3)
      if (abs(p(4)).gt.1.0e-35) go to 580
      kkk=kkk-40
      if (kkk.gt.imatch) go to 550
  570 write (9,57)z,nn,lam,kkk,imatch
   57 format (6h0at z=,f6.0,6h   nl=,i3,i1,10h,     kkk=,i3,
     1    22h  is less than imatch=,i3,
     2    48h.    inward integration will be tried at kkk+40.)
      kkk=kkk+40
      p(3)=1.0e-35
      p(4)=p(3)*exp(eras2-eras3)
  580 if (pmatch.gt.0.0) go to 582
      p(3)=-p(3)
      p(4)=-p(4)
  582 pnlo(i+1)=p(3)
      pnlo(i)=p(4)
      t(3)=p(3)*(1.0-h2*q(3))
      t(4)=p(4)*(1.0-h2*q(4))
      d(4)=t(4)-t(3)
      sum3=p(3)*p(3)/2.0/eras1
c
c          integrate inward
c
      m2=kkk-idb
      if (m2.lt.40) m2=40
  600 do 620 m=2,m2
      i=i-1
      q(5)=qq(i)
      d(5)=h1*q(4)*p(4)+d(4)
      t(5)=d(5)+t(4)
      p(5)=t(5)/(1.0-h2*q(5))
      if (i.eq.imatch-1) go to 700
      pnlo(i)=p(5)
      do 620 k=1,4
      p(k)=p(k+1)
      t(k)=t(k+1)
      d(k)=d(k+1)
  620 q(k)=q(k+1)
      m2=40
      q(5)=qq(i-2)
      d(5)=h1*q(4)*p(4)+d(4)
      t(5)=d(5)+t(4)
      p(5)=t(5)/(1.0-h2*q(5))
      p(5)=1.09375*p(4)+0.2734375*p(5)-0.546875*p(3)+0.21875*p(2)
     1     -0.0390625*p(1)
      i=i-1
      dx=dx/2.0
      q(5)=qq(i)
      h=-dx
      h1=h*h
      h2=h1/12.0
      t(5)=p(5)*(1.0-h2*q(5))
      t(4)=p(4)*(1.0-h2*q(4))
      d(5)=t(5)-t(4)
      pnlo(i)=p(5)
      do 630 l1=1,4
      p(l1)=p(l1+1)
      t(l1)=t(l1+1)
      d(l1)=d(l1+1)
  630 q(l1)=q(l1+1)
      go to 600
c
c     matching radius has been reached coming in
c            integrate pnlo**2 and calc logarithmic derivative
c
  700 k=kkk
      pnlo(k)=0.0
      call quad5(pnlo,2,imatch,kkk,sum4)
      sum3=sum3+sum4
  760 s3=sum3/p(4)**2
      ppin=(t(5)-t(3)-0.5*(p(5)-p(3)))/h
      s4=ppin/p(4)
c
c     improve trial eigenvalue by perturbation theory if necessary
c
  800 de=(s6-s4)/(s5+s3)
      if (de) 801,900,802
  801 emore=min(emore,e)
      de=max(de,0.8*(eless-e))
      go to 805
  802 eless=max(eless,e)
      de=min(de,0.8*(emore-e))
  805 e=e+de
      if (abs(de/e).gt.thresh) go to 200
      if (kkk.eq.ib) go to 912
      go to 900
c
c     improve trial eigenvalue for finite-boundary case
c
  850 kkk=ib
      ppout=(pnlo(kkk+1)-pnlo(kkk-1))/(r(kkk+1)-r(kkk-1))
      if (mod(kkk-1,40).ne.0) go to 860
      ppout=0.5*(ppout+(pnlo(kkk)-pnlo(kkk-1))/(r(kkk)-r(kkk-1)))
  860 call quad5(pnlo,2,1,kkk,sum1)
      s6=ppout/pnlo(kkk)
      s5=sum1/pnlo(kkk)**2
      s4=1.0/r(kkk)
      s3=0.0
      go to 800
c
c     determine normalization constant for bound wavefunction
c
  900 pop=pmatch/p(4)
      kkk=min0(kkk,ib)
      do 910 j=imatch,kkk
  910 pnlo(j)=pnlo(j)*pop
  912 call quad5(pnlo,2,1,kkk,sum1)
  950 c1=sqrt(sum1)
      if (pnlo(3).gt.0.0) go to 980
      c1=-c1
      go to 980
c
c     determine normalization constant for unbound wavefunction
c         ( pmax=1.0/sqrt(pi*sqrt(e)) )
c
  960 i=i-1
      nmin=i
  962 i=i-1
      inmin=i
      if (pnlo(i+1).eq.0.0) go to 962
      if (pnlo(i)/pnlo(i+1).gt.0.0) go to 962
  964 i=i-1
      if (pnlo(i)/pnlo(i+1).gt.1.0) go to 964
  972 a1=(pnlo(i)-pnlo(i+1))/(r(i)-r(i+1))
      b1=(pnlo(i+1)-pnlo(i+2))/(r(i+1)-r(i+2))
      a1=(a1-b1)/(r(i)-r(i+2))
      b1=b1-a1*(r(i+1)+r(i+2))
      c1=pnlo(i+1)-r(i+1)*(a1*r(i+1)+b1)
      c1=abs(c1-b1*b1/4.0/a1)
      zst=zzz-1.0+1.0e-7
      a1=0.5*zst/(e*r(i+1))
      b0=lam*lamp
      b1=1.0-b0/(2.0*zst*r(i+1))-2.5*a1
      f=1.0-a1*b1
      if ((kut.eq.-1.and.delta.lt.(20.0*tolend)).and.abs(delp(mm)).
     1  lt.4.0e-6) then
         write (9,97) zst,e,r(i+1),a1,b1,f,i,inmin,nmin
   97 format (1h0,9x,23h*****zst,e,r,a1,b1,f,i=,3f13.6,3f13.8,3i7)
      if (nconf.ne.nconf6.and.inmin-i.lt.2) 
     1     write (6,98) nlbcd(ncspvs),i,inmin,
     2     nlbcd(ncspvs),mesh-10,mesh-1,(pnlo(nn),nn=mesh-10,mesh-1)
   98 format (/' *****EMX too small for accurate calc. of continuum ',
     1  'function  P(',a3,') :  i,inmin=',2i6,' must differ by two',
     2   ' or more.' /' P/',a3,'\(',i4,'-',i4,')=   ',10(f9.5)/)
      nconf6=nconf
      end if
      c1=c1/f
      c1=c1*sqrt(3.141592654*sqrt(e))
      kkk=mesh
      imatch=i+1
c
c     calculate the normalized wave function
c
  980 j=1
      if (mod(niter,5).eq.0.and.irel.gt.3) write (9,9872) m,lam,erasx
 9872 format (10x,2i6,6herasx=,f10.6)
      k=2
      c1=1.0/c1
      do 985 i=1,kkk-1
      if (pnlo(j)*pnlo(j+1).ge.0.0) go to 984
      if (k.ge.nodf) go to 985
      k=k+1
  984 j=j+1
  985 pnlo(i)=pnlo(i)*c1
      if (niter.ne.17.or.nn.gt.2) go to 8790
      if (irel.le.3) go to 8790
      qq0(1)=0.0
      qq(1)=0.0
      do 8787 i=2,360
      eras=pnlo(i-1)
      if (mod(i-1,40).eq.0) eras=pnlo(i-2)
      eras=pnlo(i+1)-2.0*pnlo(i)+eras
      qq0(i)=eras/(r(i+1)-r(i))**2
 8787 qq(i)=qq(i)*pnlo(i)
      write (9,8788)
 8788 format (1h1)
      write (9,8786) (i,pnlo(i),qq0(i),qq(i), i=1,51)
 8786 format (3(i4,f11.7,2f14.2))
 8790 continue
      if (nodf.gt.1) go to 990
      j=1
  990 return
      end
      subroutine outpt
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      dimension pqnl(kmsh,ko)
c
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      integer fg
      data hdg1,hdg2,hdg3/6hrc(a0),6hr (a0),6h  pnl(/
c
  800 rewind 4
      read (4) ru,ruee,ruexch,rsatom,recorr,rucorr
c
c     compute first term of series. (pnl(r)/r**(lam+1) at r=0)
c
      do 807 m=1,ncspvs
      lp=l(m)+1
      do 806 i=1,4
      a(i,1)=1.0
      a(i,2)=r(i+1)
      a(i,3)=r(i+1)*r(i+1)
      a(i,4)=r(i+1)*a(i,3)
  806 a(i,5)=pnl(i+1,m)/(r(i+1)**lp)
      call crosym (4)
      az(m)=a(1,5)
  807 continue
      etotrl=0.0
      if (tolend.gt.1.0e-4) go to 987
      if (ihf.eq.0.or.ncspvs.le.1) go to 808
      if (irel.ne.0) go to 990
      go to 982
c
c          calculate one-electron energy integrals
c
  808 twooalf=2.0*137.0360
      froasq=twooalf**2
      do 850 m=1,ncspvs
      read (4) (v(i), i=1,mesh)
      if (m.lt.n1sc) go to 850
  810 kkk=nkkk(m)
      qnl(1,m)=0.0
      xi(1)=0.0
      do 812 i=2,kkk
      xj(i)=pnl(i,m)**2
      eras1=r(i)
      eras=rscore(i)-xj(i)
      if (ee(ncspvs).lt.0.0) eras=eras+rsvale(i)
      call subcor
      recorr(i)=b1
      rucorr(i)=b2
      if (irel.lt.2.or.ee(m).ge.0.0) go to 812
      eras=twooalf/(froasq-v(i)+ee(m))
      qnl(i,m)=eras*((pnl(i+1,m)-pnl(i-1,m))/(r(i+1)-r(i-1))
     1  -pnl(i,m)/r(i))
  812 pnlo(i)=xj(i)/r(i)
      if (irel.ge.2.and.ee(m).lt.0.0) qnl(kkk,m)=0.0
      a1=0.0
  815 do 840 n=1,6
      if (ee(m)) 820,820,840
  820 do 835 i=2,kkk
      go to (821,822,823,824,825,826), n
  821 xi(i)=xj(i)*twoz/r(i)
      go to 835
  822 xi(i)=xj(i)*v(i)
      go to 835
  823 xi(i)=pnlo(i)*ruee(i)
      go to 835
  824 xi(i)=pnlo(i)*ruexch(i)
      go to 835
  825 xi(i)=pnlo(i)*rucorr(i)
      go to 835
  826 xi(i)=pnlo(i)*recorr(i)
  835 continue
      call quad5(xi,1,1,kkk,a1)
  840 b(n)=a1
      ekin(m)=ee(m)-b(2)
      een(m)=-b(1)
      ei(m)=ekin(m)+een(m)
      uee(m)=b(3)
      uex(m)=b(4)*exf10
      ucorr(m)=b(5)
      ecorr(m)=b(6)
      epss(m)=ei(m)+uee(m)+1.5*b(4)+ecorr(m)
  850 continue
      do 870 i=1,mesh
      if (ee(ncspvs).le.0.0) go to 869
      rsatom(i)=rscore(i)
      go to 870
  869 rsatom(i)=rscore(i)+rsvale(i)
  870 continue
      rewind 4
c
c       output
c
  900 eitot=0.0
      ekint=0.0
      eentot=0.0
      eeetot=0.0
      eextot=0.0
      ucorrt=0.0
      ecorrt=0.0
      do 907 m=1,ncspvs
      if (ee(m).ge.0.0) go to 907
      eitot=eitot+ei(m)*wwnl(m)
      ekint=ekint+ekin(m)*wwnl(m)
      eentot=eentot+een(m)*wwnl(m)
      eeetot=eeetot+uee(m)*wwnl(m)
      eextot=eextot+uex(m)*wwnl(m)
      ucorrt=ucorrt+ucorr(m)*wwnl(m)
      ecorrt=ecorrt+ecorr(m)*wwnl(m)
  907 continue
c
      eeetot=0.5*eeetot
      eextot=0.75*eextot/exf10
      etot=eitot+eeetot+eextot+ecorrt
      ecorrk=-ecorrt+3.0*(ucorrt-ecorrt)
      ecorrp=ecorrt-ecorrk
      etotk=ekint+ecorrk
      etotp=etot-etotk
c
      if (ipteb-1) 917,912,908
  908 write (9,84) ncelec,nvelec,noelec,rsint
   84 format (8h1ncelec=,i3,11h    nvelec=,i2,11h    noelec=,i3,
     1  10h    rsint=,f10.6///
     2  '0 m   nl    w  imt kkk    e (ryd)', 6x,5he kin,7x,5he e-n,
     3  6x,5hu e-e,5x,6hu exch,4x,6hu corr,6x,5heps s,
     4  5x,6he corr,18h   no scheq cycles//)
      i1=min0(niter,10)
      do 910 m=1,ncspvs
      if (ifrac.eq.0) then
        write (9,85) m,nlbcd(m),wwnl(m),imat(m),nkkk(m),ee(m),ekin(m),
     1 een(m),uee(m),uex(m),ucorr(m),epss(m),ecorr(m),(nsch(i,m),i=1,i1)
   85 format (i3,1x,a3,f6.0,1x,2i4,3f12.5,f11.5,2f10.5,f12.5,f10.5,2i3,
     1  8i2)
      else
        write (9,86) m,nlbcd(m),wwnl(m),imat(m),nkkk(m),ee(m),ekin(m),
     1 een(m),uee(m),uex(m),ucorr(m),epss(m),ecorr(m),(nsch(i,m),i=1,i1)
   86 format (i3,1x,a3,f7.3,2i4,3f12.5,f11.5,2f10.5,f12.5,f10.5,2i3,8i2)
      end if
  910 continue
      write (9,87) etot,ekint,eentot,eeetot,eextot,ecorrt,etotk,etotp,
     1  ecorrk,ecorrp
   87 format (/1h0,15x,5hetot=,f13.5,1h=,f11.5,f12.5,1h+,2f10.5,f32.5///
     1    22x,10he tot k,p=,f11.5,f12.5,27x,11he corr k,p=,2f10.5)
c
  912 write (9,77)
   77 format (//1h2/'   nl   w',10x,2hee,9x,1hj,8x,4hr(j),9x,2haz/)
      do 915 m=1,ncspvs
      j=jjj(m)
      if (ifrac.eq.0) then
        write (9,78) nlbcd(m),wwnl(m),ee(m),jjj(m),r(j),az(m)
   78 format (1x,a3,f6.0,f15.5,i8,f12.5,f12.3)
      else
        write (9,79) nlbcd(m),wwnl(m),ee(m),jjj(m),r(j),az(m)
   79 format (1x,a3,f7.3,f14.5,i8,f12.5,f12.3)
      end if
  915 continue
  917 read (4) ru,ruee,ruexch,rsatom,recorr,rucorr
      kkk=0
      do 920 m=1,ncspvs
      if (nkkk(m)-kkk) 920,920,919
  919 kkk=nkkk(m)
  920 continue
      do 930 m=1,ncspvs
      k=nkkk(m)+5
      if (k-kkk) 925,925,930
  925 do 926 i=k,kkk,5
  926 pnl(i,m)=0.0
  930 continue
      kkk=min(kkk,mesh-5)
c
      if (iptvu.lt.4) go to 945
      write (9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
   89 format (1h1,3a6,4x,6hnconf=,i3,6x,2hz=,i3,5x,4hion=,i2,5x,4hkut=,
     1  i2,6x,4hexf=,f5.3,5x,6hcorrf=,f5.3,7h   ca1=,f5.3,7h   ca0=,
     2  f5.3///)
  933 write (9,90)
   90 format (25x,12houtput ru(i),10x,9hru e-e(i),11x,9hruexch(i),
     1  11x,9hrucorr(i),11x,9hrecorr(i))
      write (9,91) hdg1
   91 format (4x,a6,4x,2hnc,1x,5(9x,1h1,9x,1h6)//)
      nc=-1
      do 940 i=1,kkk,10
      nc=nc+1
  940 write (9,29) rc(i),nc,ru(i),ru(i+5),ruee(i),ruee(i+5), ruexch(i),
     1  ruexch(i+5),rucorr(i),rucorr(i+5),recorr(i),recorr(i+5)
   29 format (1h ,f9.4,i6,4x,10f10.4)
c
  945 if (norbpt.le.0) go to 980
      do 947 m=1,ncspvs
      do 947 i=6,mesh,10
      pqnl(i,m)=pnl(i,m)
      if (irel.ge.3) pqnl(i,m)=qnl(i-5,m)
  947 continue
  948 write (9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
      k=min0(ncspvs,2)
      write (9,92) (hdg3,nlbcd(m), m=1,k)
   92 format (27x,9hrscore(i),11x,9hrsvale(i),11x,9hrsatom(i),
     1  2(9x,a6,a3,2h) ))
      write (9,91) hdg2
      nc=-1
      do 950 i=1,kkk,10
      nc=nc+1
  950 write (9,93) r(i),nc,rscore(i),rscore(i+5),rsvale(i),rsvale(i+5),
     1  rsatom(i),rsatom(i+5), (pnl(i,m),pqnl(i+5,m), m=1,k)
   93 format (1h ,f9.4,i6,4x,6f10.5,4f10.6)
c
  955 if (k.ge.ncspvs) go to 961
      k=max0(k,ncspvs-norbpt)
      j=k+1
      k=min0(k+5,ncspvs)
      write (9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
      write (9,94) (hdg3,nlbcd(m), m=j,k)
   94 format (16x,5(9x,a6,a3,2h) ))
      write (9,91) hdg2
      nc=-1
      do 960 i=1,kkk,10
      nc=nc+1
  960 write (9,95) r(i),nc, (pnl(i,m),pqnl(i+5,m), m=j,k)
   95 format (1h ,f9.4,i6,4x,10f10.6)
      go to 955
c
  961 if (ee(ncspvs).le.0.0) go to 980
      if (norbpt.ge.6) go to 963
      if (nconft.gt.1.and.mod(nconft,5).ne.0) go to 980
  963 write(9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
      write (9,94) hdg3,nlbcd(ncspvs)
      write (9,28)
   28 format (16h    r (a0)    nc,10x,1h1,9x,1h2,9x,1h3,9x,1h4,9x,1h5,
     1  9x,1h6,9x,1h7,9x,1h8,9x,1h9,8x,2h10//)
      nc=-1
      do 970 imin=1,kkk,10
      nc=nc+1
      imax=imin+9
  970 write (9,95) r(imin),nc, (pnl(i,ncspvs), i=imin,imax)
c
c
c          calc expectation values of r**n, spin-orbit pars,
c            relativistic corrs, and slater integrals
c
  980 lines=54
      if (itpow.eq.1.or.itpow.ge.3) lines=lines-3-niter
      if (ipteb.gt.0) lines=60-13-ncspvs
      if (ipteb.gt.1) lines=30-2-ncspvs
      if (iptvu.ge.4.or.norbpt.gt.0) lines=0
      if (itpow.lt.2) go to 983
      if (lines.gt.ncspvs+9) go to 981
      write (9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
      lines=60
  981 lines=lines-9-ncspvs
  983 call power
      if (lines.gt.ncspvs+6) go to 982
      write (9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
      lines=60
  982 call zeta1
      if (ihf.eq.1.and.ncspvs.gt.1) go to 987
      lines=lines-6-ncspvs
      if (iptvu.gt.0) go to 985
      if (lines.gt.ncspvs+19) go to 985
      write (9,89) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
  985 call sli1
  987 rewind 4
      read (4) ru,ruee,ruexch,rsatom
      rewind 4
c
      do 988 i=1,ncspvs
  988 nnlz(i)=100*nnn(i)+l(i)
      norb=1+ncspvs-ns
      if (nnn(ncspvs).gt.98) ion=ion-1
      if (iz.eq.izo.and.ion.eq.iono) go to 989
      izo=iz
      iono=ion
      ieras=iabs(norbpt)
      norb=min0(ncspvs,max0(2,ieras))
      if (ieras.eq.9) norb=ncspvs
      ns=1+ncspvs-norb
  989 if (ihf.ne.0.and.ncspvs.gt.1.and.(ncspvs.gt.2.or.nnn(2).le.90))
     1  go to 990
      if (nconf.eq.0.or.ifrac.eq.1) go to 995
      write (2) nnn,l,wwnl,ncspvs,vprm,conf,nconf,iz,z,ion,irel,hxid,
     1  mesh,c,idb,exf,corrf,kut,npar,nlbcd,nnlz,nkkk,ee,ee8,norb,
     2  r,ru,((pnl(i,m),i=1,kmsh),m=ns,ncspvs),iw6
      go to 995
  990 call hfwrtp
c
  995 return
      end
      subroutine vinti(ivinti)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      common/lc2/ra(kmsh),rb(kmsh)
c
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
      integer fg
c
      if (ncspvs.eq.1) return
      write (9,10) conf
   10 format ('1   vinti integrals for  ',3a6////
     1  10x,'vinti integral',23x,'j**2',17x,'wt*(j**2)'/)
      sum=0.0
c
      do 390 n=1,ncspvs
      fl=l(n)
      lm1=l(n)-1
      flm1=lm1
      do 380 np=1,ncspvs
      if (l(np).ne.lm1) go to 380
      kmax=min0(nkkk(n),nkkk(np))
      kmxm2=kmax-2
      h12=1.0/(12.0*r(2))
c
      do 290 i=1,kmxm2
      if (i.ge.kmax-1) go to 130
      if (i.ge.idb+2) go to 150
      i1=i-40*((i-2)/40)
      if (i1.gt.2) go to 120
      if (i1.eq.2) go to 110
      der=-25.0*pnl(i,np)+48.0*pnl(i+1,np)-36.0*pnl(i+2,np)
     1  +16.0*pnl(i+3,np)-3.0*pnl(i+4,np)
      go to 200
  110 der=-3.0*pnl(i-1,np)-10.0*pnl(i,np)+18.0*pnl(i+1,np)
     1  -6.0*pnl(i+2,np)+pnl(i+3,np)
      go to 200
  120 if (i1.lt.40) go to 150
      if (i1.eq.41) go to 130
      der=-pnl(i-3,np)+6.0*pnl(i-2,np)-18.0*pnl(i-1,np)
     1  +10.0*pnl(i,np)+3.0*pnl(i+1,np)
      go to 200
  130 der=3.0*pnl(i-4,np)-16.0*pnl(i-3,np)+36.0*pnl(i-2,np)
     1  -48.0*pnl(i-1,np)+25.0*pnl(i,np)
      go to 200
  150 der=pnl(i-2,np)-8.0*pnl(i-1,np)+8.0*pnl(i+1,np)-pnl(i+2,np)
  200 der=h12*der
      xi(i)=pnl(i,n)*(der-fl*pnl(i,np)/r(i))
      if (ivinti.lt.2) go to 280
      if (n.eq.3.and.np.eq.2) write (9,9920) i,r(i),pnl(i,np),
     1  der,i1,h12,pnl(i,n),xi(i)
 9920 format (i5,3f14.8,i5,3f14.8)
  280 if (i1.eq.41) h12=0.5*h12
  290 continue
      xi(kmax-1)=xi(kmxm2)
      xi(kmax)=xi(kmxm2)
c
      call quad5(xi,1,1,kmax,fj)
      fjsq=fj**2
      aa=(4.0*fl+2.0)*(4.0*flm1+2.0)
      wt=2.0*fl*wwnl(n)*wwnl(np)/aa
      bb=wt*fjsq
      sum=sum+bb
      write (9,30) nlbcd(n),nlbcd(np),fj,fjsq,bb
   30 format (5x,'j(',a3,',',a3,')=',f14.8,8x,f15.8,8x,f15.8)
  380 continue
  390 continue
c
      write (9,40) sum
   40 format (/10x,'total k-factor=',f16.8)
      return
      end
      subroutine phshift
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      dimension fmz(5,2),rz(5,2),phase(5,2),shift(5)
c
      pi=3.141592653590
      twopi=2.0*pi
      twozz=twozzz-2.0
      mm=ncspvs
      eras=l(mm)
      aa=0.5*pi*eras
      fk=sqrt(ee(mm))
      n1=5
      i1=1
      do 120 i=1000,mesh
  120 pnlo(i)=pnl(i,mm)
      go to 300
c
c          calculate (non-relativistic) coulomb function
c
  200 do 210 i=2,mesh
      v(i)=-twozz/r(i)
  210 x2(i)=v(i)
      nn=nnn(mm)
      lam=l(mm)
      e=ee(mm)
      irel1=irel
      irel=0
      call scheq
      irel=irel1
c
c          calc distorted-wave or coulomb phase shift
c
  300 m=mesh-5
      do 330 n=1,n1
  310 m=m-1
      if (pnlo(m).lt.0.0) go to 310
      if (pnlo(m-1).gt.0.0) go to 310
      m=m-1
      ll=m
      if (abs(pnlo(m)).gt.pnlo(m+1)) ll=m+1
c          find node via r=c1+c2*p+c3*(p**3)
      a1=pnlo(ll-1)
      a2=pnlo(ll)
      a3=pnlo(ll+1)
      a13=a1**3
      a23=a2**3
      a33=a3**3
      b1=r(ll-1)
      b2=r(ll)
      b3=r(ll+1)
      b0=b1*(a2*a33-a3*a23)+b2*(a3*a13-a1*a33)+b3*(a1*a23-a2*a13)
      a0=a1*(a23-a33)+a2*(a33-a13)+a3*(a13-a23)
      rz(n,i1)=b0/a0
      f1=m
      fmz(n,i1)=f1+(rz(n,i1)-r(m))/(r(m+1)-r(m))
      eras=fk*rz(n,i1)+0.5*twozz*log(2.0*fk*rz(n,i1))/fk-aa
      eras=mod(eras,twopi)
      if (eras.le.0.0) eras=eras+twopi
      phase(n,i1)=eras
      if (i1.ne.2) go to 330
      shift(n)=phase(n,1)-phase(n,2)
      if (shift(n).lt.0.0) shift(n)=shift(n)+twopi
  330 continue
      if (i1.eq.2) go to 400
      i1=2
      go to 200
c
  400 write (9,40) (fmz(n,1),rz(n,1),phase(n,1),fmz(n,2),rz(n,2),
     1  phase(n,2),shift(n), n=1,n1)
   40 format (//2x,7hmesh pt,5x,1hr,5x,8hdw phase,7x,7hmesh pt,5x,
     1  1hr,6x,7hcoul ph,4x,8hph shift/(f9.2,f9.3,f10.5,f14.2,f9.3,
     2  f10.5,f12.5))
      return
      end
      subroutine subcor
c
c          calculate correlation energy density
c
      implicit real*8 (a-h,o-z)
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
c
  100 if (eras.gt.1.0e-20) go to 120
      b0=9999.9999
      b1=0.0
      b2=0.0
      b3=0.0
      go to 200
  120 b0=(3.0*eras1*eras1/eras)**.333333333e+00
      a1=sqrt(b0+9.0)
      a2=4.0*a1+1.142*b0
      ec=-1.0/a2
      b1=eras1*ec
      a2=ec*ec*(2.0/a1+1.142)
      b2=eras1*(ec-b0*a2/3.0)
      b3=-1.22177412*eras1/b0
      ec=b2/eras1
      b1=b2
      ec1=-1.0/(4.0*a1+0.75*1.142*b0)
      b2=eras1*ec1
      b1=b2
  200 return
      end
      subroutine crosym(m)
c
c     simultaneous equation solver
c     written by i.c. hanson, scientific computation department,
c     lockheed missles and space company, sunnyvale, california
c     solve m simultaneous equations by the method of crout
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
c
      n=m+1
      i1=1
    1 i3=i1
      sum=abs(a(i1,i1))
      do 3 i=i1,m
      if (sum-abs(a(i,i1))) 2,3,3
    2 i3=i
      sum=abs(a(i,i1))
    3 continue
      if (i3-i1) 4,6,4
    4 do 5 j=1,n
      sum=-a(i1,j)
      a(i1,j)=a(i3,j)
    5 a(i3,j)=sum
    6 i3=i1+1
      do 7 i=i3,m
    7 a(i,i1)=a(i,i1)/a(i1,i1)
   14 j2=i1-1
      i3=i1+1
      if (j2) 8,11,8
    8 do 9 j=i3,n
      do 9 i=1,j2
    9 a(i1,j)=a(i1,j)-a(i1,i)*a(i,j)
      if (i1-m) 11,15,11
   11 j2=i1
      i1=i1+1
      do 12 i=i1,m
      do 12 j=1,j2
   12 a(i,i1)=a(i,i1)-a(i,j)*a(j,i1)
      if (i1-m) 1,14,1
   15 do 17 i=1,m
      j2=m-i
      i3=j2+1
      a(i3,n)=a(i3,n)/a(i3,i3)
      if (j2) 16,18,16
   16 do 17 j=1,j2
   17 a(j,n)=a(j,n)-a(i3,n)*a(j,i3)
   18 return
      end
      subroutine quad5(fi,ipow,i1,i2,a1)
c
c     integrate fi (ipow=1) or fi**2 (ipow=2) from mesh point i1 to i2.
c          integration is by 8 applications of newton-cotes closed
c            quadrature for five intervals on each block.
c if i2-1 not a multiple of 5, fill out remainder
c   of mesh by simpsons and/or trapezoidal rule.
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      dimension fi(kmsh)
c
      mmsave=mm
      xif=(r(i1+1)-r(i1))*2.5/288.0
      a1=0.0
      j=i1
      k=j-40*(j/40)
      mm=8-k/5
      i3=5*((i2-1)/5)+1
      if (ipow.eq.2) go to 200
      go to 106
c
  105 mm=8
  106 mm1=(i2-j)/5
      mm=min0(mm,mm1)
      if (j.eq.idb) mm=mm1
      xif=xif+xif
      a2=0.0
  110 a2=a2+19.0*(fi(j)+fi(j+5))+75.0*(fi(j+1)+fi(j+4))
     1   +50.0*(fi(j+2)+fi(j+3))
      j=j+5
      mm=mm-1
      if (ib.ge.0) go to 115
      a3=a1+xif*a2
      write (9,10) ib,j,r(j),fi(j),xif,a2,a1,a3
   10 format (2i6,f12.5,1p,5e16.6)
  115 if (mm.gt.0) go to 110
      a1=a1+xif*a2
      if (j.lt.i3) go to 105
      if (j.eq.i2) go to 150
      xif=(r(i2)-r(i2-1))/3.0
  120 if (j.eq.i2-1) go to 130
      a1=a1+xif*(fi(j)+4.0*fi(j+1)+fi(j+2))
      j=j+2
      if (j.eq.i2) go to 150
      go to 120
  130 a1=a1+1.5*xif*(fi(j)+fi(j+1))
  150 mm=mmsave
      return
c
  200 fj5=fi(j)**2
      go to 206
  205 mm=8
  206 mm1=(i2-j)/5
      mm=min0(mm,mm1)
      if (j.eq.idb) mm=mm1
      xif=xif+xif
      a2=0.0
  210 fj=fj5
      fj5=fi(j+5)**2
      a2=a2+19.0*(fj+fj5)+75.0*(fi(j+1)**2+fi(j+4)**2)
     1  +50.0*(fi(j+2)**2+fi(j+3)**2)
      j=j+5
      mm=mm-1
      if (mm.gt.0) go to 210
      a1=a1+xif*a2
      if (j.lt.i3) go to 205
      if (j.eq.i2) go to 250
      xif=(r(i2)-r(i2-1))/3.0
  220 if (j.eq.i2-1) go to 230
      fj=fj5
      fj5=fi(j+2)**2
      a1=a1+xif*(fj+4.0*fi(j+1)**2+fj5)
      j=j+2
      if (j.eq.i2) go to 250
      go to 220
  230 a1=a1+1.5*xif*(fj5+fi(j+1)**2)
  250 mm=mmsave
      return
      end
      subroutine quad2(m)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
c
c     calculate atomic coulomb potential  (integrate by simpsons rule)
c           xi=integral(rsatom),  xj=integral(rsatom/r)
c
  510 xi(1)=0.0
      xj(1)=0.0
      b2=0.0
      im=min0(mesh,ib+6)
c
      ho3=r(2)/3.0
      i=1
      do 515 j=1,nblock
      do 513 k=1,20
      i=i+2
      if (i.gt.ib) go to 520
      xi(i)=xi(i-2)+ho3*(rsatom(i-2)+4.0*rsatom(i-1)+rsatom(i))
      b0=b2
      b2=rsatom(i)/r(i)
  513 xj(i)=xj(i-2)+ho3*(b0+4.0*rsatom(i-1)/r(i-1)+b2)
      if (i.gt.idb) go to 515
      ho3=ho3+ho3
  515 continue
  520 ic=1
      if (i.eq.ib) go to 525
      i=i-1
      if (i.ne.ib) go to 525
      ic=2
      xi(i)=xi(i-1)+1.5*ho3*(rsatom(i-1)+rsatom(i))
      b0=b2
      b2=rsatom(i)/r(i)
      xj(i)=xj(i-1)+1.5*ho3*(b0+b2)
  525 if (ib.eq.mesh) go to 530
      do 527 i=ib,im
      xi(i)=xi(ib)
  527 xj(i)=xj(ib)
  530 rsint1=xi(ib)
      if (m.eq.2) rsint=rsint1
c
c           ruee=r*(electron-electron energy), xj=r*(total energy)
c
      do 535 i=1,im,2
  535 xi(i)=2.0*(xi(i)+r(i)*(xj(ib)-xj(i)))
      if (ic.eq.2) xi(ib)=2.0*xi(ib)
      do 540 i=3,ib,4
      b0=(xi(i)-xi(i-2)+xi(i)-xi(i+2))/8.0
      xi(i-1)=0.5*(xi(i)+xi(i-2))+b0
  540 xi(i+1)=0.5*(xi(i)+xi(i+2))+b0
      i1=mod(ib-1,4)
      if (i1.ne.1) go to 545
      i=ib-i1
      b0=(xi(i)-xi(i-2)+xi(i)-xi(i+2))/8.0
      xi(i+1)=0.5*(xi(i)+xi(i+2))+b0
  545 if (m.eq.1) return
      do 550 i=1,ib
      ruee(i)=xi(i)
  550 xj(i)=-twoz+ruee(i)+exf10*ruexch(i)
      do 560 i=ib,mesh
      ruee(i)=ruee(ib)
  560 xj(i)=xj(ib)
      return
      end
      subroutine clock
c
      implicit real*8 (a-h,o-z)
      common/c7/t0,t1,t2
c
      call seconds(t2)
      tr=t2-t1
      tj=t2-t0
      write (9,99) tr,tj
   99 format (//12h0    t(run)=,f6.2,1h,,5x,7ht(job)=,f7.2,8h seconds)
      return
      end
      subroutine seconds(t)
c
c          Timing routine (results for information only)
c               For other machines, set t=0.0 or add
c               appropriate code.
c
      real*8 t
c
c------------------------------------
c          Use this code for a CRAY
c     t=second(t)
c------------------------------------
c          Use this code for a VAX or Macintosh
      real t1
      t1=secnds(0.0)
      t=dble(t1)
c------------------------------------
c          Use this code for a PC
c     real*4 t1,second
c     t1=second()
c     t=dble(t1)
c------------------------------------
c          Use this code for a SUN
c     real etime,t2
c     dimension t2(2)
c     t=etime(t2)
c-------------------------------------
c          Use this code for IBM RISC
c     it=mclock()
c     t=it
c     t=t/100.0
c-------------------------------------
      return
      end
      subroutine power
c
c          evaluate expectation values of r to n
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
      dimension s1(10,ko),s2(10,ko)
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
      dimension psq(kmsh),p(kmsh,10)
c     equivalence (xj,psq,p)
c
      integer fg
c
      mm=ncspvs
      if (ee(mm).gt.0.0) mm=mm-1
      if (mm.le.0) return
  350 if (itpow.ge.2. and.end.eq.1.0) write (9,35)
   35 format (///2x,'nl   wnl  ',7x,2hee,10x,2haz,
     1    6x,5h(r-3),6x,5h(r-2),6x,5h(r-1),6x,5h(r+1),6x,5h(r+2),
     2    6x,5h(r+3),6x,5h(r+4),6x,5h(r+6)/)
c
  400 kkk=0
      do 410 m=1,mm
      if (nkkk(m).gt.kkk) kkk=nkkk(m)
  410 continue
c
c          integrate
c
      do 420 k=1,kkk
      r1=r(k)
      r2=r1*r1
      r3=r1*r2
      p(k,2)=1.0/r3
      p(k,3)=1.0/r2
      p(k,4)=1.0/r1
      p(k,5)=r1
      p(k,6)=r2
      p(k,7)=r3
      p(k,8)=r2*r2
      p(k,9)=r2*r3
  420 p(k,10)=r3*r3
      do 599 m=1,mm
      kkk=nkkk(m)
      do 510 k=1,kkk
  510 psq(k)=pnl(k,m)**2
      do 530 n=2,10
      s2(n,m)=0.0
      do 520 k=1,kkk
  520 xi(k)=psq(k)*p(k,n)
      if (n.ne.3.or.l(m).ne.0) go to 530
      xi(1)=4.0*(xi(2)+xi(4))-6.0*xi(3)-xi(5)
  530 call quad5(xi,1,1,kkk,s1(n,m))
      rm3(m)=s1(2,m)
  599 continue
c               print <r> as a diagnostic (npr.ne.0.and.niter.ge.np)
      if (end.eq.1.0) go to 700
      n1=max(1,ncspvs-6)
      write (6,70) (nlbcd(n),s1(5,n), n=n1, ncspvs)
   70 format (' <r>=',7(4x,a3,f7.3))
      go to 900
c
  700 if (itpow.lt.2) go to 900
c
      fn=1.0e-20
      do 721 m=1,mm
  710 do 711 n=2,10
  711 s2(n,1)=s2(n,1)+wwnl(m)*s1(n,m)
      fn=fn+wwnl(m)
      s1(9,m)=s1(10,m)
      if (l(m).le.0) s1(2,m)=0.0
      if (ifrac.eq.0) then
        write (9,71) nlbcd(m),wwnl(m),ee(m),az(m),(s1(n,m),n=2,9)
   71 format (1x,a3,f6.0,f14.5,f11.3,1p,8e11.3)
      else
        write (9,72) nlbcd(m),wwnl(m),ee(m),az(m),(s1(n,m),n=2,9)
   72 format (1x,a3,f7.3,f13.5,f11.3,1p,8e11.3)
      end if
  721 continue
      if (mm.eq.1) return
      do 730 n=3,10
  730 s2(n,2)=s2(n,1)/fn
      r1=0.1
      s2(9,1)=s2(10,1)
      s2(9,2)=s2(10,2)
      if (ifrac.eq.0) then
        write (9,73)fn,(s2(n,1), n=3,9), r1,(s2(n,2),n=3,9)
   73 format (/4x,f6.0,36x,1p,7e11.3)
      else
        write (9,74)fn,(s2(n,1), n=3,9), r1,(s2(n,2),n=3,9)
   74 format (/4x,f7.3,35x,1p,7e11.3)
      end if
c
  900 return
      end
      subroutine zeta1
c
c          calculate spin-orbit parameters (via both central-field
c               formula and blume-watson method)
c                     and relativistic energy corrections
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
      dimension rsq(kmsh),zeta(ko,3)
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      integer fg
      dimension fi(kmsh,5),d(kmsh,3),drudr(kmsh),p(kmsh)
c     equivalence (rscore,fi),(pnlo,d),(ruee,p)
c
  100 npar=1
      nzi=0
      nfij=0
      ngij=0
      do 101 i=1,45
      vprm(i)=0.0
  101 kpar(i)=0
c
      if (ihf.eq.1) go to 360
  350 write (9,35)
   35 format(///1h ,11x,'--------------------zeta---------------------'
     1/'  nl   wnl ',2x,21h ----blume-watson----,5x,16h-----(r*vi)-----,
     2  5x,7hi (ryd),6x,5heps s,7x,4hevel,7x,24hedar   evel+edar    erel
     3  /16x,6h (ryd),6x,6h(cm-1),5x,5h(ryd) ,6x,6h(cm-1))
  360 do 380 i=1,mesh
      rsq(i)=r(i)*r(i)
      d(i,1)=-z-z+ruee(i)
  380 d(i,2)=ru(i)
c
c           calculate evel, edar, and zeta for each orbital
c
  400 etotrl=0.0
      do 899 m=1,ncspvs
      if (ee(m)) 405,405,401
  401 evel(m)=0.0
      edar(m)=0.0
      zeta(m,1)=0.0
      zeta(m,2)=0.0
      zeta(m,3)=0.0
      go to 700
c
c          set up potentials for individual orbitals and total atom
c
  405 read (4)(v(i), i=1,mesh)
      do 406 i=1,mesh
  406 d(i,3)=v(i)*r(i)
c
c          calculate derivatives of potentials
c
      n1=3
      nblock=mesh/40
      do 460 n=n1,3
      dr=3.0*r(2)
      i=1
      do 420 m123=1,nblock
      do 410 j=1,10
      i=i+4
      b0=d(i-3,n)-d(i-2,n)-d(i-2,n)+d(i-1,n)
      b1=d(i-2,n)-d(i-1,n)-d(i-1,n)+d(i,n)
      b2=d(i-3,n)-d(i,n)
      drudr(i-3)=(-4.5*b0-b2)/dr
      drudr(i-2)=(-1.5*b1-b2)/dr
      drudr(i-1)=(1.5*b0-b2)/dr
  410 drudr(i)=(4.5*b1-b2)/dr
      if (i.gt.idb) go to 420
      dr=dr+dr
  420 continue
      do 450 i=2,mesh
      v(i)=d(i,n)/r(i)
  450 d(i,n)=(drudr(i)-v(i))/rsq(i)
  460 continue
c
c          calc p=dpnl/dr and all integrands
c
      dr=3.0*r(2)
      kkk=nkkk(m)
      nblokk=kkk/40
      i=1
      p(1)=0.0
      do 515 nb=1,nblokk
      do 514 j=1,10
      i=i+4
      b0=pnl(i-3,m)-2.0*pnl(i-2,m)+pnl(i-1,m)
      b1=pnl(i-2,m)-2.0*pnl(i-1,m)+pnl(i,m)
      b2=pnl(i-3,m)-pnl(i,m)
      p(i-3)=(-4.5*b0-b2)/dr
      p(i-2)=(-1.5*b1-b2)/dr
      p(i-1)=(1.5*b0-b2)/dr
  514 p(i)=(4.5*b1-b2)/dr
      if (i.gt.idb) go to 515
      dr=dr+dr
  515 continue
c
      eras=(0.5/137.036)**2
      do 520 i=1,kkk
      psq=pnl(i,m)*pnl(i,m)
      fi(i,1)=psq*d(i,1)
      fi(i,2)=psq*d(i,2)
      fi(i,3)=psq*d(i,3)
      fi(i,4)=psq*((ee(m)-v(i))**2)
      fi(i,5)=pnl(i,m)*(r(i)*p(i)-pnl(i,m))*d(i,3)
      eras1=1.0+eras*(ee(m)-v(i))
      fi(i,2)=fi(i,3)/eras1
      if (irel.ne.0) fi(i,5)=fi(i,5)/eras1
  520 continue
      fi(1,4)=4.0*(fi(2,4)-1.5*fi(3,4)+fi(4,4))-fi(5,4)
      fi(1,5)=4.0*(fi(2,5)-1.5*fi(3,5)+fi(4,5))-fi(5,5)
c
c          integrate
c
      p(1)=0.0
      p(2)=0.0
      p(3)=0.0
      ibo=ib
      n1=3
      if (l(m).le.0) n1=4
      do 530 n=n1,5
      if (itpow.gt.5.and.n.ge.4) ib=-m
  530 call quad5(fi(1,n),1,1,kkk,p(n))
      ib=ibo
  550 zeta(m,1)=2.66246e-5*p(1)
      zeta(m,2)=2.66246e-5*p(2)
      zeta(m,3)=2.66246e-5*p(3)
      if (m.lt.n1sc) go to 700
      evel(m)=-1.33123e-5*p(4)
      edar(m)=-1.33123e-5*p(5)
c
  700 zetak1=zeta(m,1)*109737.31
      zetak1=zeta(m,2)
      zetak2=zeta(m,2)*109737.31
      zetak3=zeta(m,3)*109737.31
      ec=evel(m)
      if (l(m).eq.0) ec=ec+edar(m)
      etotrl=etotrl+wwnl(m)*ec
      if (ihf.eq.1) go to 705
      erll=epss(m)+ec
c          call zetabw to calc spin-orbit integral by the blume-watson m
  701 if (ee(m).lt.0.0) zetak2=zetabw(m,rm3(m))
  702 zetak1=zetak2/109737.3
      if (ifrac.eq.0) then
        write (9,70) nlbcd(m),wwnl(m),zetak1,zetak2,zeta(m,3),zetak3,
     1  ei(m),epss(m),evel(m),edar(m),ec,erll
   70 format (1x,a3,f6.0,f13.5,f12.3,f10.5,f12.3,f11.5,f12.5,2f11.5,
     1  f9.3,f12.5)
      else
        write (9,71) nlbcd(m),wwnl(m),zetak1,zetak2,zeta(m,3),zetak3,
     1  ei(m),epss(m),evel(m),edar(m),ec,erll
   71 format (1x,a3,f7.3,f12.5,f12.3,f10.5,f12.3,f11.5,f12.5,2f11.5,
     1  f9.3,f12.5)
      end if
c
c          store zeta for rcg input
c
  705 if (l(m).le.0) go to 899
      if (wwnl(m).eq.0.0) go to 899
      eras=4*l(m)+2
      if (wwnl(m)-eras) 720,899,899
  720 nzi=nzi+1
      zi(nzi)=zetak3*0.001
  721 zi(nzi)=zetak2*0.001
c
  899 continue
      if (irel.ne.0) etotrl=0.0
      return
      end
      subroutine sli1
c
c          evaluate slater radial coulomb integrals fk and gk
c               and radial overlap integrals
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
      common/sli/ xa(kmsh),xb(kmsh),f1(kmsh),f2(kmsh),rkp1(kmsh),
     1  vpar(7,ko,ko),emn(ko,ko),rryd(10,2),rkay(10,2),k1(10,2),
     2  frac(10,2)
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      integer fg
c     equivalence (v,rkp1),(rscore,vpar)
      data jhlf,jhlg/102,103/
      data jhuF,jhuG/70,71/
c
      npr1=0
      kdmin=0
      eterm=0.0
  120 do 969 ma=1,ncspvs
      la=l(ma)
      i1max=nkkk(ma)
      do 959 mb=ma,ncspvs
      lb=l(mb)
      i2max=nkkk(mb)
      imax=max0(i1max,i2max)
      kemin=iabs(la-lb)
      nokmax=min0(la,lb)+1
      ne=nokmax
      nd=nokmax
      if (ma.ne.mb) go to 400
      ne=0
      if (wwnl(ma).le.1.0) nd=0
c
c          calc pnl-products for direct or exchange case
c
  400 do 700 n1=1,2
      if (n1.le.1) go to 420
      if (ne.le.0) go to 500
      do 415 i=1,imax
      f1(i)=pnl(i,ma)*pnl(i,mb)
  415 f2(i)=f1(i)
      k=kemin-2
      go to 500
  420 do 425 i=1,imax
      f1(i)=pnl(i,ma)**2
  425 f2(i)=pnl(i,mb)**2
      k=kdmin-2
      if (nd.le.0) go to 800
c
c          calc for each value of k
c
  500 do 700 n=1,nokmax
      a1=0.0
      b1=1.0
      k=iabs(k+2)
      k1(n,n1)=k
      if (ee(ma).ge.0.0.or.ee(mb).ge.0.0) go to 650
      if (n1.lt.ne+2) go to 520
      k1(n,n1)=0
      go to 650
  520 continue
c
      rkp1(1)=0.001
      do 530 i=2,imax
      if (n.gt.1) go to 529
      rkp1(i)=r(i)
      if (k.le.0) go to 530
      do 527 m=1,k
  527 rkp1(i)=rkp1(i)*r(i)
      go to 530
  529 rkp1(i)=rkp1(i)*r(i)*r(i)
  530 continue
c
c          calc r1 integral
c
      xi(1)=0.0
      xj(1)=0.0
      xa(1)=0.0
      xb(1)=0.0
      a2=0.0
      b2=0.0
      ho12=r(2)/12.0
      ni=40
      do 550 i=3,i1max,2
      a0=a2
      b0=b2
      eras=8.0*f1(i-1)
      a1=eras*rkp1(i-1)/r(i-1)
      b1=eras/rkp1(i-1)
      a2=f1(i)*rkp1(i)/r(i)
      b2=f1(i)/rkp1(i)
      eras=ho12*(5.0*a0+a1-a2)
      xi(i-1)=xi(i-2)+eras
      xa(i-1)=xa(i-2)+abs(eras)
      eras=ho12*(5.0*b0+b1-b2)
      xj(i-1)=xj(i-2)+eras
      xb(i-1)=xb(i-2)+abs(eras)
      eras=ho12*(-a0+a1+5.0*a2)
      xi(i)=xi(i-1)+eras
      xa(i)=xa(i-1)+abs(eras)
      eras=ho12*(-b0+b1+5.0*b2)
      xj(i)=xj(i-1)+eras
      xb(i)=xb(i-1)+abs(eras)
      ni=ni-2
      if (ni) 540,540,550
  540 ho12=ho12+ho12
      ni=40
      if (i.eq.idb) ni=10000
  550 continue
c
      do 560 i=2,i1max
      eras=rkp1(i)/r(i)
      xi(i)=xi(i)/rkp1(i)+(xj(i1max)-xj(i))*eras
  560 xa(i)=xa(i)/rkp1(i)+(xb(i1max)-xb(i))*eras
      if (i2max-i1max) 600,600,569
  569 a0=xi(i1max)*rkp1(i1max)
      b0=xa(i1max)*rkp1(i1max)
      do 570 i=i1max,i2max
      xi(i)=a0/rkp1(i)
  570 xa(i)=b0/rkp1(i)
c
c          calc r2 integral
c
  600 a0=0.0
      b0=0.0
      a1=0.0
      b1=0.0
      ho3=r(2)/1.5
      ni=40
      do 620 i=3,i2max,2
      a2=f2(i)*xi(i)
      b2=abs(f2(i))*xa(i)
      eras=4.0*f2(i-1)
      a1=a1+ho3*(a0+eras*xi(i-1)+a2)
      b1=b1+ho3*(b0+abs(eras)*xa(i-1)+b2)
      a0=a2
      b0=b2
      ni=ni-2
      if (ni.gt.0) go to 620
      ho3=ho3+ho3
      ni=40
      if (i.eq.idb) ni=10000
  620 continue
c
  650 frac(n,n1)=a1/b1
      rryd(n,n1)=a1
      rkay(n,n1)=a1*109737.31
      if (nfg.eq.0) go to 700
c          calc ls-term contribution to e(av)
      nla=nlbcd(ma)
      nlb=nlbcd(mb)
      do 670 j=1,nfg
      if ((n1.eq.1.and.fg(j,ii).ne.jhlf.and.fg(j,ii).ne.jhuF).or.
     1    (n1.eq.2.and.fg(j,ii).ne.jhlg.and.fg(j,ii).ne.jhuG))
     1  go to 670
      if (k1(n,n1).ne.kfg(j,ii)) go to 670
      if ((ifg(j,ii).ne.nla.or.jfg(j,ii).ne.nlb).and.
     1    (jfg(j,ii).ne.nla.or.ifg(j,ii).ne.nlb)) go to 670
      eterm=eterm+cfg(j,ii)*a1
  670 continue
  700 continue
      if (iptvu.lt.3) go to 840
c
c          calc overlap integrals
c
  800 i=1
      xif=(r(imax)-r(imax-1))*5.0/288.0
      eras=0.0
      if (ee(ma).ge.0.0.or.ee(mb).ge.0.0) go to 830
  810 if (i.le.idb) eras=0.5*eras
      do 820 n=1,8
      i=i+5
  820 eras=eras+19.0*(f1(i-5)+f1(i))+75.0*(f1(i-4)+f1(i-1))
     1         +50.0*(f1(i-3)+f1(i-2))
      if (i-imax) 810,830,830
  830 vpar(7,mb,ma)=eras*xif
      if (nd.le.0) go to 959
  840 if (iptvu.eq.0) go to 915
      if (ma.lt.ncspvs-iabs(norbpt)+1) go to 915
c
c          output
c
  899 if (npr1) 900,900,905
  900 write (9,90) conf,nconf,iz,ion,kut,exf,corrf,ca1,ca0
   90 format (17h1slater integrals///1x,3a6,6x,6hnconf=,i3,6x,2hz=,i3,
     1  5x,4hion=,i2,5x,4hkut=,i2,6x,4hexf=,f5.3,5x,6hcorrf=,f5.3,3x,
     2  4hca1=,f5.3,7h   ca0=,f5.3//)
  901 write (9,91)
   91 format (1h0,16x,1hk,19x,2hfk,21x,4hfrac,9x,1hk,19x,2hgk,21x,
     1  4hfrac)
      write (9,92)
   92 format (1h ,16x,39h-   -----------------------------------,
     1  8h   -----,9x,39h-   -----------------------------------,
     2  8h   -----//)
  905 npr1=npr1+nokmax+1
      if (npr1-50) 910,910,906
  906 npr1=0
      go to 900
  910 do 911 n=1,nokmax
  911 write (9,93) nlbcd(ma),nlbcd(mb), (k1(n,n1),rryd(n,n1),rkay(n,n1),
     1  frac(n,n1), n1=1,2)
   93 format (2h (,a3,1h,,a3,1h),i8,f13.7,7h ryd  =,f13.3,5h cm-1,f8.3,
     1  i10,f13.7,7h ryd  =,f13.3,5h cm-1,f8.3)
      write (9,94)
   94 format (1h )
c
c          store slater integrals for rcg input
c
  915 a1=4*l(ma)+1
      a2=4*l(mb)+1
      if (wwnl(ma).gt.a1) go to 950
      if (ma.ne.mb) go to 930
      if (wwnl(ma).le.1.0) go to 950
      if (wwnl(ma).ge.a1) go to 950
      do 924 i=2,nokmax
      npar=npar+1
      kpar(npar)=1
  924 vprm(npar)=rkay(i,1)*0.001
      go to 950
  930 if (wwnl(mb).gt.a2) go to 950
      if (wwnl(ma).eq.0.0.or.wwnl(mb).eq.0.0) go to 950
      do 934 i=1,nokmax
      if (k1(i,1).le.0) go to 933
      nfij=nfij+1
      fij(nfij)=rkay(i,1)*0.001
  933 ngij=ngij+1
  934 gij(ngij)=rkay(i,2)*0.001
c
  950 if (ma.ne.mb) go to 955
      do 953 i=1,nokmax
  953 vpar(i,ma,mb)=rryd(i,1)
      go to 959
  955 vpar(1,ma,mb)=rryd(1,1)
      do 956 i=1,nokmax
  956 vpar(i+1,ma,mb)=rryd(i,2)
  959 continue
c
  969 continue
c
      call rcn3s
  980 return
      end
      subroutine rcn3s
c
c          calculate ionization potentials and total energies
c           from values of slater coulomb integrals
c
c
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/brt/emag,eret,fmi(10,3),fni(10,4),fnorm(ko)
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
      common/sli/ xa(kmsh),xb(kmsh),f1(kmsh),f2(kmsh),rkp1(kmsh),
     1  vpar(7,20,ko),emn(ko,ko),rryd(10,2),rkay(10,2),k1(10,2),
     2  frac(10,2)
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      integer fg
c     equivalence (v,rkp1),(rscore,vpar)
c
      rdb=r(idb)
      if (iptvu.lt.3) go to 200
c
c          print overlap integrals
c
  100 write (9,30) conf,nconf,iz,ncores,nvales,ion,c,ca1
      write (9,10)
   10 format (18h0overlap integrals)
      n2=min0(ncspvs,12)
      write (9,33) (nlbcd(n), n=1,n2)
      write (9,35)
      do 120 m=1,ncspvs
      n2=min0(m,12)
  120 write (9,12) nlbcd(m), (vpar(7,m,n), n=1,n2)
   12 format (1h ,a3,f12.5,11f9.5)
      if (ncspvs.le.12) go to 200
      write (9,33) (nlbcd(n), n=13,ncspvs)
      write (9,35)
      do 140 m=13,ncspvs
  140 write (9,12) nlbcd(m), (vpar(7,m,n), n=13,m)
c
c          calc coefs of slater parameters, and interaction energies
c
  200 do 250 m=1,ncspvs
      b(m)=0.0
      flm=l(m)
      if (wwnl(m).gt.1.0) go to 202
      emn(m,m)=0.0
      go to 220
  202 emn(m,m)=vpar(1,m,m)
      if (l(m).le.0) go to 220
      im=l(m)+1
      fk=0.0
      do 206 i=2,im
      fk=fk+2.0
      coef=-s3j0sq(flm,fk,flm)*(2.0*flm+1.0)/(4.0*flm+1.0)
  206 emn(m,m)=emn(m,m)+coef*vpar(i,m,m)
c
  220 n1=m+1
      if (n1.gt.ncspvs) go to 250
      do 240 n=n1,ncspvs
      fln=l(n)
      k=iabs(l(m)-l(n))-2
      fk=k
      emn(m,n)=vpar(1,m,n)
      im=(l(m)+l(n)-k)/2+1
      do 225 i=2,im
      fk=fk+2.0
      coef=-s3j0sq(flm,fk,fln)/2.0
  225 emn(m,n)=emn(m,n)+coef*vpar(i,m,n)
      emn(n,m)=emn(m,n)
  240 continue
  250 continue
c
      n=ncspvs
      if (ee(n).lt.0.0) go to 300
      do 265 m=1,n
      emn(m,n)=0.0
  265 emn(n,m)=0.0
      evel(n)=0.0
      edar(n)=0.0
      ecorr(n)=0.0
      ucorr(n)=0.0
      epss(n)=0.0
c
  300 if (iptvu.eq.0) go to 400
      write (9,30) conf,nconf,iz,ncores,nvales,ion,c,ca1
   30 format (3h1  ,3a6,2x,10h    nconf=,i3,6h    z=,i3,11h    ncores=,
     1  i2,11h    nvales=,i2, 8h    ion=,i2, 14h    r(i)/x(i)=,f10.8,
     2  6x,4hca1=,f5.3)
      write (9,31) tolstb,tolkm2,tolend,thresh,kut, exf,corrf,ca0,ihf1,
     1  mesh,idb,rdb,emx,r(mesh),irel
   31 format (/8h0tolstb=,f5.3,3x,7htolk-2=,f7.5,3x,7htolend=,e8.1,
     1  3x,7hthresh=,e8.1,3x,4hkut=,i2,4x,4hexf=,f5.3,9h   corrf=,f5.3,
     2  6x,4hca0=,f6.3//18h0rcn mod 36  ihf1=,i2,
     3  10x,5hmesh=,i4,3x,4hidb=,i4,3x,4hrdb=,f8.3,3x,4hemx=,f7.2,
     4  3x,8hr(mesh)=,f8.3,5x,5hirel=,i2)
      if (iptvu.eq.1) go to 400
      write (9,32)
   32 format (21h0interaction energies)
      n2=min0(ncspvs,12)
      write (9,33) (nlbcd(n), n=1,n2)
   33 format (///5h     ,12(6x,a3))
      write (9,35)
      do 320 m=1,ncspvs
      n2=min0(m,12)
  320 write (9,34) nlbcd(m), (emn(m,n), n=1,n2)
   34 format (1h ,a3,f12.4,11f9.4)
      if (ncspvs.le.12) go to 400
      write (9,33) (nlbcd(n), n=13,ncspvs)
      write (9,35)
   35 format (1h )
      do 340 m=13,ncspvs
  340 write (9,34) nlbcd(m), (emn(m,n), n=13,m)
      if (ncspvs.le.12) go to 400
      write (9,30) conf,nconf,iz,ncores,nvales,ion,c,ca1
      write (9,31) tolstb,tolkm2,tolend,thresh,kut, exf,corrf,ca0,ihf1,
     1  mesh,idb,rdb,emx,r(mesh),irel
c
c          calc ionization potentials, eps(m), and total correlation ene
c
  400 write (9,40)
   40 format (  //1h0,15x,19hwithout correlation,16x,16hwith correlation
     1  /14x,23h-----------------------,11x,23h-----------------------/
     2  '  nl   wnl ',6x,6heps fg,6x,7heps fgr,4x,4hcorr,7x,6heps fg,6x,
     3  7heps fgr,5x,3hn*r,5x,4hn*rc,5x,3hect,5x,19hec-----------------,
     4  /)
      ect=0.0
      xi(1)=0.0
      fion=ion+1
      do 405 i=1,mesh
  405 xj(i)=0.0
      do 460 m=1,ncspvs
  409 a0=0.0
      do 412 np=1,ncspvs
      n=ncspvs+1-np
      eras=wwnl(n)
      if (m.ne.n) go to 411
      eras=eras-1.0
  411 a0=a0+eras*emn(m,n)
  412 a(n,1)=a0
      eps2=ei(m)+a0
      evd=0.0
      if (irel.ne.0) go to 415
      evd=evel(m)
      if (l(m).eq.0) evd=evd+edar(m)
  415 epscor=ecorr(m)
      eps1c=epss(m)
      eps1=eps1c-epscor
      eps2c=eps2+epscor
      eps2r=eps2+evd
      eps2cr=eps2c+evd
      ehf(m)=abs(eps2)
      fnstr=0.0
      fnstrc=0.0
      if (eps2.ge.0.0) go to 420
      fnstr=fion/sqrt(-eps2r)
      fnstrc=fion/sqrt(-eps2cr)
  420 k=wwnl(m)+0.999999
      kkk=nkkk(m)
      do 426 i=1,kkk
      xj(i)=0.0
      do 424 mm=1,ncspvs
      if (mm.eq.m) go to 424
      if (ee(mm).ge.ee(m)) go to 424
      xj(i)=xj(i)+wwnl(mm)*(pnl(i,mm)**2)
  424 continue
  426 continue
      do 440 lp=1,k
      ax=0.0
      if (ee(m).ge.0.0) go to 439
      do 435 i=2,kkk
      eras2=pnl(i,m)**2
      eras1=r(i)
      eras=xj(i)
      call subcor
      xi(i)=b1*eras2/r(i)
  435 xj(i)=xj(i)+eras2
      call quad5(xi,1,1,kkk,ax)
  439 kk=min0(lp,4)
      ecm(kk)=ax
      occ=1.0
      if (lp.eq.k) occ=mod(wwnl(m),1.0)
      if (occ.lt.0.000001) occ=1.0
      ect=ect+occ*ax
  440 continue
      if (ifrac.eq.0) then
        write (9,44) nlbcd(m),wwnl(m),eps2,eps2r,epscor,eps2c,eps2cr,
     1  fnstr,fnstrc, ect, (ecm(k), k=1,kk)
   44 format (1x,a3,f6.0,1x,2f13.5,f8.5,2f13.5,2f8.4,2f10.5,4f7.4)
      else
        write (9,45) nlbcd(m),wwnl(m),eps2,eps2r,epscor,eps2c,eps2cr,
     1  fnstr,fnstrc, ect, (ecm(k), k=1,kk)
   45 format (1x,a3,f7.3,2f13.5,f8.5,2f13.5,2f8.4,2f10.5,4f7.4)
      end if
      etvd=etotrl
      do 450 mp=1,ncspvs
      eras=a(mp,1)*wwnl(m)
      if (m.lt.mp) go to 450
      eras=0.5*eras
  450 b(mp)=b(mp)+eras
      eeo(m)=ee(m)
      if (m.ge.n1sc) ee(m)=eps2cr
  460 continue
      etvd=etotrl
      emag=0.0
      eret=0.0
      if (irel.ge.2) call ebreit
c
c          calc total binding energy
c
  500 a0=0.0
      b0=0.0
      eras=0.0
      do 501 mp=1,ncspvs
      m=ncspvs+1-mp
      temp=1.0
      if (irel.ne.0) temp=0.0
      eras=eras+temp*(evel(m)+edar(m))*wwnl(m)
      erel(m)=eras
      a0=a0+ei(m)*wwnl(m)
      a(m,1)=a0+b(m)
      b0=b0+ecorr(m)*wwnl(m)
  501 b(m)=b0
      et2=a(1,1)+eterm
      et1c=etot+eterm
      et1=et1c-ecorrt
      et2c=et2+ecorrt
      et3c=et2+ect
      et2r=et2+etvd
      et2cr=et2c+etvd
      et3cr=et3c+etvd
      et4cr=et3cr+emag+eret
      write (9,49) ecorrt,ect,etotrl, emag,eret
   49 format (/1h0,27x,7hecorrt=,f10.5,6h (old),14x,4hect=,f10.5,
     1  6h (new),5x,7hetotrl=,f13.5,2h +,f9.5,2h +,f9.5)
      write (9,50)    et2,et2r,et2c,et2cr,et3c,et3cr,et4cr
   50 format( /'      etot=',2f13.5,f21.5,4f13.5)
c          store eav for rcg input
c
      vprm(1)=et4cr-eterm
c
c
c          print rcg input parameter values
c
  550 if (nzi.le.0) go to 560
      do 552 i=1,nzi
      npar=npar+1
      kpar(npar)=2
  552 vprm(npar)=zi(i)
  560 if (nfij.le.0) go to 570
      do 562 i=1,nfij
      npar=npar+1
      kpar(npar)=3
  562 vprm(npar)=fij(i)
  570 if (ngij.le.0) go to 580
      do 572 i=1,ngij
      npar=npar+1
      kpar(npar)=4
  572 vprm(npar)=gij(i)
c
  580 write (9,57) conf,npar, (vprm(i),kpar(i), i=1,npar)
   57 format (//1h ,3a6,i6,f12.4,i2,4(f10.4,i2)/1x,7(f10.4,i2)/
     1  1x,7(f10.4,i2)/1x,7(f10.4,i2)/1x,7(f10.4,i2))
      if (iw6.lt.0) write (6,58) npar, (vprm(i),kpar(i), i=1,npar)
   58 format (' Eav & pars',i6,f12.4,i2,4(f8.4,i2)/1x,7(f8.4,i2)/
     1  1x,7(f8.4,i2)/1x,7(f8.4,i2)/1x,7(f8.4,i2))
c
  600 return
      end
      function s3j0sq(fj1,fj2,fj3)
c
c          calc square of 3-j symbol with zero magnetic quantum numbers
c
      implicit real*8 (a-h,o-z)
c
      fj=fj1+fj2+fj3
      a=fj-fj1-fj1
      b=fj-fj2-fj2
      c=fj-fj3-fj3
      s3j0sq=fctrl(a)*fctrl(b)*fctrl(c)/fctrl(fj+1.0)
      a=fctrl(a/2.0)*fctrl(b/2.0)*fctrl(c/2.0)/fctrl(fj/2.0)
      s3j0sq=s3j0sq/a/a
      return
      end
      function fctrl(a)
c
c          calculate factorial of a
c
      implicit real*8 (a-h,o-z)
c
      fctrl=1.0
      if (abs(a)-0.1) 10,10,11
   10 a=0.0
   11 if (a) 12,15,13
   12 fctrl=0.0
      go to 15
   13 imax=a+0.1
      do 14 i=1,imax
   14 fctrl=fctrl*float(i)
   15 return
      end
      subroutine hfwrtp
c
c     subroutine written by d.c.griffin to convert hx output to hf input
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c5/izo,iono,lbcd(21),morb(7),eeo(ko),
     1  ee8,noelec,ierror,wwnlo(ko),n1sc,n1scf,iblank
      common/c6/a0,a1,a2,a3,a4,ab1,ab3,b0,b1,b2,b3,eras,eras1,eras2
      common/c7/t0,t1,t2
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
c
c     equivalence (ru,      f1), (ruee,  f2), (ruexch,rsq,xa,emn),
c    1 (rsatom,qq,drudr,xb), (a,delp,rryd), (b,alfm,rkay)
c     equivalence (rsatom,s1),(rsatom(251),s2),(rsatom(501),zeta)
c
      integer fg
      dimension maxi(ko),rhf(200),phf(ko,200)
c     equivalence (ruexch,maxi),(xi,rhf),(v,phf)
c
c
c
  200 rhohf=-3.0
      h=0.0625
      no=200
      rhoi=rhohf
      do 250 i=1,no
      rhf(i)=exp(rhoi)/z
  250 rhoi=rhoi+h
      notemp=no
      if (rhf(200).gt.r(mesh)) notemp=log(r(mesh)/rhf(1))/h
      nufsh=0
      do 270 i=1,ncspvs
      if(wwnl(i).eq.(4.*float(l(i))+2.)) go to 270
      nufsh=nufsh+1
  270 continue
      write (7) conf,nconf,iz,z,ncspvs,rhohf,h,no,nhftrm,irel,etotrl,
     1  nufsh,izhxbw,ns,iphfwf
c
      if (nhftrm.eq.0) go to 502
      do 501 i=1,nhftrm
      nfg=nf(i)+ng(i)
  501 write (7) t2hf(i),nf(i),ng(i),nfg,
     1  (cfg(j,i),fg(j,i),kfg(j,i),ifg(j,i),jfg(j,i),j=1,nfg)
c
  502 do 503 i=1,ncspvs
      j=imat(i)
      rhfmtc(i)=r(j)
      msend=nkkk(i)
      rend=r(msend)
      max=log(rend/rhf(1))/h
      if(max.gt.notemp) max=notemp
      maxi(i)=max
  503 write (7) nlbcd(i),wwnl(i),az(i),ehf(i),rhfmtc(i),max
c
      jmin=2
      do 510 i=1,no
      do 504 m=1,ncspvs
  504 phf(m,i)=0.0
      do 508 j=jmin,mesh
      if (rhf(i).gt.r(j)) go to 508
      jm2=j-2
      jm1=j-1
      if (mod(jm2,40).eq.0.and.jm2.lt.idb) jm2=j-3
      u=(rhf(i)-r(jm2))/(r(jm1)-r(jm2))
      u2=0.5*u*(u-1.0)
      rhm1=1.0/sqrt(rhf(i))
      do 506 m=1,ncspvs
      if (i.gt.maxi(m).or.nnn(m).gt.90) go to 506
      del1=pnl(jm1,m)-pnl(jm2,m)
      del2=pnl(j,m)-pnl(jm1,m)
      phf(m,i)=(pnl(jm2,m)+del1*u+(del2-del1)*u2)*rhm1
  506 continue
      go to 509
  508 continue
  509 jmin=j-1
  510 continue
  520 write (7) (rhf(i),i=1,no), ((phf(i,j),j=1,no), i=1,ncspvs)
c
      meshi=min0(mesh,1001)
      if (nnn(ncspvs).gt.10) write (7) meshi,idb,
     1  (r(j),pnl(j,ncspvs), j=1,meshi)
c
      if (ihf.lt.3) go to 700
      idel=1
      if (iphfwf.lt.0) idel=4
      write(9,103)
  103 format (1h1,50x,22hhf input wavefunctions)
      npwftp=ncspvs-10
      if (npwftp.lt.1) npwftp=1
      write(9,104) (nlbcd(i),i=npwftp,ncspvs)
  104 format (13x,3hrhf,6x,a3,10(7x,a3))
      do 14 i=1,no,idel
      eras=sqrt(rhf(i))
      do 15 j=npwftp,ncspvs
      phf(j,i)=phf(j,i)*eras
   15 continue
   14 continue
      do 13 i=1,no,idel
      write (9,102) i,rhf(i),(phf(j,i),j=npwftp,ncspvs)
  102 format (i5,f12.5,11f10.5)
   13 continue
  700 return
      end
      subroutine hfpot(m)
c
c          calculate hartree-fock potential terms over and above
c               the classical potential
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/c8/nhftrm,nfg,etotrl,ehf(ko),t1hf(9),t2hf(9),eterm,
     1  nf(9),ng(9),cfg(10,9),fg(10,9),kfg(10,9),rhfmtc(ko),ii
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      common/lc2/ra(kmsh),rb(kmsh)
c
      integer fg
      data jhlf,jhlg/102,103/
      data jhuF,jhuG/70,71/
c
      kkk1=min0(nkkk(m)+40,mesh)
      kkk=kkk1
      x1(1)=0.0
      x2(1)=0.0
      do 110 i=1,mesh
      ra(i)=0.0
  110 rb(i)=0.0
      if (niter.eq.1.and.ee(m).gt.0.0) go to 455
      if ((wwnl(m).lt.2.0.and.nhftrm.eq.0).or.l(m).eq.0) go to 300
c          calc direct terms
      fl=l(m)
      c1=(wwnl(m)-1.0)*(4.0*fl+2.0)/(4.0*fl+1.0)
      do 149 m1=1,ncspvs
      if (ee(m1).gt.0.0) go to 149
      kx=2*l(m1)
      if (kx.lt.2) go to 149
      full=2*kx+2
      kkk=nkkk(m1)
      do 120 i=1,kkk
  120 xj(i)=pnl(i,m1)**2
      do 140 k=2,kx,2
      fk=k
      c2=0.0
      if (nhftrm.eq.0.or.wwnl(m1).eq.full) go to 125
      if (nfg.le.0) go to 125
c               ls-dependent terms
      do 123 j=1,nfg
      if ((fg(j,ii).ne.jhlf.and.fg(j,ii).ne.jhuF).or.kfg(j,ii).ne.k)
     1  go to 123
      if ((ifg(j,ii).eq.nlbcd(m1).and.jfg(j,ii).eq.nlbcd(m)).or
     1   .(jfg(j,ii).eq.nlbcd(m1).and.ifg(j,ii).eq.nlbcd(m)))
     2    c2=-2.0*cfg(j,ii)/wwnl(m)
  123 continue
c
  125 if (m1.ne.m) go to 128
      c2=2.0*c2
      c2=c2+c1*s3j0sq(fl,fk,fl)
  128 continue
      if (c2.eq.0.0) go to 140
      call quadk(m,1)
      do 130 i=3,kkk,2
      eras=r(i)**k
      eras=x1(i)/eras+eras*r(i)*(x2(kkk)-x2(i))
  130 ra(i)=ra(i)+c2*eras
  140 continue
  149 continue
      kkk=kkk1
      do 150 i=3,kkk,4
      b0=(ra(i)-ra(i-2)+ra(i)-ra(i+2))*0.125
      ra(i-1)=0.5*(ra(i)+ra(i-2))+b0
  150 ra(i+1)=0.5*(ra(i)+ra(i+2))+b0
      if (ra(5).eq.0.0.or.ra(9).eq.0.0) go to 300
      ra(2)=(ra(3)**2)/ra(5)
      ra(4)=ra(5)*ra(7)/ra(9)
c          calc exchange terms
  300 if (ncspvs.eq.1) go to 455
      fl2=l(m)
      do 400 m1=1,ncspvs
      if (m1.eq.m) go to 400
      if (ee(m1).gt.0.0) go to 400
      full=4*l(m1)+2
      kkk=min0(nkkk(m1),nkkk(m))
      fl1=l(m1)
      do 320 i=1,kkk
  320 xj(i)=pnl(i,m1)*pnl(i,m)
      kn=iabs(l(m1)-l(m))
      kx=l(m1)+l(m)
      do 340 k=kn,kx,2
      fk=k
      c2=0.0
      if (nhftrm.eq.0.or.wwnl(m1).eq.full) go to 325
      if (nfg.le.0) go to 325
c               ls-dependent terms
      do 323 j=1,nfg
      if ((fg(j,ii).ne.jhlg.and.fg(j,ii).ne.jhuG).or.kfg(j,ii).ne.k) 
     1  go to 323
      if ((ifg(j,ii).eq.nlbcd(m1).and.jfg(j,ii).eq.nlbcd(m)).or
     1   .(jfg(j,ii).eq.nlbcd(m1).and.ifg(j,ii).eq.nlbcd(m)))
     2    c2=-2.0*cfg(j,ii)/wwnl(m)
  323 continue
c
  325 c2=c2+wwnl(m1)*s3j0sq(fl1,fk,fl2)
      call quadk(m,1)
      do 330 i=3,kkk,2
      eras=r(i)**k
      eras=x1(i)/eras+eras*r(i)*(x2(kkk)-x2(i))
  330 rb(i)=rb(i)+c2*eras*pnl(i,m1)
  340 continue
  400 continue
      kkk=kkk1
      do 450 i=3,kkk,4
      b0=(rb(i)-rb(i-2)+rb(i)-rb(i+2))*0.125
      rb(i-1)=0.5*(rb(i)+rb(i-2))+b0
  450 rb(i+1)=0.5*(rb(i)+rb(i+2))+b0
      if (rb(5).eq.0.0.or.rb(9).eq.0.0) go to 455
      rb(2)=(rb(3)**2)/rb(5)
      rb(4)=rb(5)*rb(7)/rb(9)
c
  455 frac=0.01
      pmax=0.0
      do 458 i=2,kkk
      eras=pnl(i,m)
      pmax=max(pmax,abs(frac*eras))
      if (eras.ge.0.0.and.eras.lt.pmax) eras=pmax
      if (eras.lt.0.0.and.eras.gt.-pmax) eras=-pmax
  458 rb(i)=rb(i)/eras
      idel=1
      idel2=idel*2
      km10=kkk-10
      do 465 i=10,km10
      if (pnl(i,m)*pnl(i+1,m).gt.0.0) go to 465
      eras=(rb(i+1+idel)-rb(i-idel))/(r(i+1+idel)-r(i-idel))
      do 462 j=1,idel2
  462 rb(i-idel+j)=rb(i-idel)+eras*(r(i-idel+j)-r(i-idel))
  465 continue
      do 468 i=2,kkk
  468 xj(i)=ru(i)-xi(i)-exf10*ruexch(i)-ra(i)-rb(i)
      if (kkk.ge.mesh) go to 480
      do 470 i=kkk,mesh
  470 xj(i)=xj(kkk)
  480 continue
      if (end.gt.0.0.and.iptvu.ge.5) write (9,46) m,k,kkk,
     1  (ra(i),i=2,kkk),(rb(i),i=2,kkk),(xj(i),i=2,kkk)
   46 format (1h0/3i10/(10f12.4))
      return
      end
      subroutine quadk(m,kk)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c3/npar,nzi,nfij,ngij,vprm(50),kpar(50),zi(8),fij(30),
     1  gij(35),norb,norbpt,dh(5),delta,npr
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/lc2/ra(kmsh),rb(kmsh)
c
c       calculate  const*yk function   (integrate by simpsons rule)
c
  510 do 511 i=2,kkk
      eras=r(i)**k
      x1(i)=eras*xj(i)
  511 x2(i)=xj(i)/(eras*(r(i)**kk))
  512 x1(1)=0.0
      x2(1)=2.0*x2(2)-x2(3)
      if (abs(x2(1)).lt.0.2*abs(x2(2))) x2(1)=0.0
c
      ho3=r(2)/3.0
      i=1
      b12=x1(1)
      b22=x2(1)
      x2(1)=0.0
      do 515 j=1,nblock
      do 513 k1=1,20
      i=i+2
      if (i.gt.kkk) go to 520
      b10=b12
      b12=x1(i)
      x1(i)=x1(i-2)+ho3*(b10+4.0*x1(i-1)+b12)
      b20=b22
      b22=x2(i)
  513 x2(i)=x2(i-2)+ho3*(b20+4.0*x2(i-1)+b22)
      if (i.gt.idb) go to 515
      ho3=ho3+ho3
  515 continue
  520 kkk1=kkk
      if (m.le.25) kkk=min0(nkkk(m)+40,mesh)
      if (m.gt.25) kkk=min0(kkk+40,mesh)
      do 527 i=kkk1,kkk,2
      x1(i)=x1(kkk1)
  527 x2(i)=x2(kkk1)
      kkk=kkk1
  530 n=2
      if (m.lt.31) go to 600
c
      do 550 i=3,kkk,4
      b0=(x1(i)-x1(i-2)+x1(i)-x1(i+2))*0.125
      x1(i-1)=0.5*(x1(i)+x1(i-2))+b0
      x1(i+1)=0.5*(x1(i)+x1(i+2))+b0
      b0=(x2(i)-x2(i-2)+x2(i)-x2(i+2))*0.125
      x2(i-1)=0.5*(x2(i)+x2(i-2))+b0
  550 x2(i+1)=0.5*(x2(i)+x2(i+2))+b0
      n=1
  600 if (end.gt.0.0.and.iptvu.ge.7) write (9,60) m,k,kkk1,kkk,
     1  (xj(i),i=1,kkk),(x1(i),i=1,kkk,n),(x2(i),i=1,kkk,n)
   60 format (1h0/4i10/(10f12.4))
      return
      end
      function zetabw(m,rm3)
c
c      calculate spin-orbit parameter via blume-watson method
c       proc. roy. soc. (london) a270, 127 (1962); a271, 565 (1963)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      character*8 nlbcd(ko),nlbcdo(ko),nlbcd8(ko),ifg(10,9),jfg(10,9),
     1  nla,nlb
      character conf(3)*6,inpcard*70
      common/char/conf,nlbcd,nlbcdo,nlbcd8,ifg,jfg,nla,nlb,inpcard
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,dummy(ko),iw6,nconft,ifrac
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      dimension snkk(ko),vkk(ko)
c
      sp=0.0
      if (l(m).eq.0) go to 400
      flm=l(m)
      e2=-12*l(m)-6
      e3=3.0*sqrt((2.0*flm+1.0)/(flm*(flm+1.0)))
      sp=z*rm3*5.843660
      if (ihf1.gt.5) write (9,8810) m,flm,e2,e3,sp
 8810 format (/i10,4f15.6)
c
      do 290 mp=1,ncspvs
      if (ee(mp).gt.0.0) go to 290
      if (mp.eq.m) go to 290
      lp=l(mp)
      flp=lp
      wp=wwnl(mp)
      k=0
      sp=sp-2.0*wp*sm(m,mp)
c
      kn=iabs(l(m)-lp)-2
      kx=l(m)+lp
      inv=-1
      do 130 k=kn,kx
      snkk(k+3)=0.0
      vkk(k+3)=0.0
      inv=-inv
      if (inv.lt.0) go to 120
      if (k.lt.-1) go to 130
      if (k.eq.kx.and.l(m).eq.lp) go to 130
      snkk(k+3)=sn(m,mp)
      go to 130
  120 if (k.lt.0) go to 130
      vkk(k+3)=vk(m,mp)
  130 continue
      if (ihf1.gt.5) write (9,8830) m,mp,nlbcd(m),nlbcd(mp),flp,wp,sp,
     1  (k,snkk(k),vkk(k), k=1,9)
 8830 format (/2i10,3x,a4,3x,a4,3f15.6,/(i10,2f30.6))
      kn=kn+2
      kxm2=kx-2
      if (kxm2.lt.kn) go to 200
      sum2=0.0
      do 160 k=kn,kxm2
      fk=k
      fkp1=fk+1.0
      eras=s6j(fk,1.0d0,fkp1,flm,flp,flm)**2
      eras=eras*s3j0sq(flm,fk,flp)
      eras=eras*(2.0*fk+1.0)*(2.0*fk+3.0)/(fk+2.0)
      sum2=sum2+eras*snkk(k+3)
      if (ihf1.gt.5) write (9,8840) k,fk,eras,snkk(k+3),sum2
 8840 format (/i10,5f15.6)
  160 continue
      sp=sp-e2*wp*sum2
c
  200 e13=flp*(flp+1.0)-flm*(flm+1.0)
      if (ihf1.gt.5) write (9,8850) sp,e13
 8850 format (/4f15.6)
      if (kn.eq.0) kn=2
      sum13=0.0
      do 220 k=kn,kx
      fk=k
      eras=sqrt(fk*(fk+1.0)*(2.0*fk+1.0))
      eras=eras*s3j0sq(flm,fk,flp)
      eras=eras*s6j(fk,1.0d0,fk,flm,flp,flm)
      eras1=e13*(snkk(k+1)/fk-snkk(k+3)/(fk+1.0))
      sum13=sum13+eras*(vkk(k+2)+eras1)
      if (ihf1.gt.5) write (9,8840) k,fk,eras,eras1,sum13,vkk(k+2)
  220 continue
      sp=sp-e3*wp*sum13
      if (ihf1.gt.5) write (9,8850) sp,e3,wp,sum13
  290 continue
c
  300 w=wwnl(m)
      if (w.lt.2.0) go to 400
      k=0
      sm0=sm(m,m)
      sp=sp-(2.0*w-3.0)*sm0
      if (ihf1.gt.5) write (9,8850) sp,sm0
      if (flm.lt.2.0) go to 400
      e=6.0*((2.0*flm+1.0)**2)
      kn=2
      kx=2*l(m)-2
      do 330 k=kn,kx
      fk=k
      fkp1=fk+1.0
      sum=0.0
      fkp=fk-2.0
      do 320 kp=1,2
      fkp=fkp+2.0
      eras=s6j(fkp1,1.0d0,fkp,flm,flm,flm)**2
      eras=eras*s3j0sq(flm,fkp,flm)
  320 sum=sum+(fkp1-fkp)*(2.0*fkp+1.0)*eras
      smk=sm(m,m)
      sum1=sum*(2.0*fk+3.0)*smk
      sp=sp+e*sum1
      if (ihf1.gt.5) write (9,8840) k,eras,sum,smk,sum1,sp
  330 continue
c
  400 zetabw=sp
      return
      end
      function sm(m,mp)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
c
      kkk=min0(nkkk(m),nkkk(mp))
      call zk(m,mp,m,mp)
      call quad5(xj,1,1,kkk,a1)
      sm=a1*2.921830
      return
      end
      function sn(m,mp)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
c
      kkk=min0(nkkk(m),nkkk(mp))
      call zk(m,mp,mp,m)
      call quad5(xj,1,1,kkk,a1)
      sn=a1*2.921830
      return
      end
      subroutine zk(m1,m2,m3,m4)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      common/lc2/ra(kmsh),rb(kmsh)
c
      do 120 i=1,kkk
  120 xj(i)=pnl(i,m2)*pnl(i,m4)
      call quadk(m1,3)
      kp3=k+3
      do 160 i=3,kkk,2
  160 xj(i)=pnl(i,m1)*pnl(i,m3)*x1(i)/(r(i)**kp3)
      xj(1)=0.0
      do 170 i=3,kkk,4
      b0=(xj(i)-xj(i-2)+xj(i)-xj(i+2))*0.125
      xj(i-1)=0.5*(xj(i)+xj(i-2))+b0
  170 xj(i+1)=0.5*(xj(i)+xj(i+2))+b0
      xj(2)=0.0
      if (xj(5).ne.0.0) xj(2)=(xj(3)**2)/xj(5)
      do 180 i=4,14,2
      xm=xj(i-1)
      x0=xj(i)
      xp=xj(i+1)
      if (x0.lt.xm.and.x0.lt.xp) go to 175
      if (x0.gt.xm.and.x0.gt.xp) go to 175
      if (xm.eq.0.0) go to 180
      if (xp/xm.lt.6.0) go to 180
  175 xj(i)=sqrt(xm*xp)
      if (xp.lt.0.0) xj(i)=-xj(i)
  180 continue
      return
      end
      function vk(m,mp)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
c
      kkk=min0(nkkk(m),nkkk(mp))
      call dyk(m,mp)
      call quad5(xj,1,1,kkk,a1)
      vk=a1
      call dyk(mp,m)
      call quad5(xj,1,1,kkk,a1)
      vk=(vk-a1)*2.921830
      return
      end
      subroutine dyk(m,mp)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      common/lc2/ra(kmsh),rb(kmsh)
c
      thm1=1.0/r(3)
      i=0
      eras=0.0
      eras1=pnl(1,mp)
      eras2=eras1
      do 140 j=1,nblock
      if (i.le.idb) eras1=eras
      do 130 jj=1,40
      i=i+1
      ip1=i+1
      if (ip1.gt.kkk) ip1=kkk
      eras=eras1
      eras1=eras2
      eras2=pnl(ip1,mp)
  130 xj(i)=pnl(i,m)*(thm1*(eras2-eras)-eras1/r(i))
      if (i.gt.idb) go to 140
      thm1=0.5*thm1
  140 continue
      call quadk(m,3)
      kp1=k+1
      do 160 i=3,kkk,2
      eras=r(i)**kp1
      eras=x1(i)/(eras*r(i))+eras*(x2(kkk)-x2(i))
  160 xj(i)=pnl(i,m)*pnl(i,mp)*eras
      xj(1)=0.0
      do 170 i=3,kkk,4
      b0=(xj(i)-xj(i-2)+xj(i)-xj(i+2))*0.125
      xj(i-1)=0.5*(xj(i)+xj(i-2))+b0
  170 xj(i+1)=0.5*(xj(i)+xj(i+2))+b0
      xj(2)=0.0
      if (xj(5).ne.0.0) xj(2)=(xj(3)**2)/xj(5)
      do 180 i=4,14,2
      xm=xj(i-1)
      x0=xj(i)
      xp=xj(i+1)
      if (x0.lt.xm.and.x0.lt.xp) go to 175
      if (x0.gt.xm.and.x0.gt.xp) go to 175
      if (xm.eq.0.0) go to 180
      if (xp/xm.lt.6.0) go to 180
  175 if (xm*xp.lt.0.0) go to 178
      xj(i)=sqrt(xm*xp)
      if (xp.lt.0.0) xj(i)=-xj(i)
      go to 180
  178 xj(i)=0.5*(xm+xp)
  180 continue
      return
      end
      function s6j(fj1,fj2,fj3,fl1,fl2,fl3)
c
c          calculate 6-j symbol
c
      implicit real*8 (a-h,o-z)
c
      a6j=0.0
      f1=delsq(fj1, fj2, fj3)*delsq(fj1, fl2, fl3)
     1*delsq(fl1, fj2, fl3)*delsq(fl1, fl2, fj3)
      if (f1.le.0.0) go to 214
  210 a=fj1+fj2+fj3
      b=fj1+fl2+fl3
      c=fl1+fj2+fl3
      d=fl1+fl2+fj3
      e=fj1+fj2+fl1+fl2
      f=fj2+fj3+fl2+fl3
      g=fj3+fj1+fl3+fl1
      fi=max(a,b,c,d)
      i1=fi+0.001
      i2=min(e,f,g)+0.001
      if (i2.lt.i1) go to 214
  211 do 212 i=i1,i2
      f2=fctrl(fi-a)*fctrl(fi-b)*fctrl(fi-c)*fctrl(fi-d)
      f3=fctrl(e-fi)*fctrl(f-fi)*fctrl(g-fi)
      a6j=a6j+((-1.0)**mod(i,2))*fctrl(fi+1.0)/(f2*f3)
  212 fi=fi+1.0
      s6j=a6j*sqrt(f1)
      return
  214 s6j=a6j
  215 return
      end
      function delsq(a,b,c)
c
c          short (but slower) version calling factorial function
c
      implicit real*8 (a-h,o-z)
c
      d=a+b
      e=d+c
c     if (aint(e).ne.e) go to 120
      tenm6=1.0e-6
      if (abs(aint(e)-e).gt.tenm6) go to 120
      f=a-b
      delsq=fctrl(d-c)*fctrl(c+f)*fctrl(c-f)/fctrl(e+1.0)
      return
  120 delsq=0.0
      return
      end
      subroutine ebreit
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/pq/pnl(kmsh,ko),qnl(kmshq,koq)
      common/brt/emag,eret,fmi(10,3),fni(10,4),fnorm(ko)
      dimension uab(kmsh),vab(kmsh)
      equivalence (uab,rscore),(vab,rsvale)
c
c          calculate normalization constants
c
      mmax=ncspvs
      if (nnn(ncspvs).gt.98) mmax=ncspvs-1
      do 110 m=1,mmax
      kkk=nkkk(m)
      do 105 i=1,kkk
      xi(i)=pnl(i,m)**2+qnl(i,m)**2
      call quad5(xi,1,1,kkk,a1)
  105 fnorm(m)=a1
  110 continue
c
      emag=0.0
      eret=0.0
c
c          sum over all pairs of orbitals
c
      do 690 mb=1,mmax
      if (l(mb).gt.4) go to 690
      do 680 ma=1,mb
      if (l(ma).gt.4) go to 680
      if (ma.eq.mb.and.wwnl(ma).eq.1.0) go to 680
      kkk=min0(nkkk(ma),nkkk(mb))
      do 120 i=1,kkk
      eras=pnl(i,ma)*qnl(i,mb)
      eras1=pnl(i,mb)*qnl(i,ma)
      uab(i)=eras-eras1
  120 vab(i)=eras+eras1
c
c          calculate mi and ni integrals
c
      kp1x=l(ma)+l(mb)+2
      do 150 kp1=1,kp1x
      do 125 i=1,4
      if (i.eq.4) go to 125
      fmi(kp1,i)=0.0
  125 fni(kp1,i)=0.0
      k=kp1-1
      if (ma.eq.mb) go to 135
      do 130 i=1,kkk
  130 xj(i)=uab(i)
      call quadk(31,1)
      call brtint(1,ma,mb)
  135 do 140 i=1,kkk
  140 xj(i)=vab(i)
      call quadk(32,1)
  150 call brtint(2,ma,mb)
c
c          sum over both kappas, each orbital
c
      fkapb=l(mb)
      if (fkapb.eq.0.0) fkapb=-1.0
  200 eras=2*l(mb)+1
      qb=wwnl(mb)*abs(fkapb)/eras
      fjb=abs(fkapb)-0.5
      fkapa=l(ma)
      if (fkapa.eq.0.0) fkapa=-1.0
  210 eras=2*l(ma)+1
      qa=wwnl(ma)*abs(fkapa)/eras
      fja=abs(fkapa)-0.5
      kp1n0=abs(fja-fjb)+1.0
      kp1x=fja+fjb+1.0
      flmag=0.0
      flret=0.0
c
c          calculate lambda(kapa,kapb) terms
c
  300 kp1n=kp1n0
      if (mod(l(ma)+kp1n0+l(mb),2).eq.0) kp1n=kp1n0+1
      if (ma.eq.mb.and.fkapa.eq.fkapb) go to 400
      eras=fkapb-fkapa
      do 320 kp1=kp1n,kp1x,2
      k=kp1-1
      fk=k
      fm1=0.5
      fm2=0.0
      flambda=s3j(fja,fk,fjb,fm1,fm2)**2
      t1=fk*fmi(k,1)
      if (eras.ne.0.0) t1=t1+eras*(2.0*fmi(k,2)+eras*fmi(k,3)/fk)
      t1=t1/(2.0*fk-1.0)
      fkp1=fk+1.0
      t2=fkp1*fmi(k+2,1)
      if (eras.ne.0.0) t2=t2+eras*(-2.0*fmi(k+2,2)+eras*fmi(k+2,3)/fkp1)
      t2=t2/(2.0*fk+3.0)
      x=fk*fkp1*(fni(k+2,1)-fni(k,1))
      if (eras.eq.0.0) go to 310
      x=x+eras*(fkp1*(fni(k+2,2)-fni(k,2))+fk*(fni(k,3)-fni(k+2,3)))
      x=x+(eras**2)*(fni(k,4)-fni(k+2,4))
  310 flmag=flmag+flambda*(t1+t2)
      flret=flret+flambda*(fk*t1+fkp1*t2+x)/(fk+fkp1)
      if (itpow.eq.4)
     1  write (9,3100) ma,mb,fkapa,fkapb,k,flambda,t1,t2,x,flmag,flret
 3100 format (2i4,2f5.1,i4,6f12.7)
  320 continue
c
c          calculate lambda(-kapa,kapb) terms
c
  400 eras=fkapa+fkapb
      if (eras.eq.0.0) go to 500
      kp1n=kp1n-1
      if (kp1n.lt.kp1n0) kp1n=kp1n+2
      do 420 kp1=kp1n,kp1x,2
      k=kp1-1
      fk=k
      half=0.5
      zero=0.0
      flambda=s3j(fja,fk,fjb,half,zero)**2
      flmag=flmag+flambda*(eras**2)*fmi(k+1,3)/(fk*(fk+1.0))
      if (itpow.eq.4)
     1  write (9,3100) ma,mb,fkapa,fkapb,k,flambda,flmag
  420 continue
c
c
c
  500 eras=2.0/(fnorm(ma)*fnorm(mb))
      flmag=eras*flmag
      flret=-eras*flret
      if (ma.eq.mb.and.fkapa.eq.fkapb) go to 520
      eras=qa*qb
      emag=emag+eras*flmag
      eret=eret+eras*flret
      go to 530
  520 emag=emag+qa*(qa-1.0)*flmag*(fja+0.5)/(2.0*fja)
  530 continue
      if (itpow.eq.4)
     1 write (9,5300) ma,mb,fkapa,fkapb,qa,qb,eras,flmag,flret,emag,eret
 5300 format (2i4,2f5.1,2f8.4,5f12.7)
      if (fkapa.lt.0.0) go to 600
      fkapa=-fkapa-1.0
      go to 210
c
  600 if (fkapb.lt.0.0) go to 680
      fkapb=-fkapb-1.0
      go to 200
  680 continue
  690 continue
c
  700 return
      end
      subroutine brtint(n,ma,mb)
c
      implicit real*8 (a-h,o-z)
      parameter (kmsh=1801,ko=20,kmshq=1801,koq=20)
      common/c1/r(kmsh),ru(kmsh),ruee(kmsh),ruexch(kmsh),rsatom(kmsh),
     1     nnlz(ko),nnn(ko),
     2  a0m(ko),jjj(ko),wwnl(ko),l(ko),nkkk(ko),ee(ko),
     3  ei(ko),ekin(ko),een(ko),uee(ko),uex(ko),ecorr(ko),ucorr(ko),
     4  epss(ko),evel(ko),edar(ko),erel(ko),a(ko,5),b(ko),nsch(10,ko),
     5  az(ko),imat(ko),ecm(10)
      common/xij/xi(kmsh),xj(kmsh),v(kmsh),rscore(kmsh),rsvale(kmsh),
     1  recorr(kmsh),rucorr(kmsh),rc(kmsh),pnlo(kmsh),x1(kmsh),x2(kmsh)
c
      common/c2/z,iz,nconf,ncores,nvales,ncspvs,ion,mesh,kut,icut,idb,
     1  rdb,corrf,c,twozzz,nblock,kkk,tolstb,tolkm2,tolend,thresh,
     2  etot,ecorrt,ca0,ca1,emx,ekint,eentot,etotk,etotp,zzz,
     3  alfmin,alfmax,twoion,twoz,ncelec,nvelec,itpow,iptvu,xif,
     4  exf,exf10,exfm1,ipteb,ns,ib,nconft2,rm3(ko),iw6,nconft,ifrac
      common/c4/vr(kmsh),qq0(kmsh),irel,niter,nn,lam,e,end,nodf,i,j,k,
     1  mm,imatch,rsint,rsint1,hxid,izhxbw,iphfwf,ihf,ihf1,nprint
      common/brt/emag,eret,fmi(10,3),fni(10,4),fnorm(ko)
      dimension uab(kmsh),vab(kmsh)
      equivalence (uab,rscore),(vab,rsvale)
c
      save
c
      kp1=k+1
  120 xi(1)=0.0
      xj(1)=0.0
      do 130 i=2,kkk
  130 xi(i)=x1(i)/(r(i)**kp1)
      if (ma.eq.mb) go to 300
c
      do 210 i=1,kkk
  210 v(i)=uab(i)*xi(i)
      call quad5(v,1,1,kkk,a1)
      fni(kp1,n)=a1
c
  300 do 310 i=1,kkk
  310 v(i)=vab(i)*xi(i)
      call quad5(v,1,1,kkk,a1)
c     write (9,3100) n,a1,(xi(i),i=1,50),(xj(i),i=1,50),(v(i),i=1,240)
c3100 format (//i5,f14.6/(10f12.5))
      fni(kp1,n+2)=a1
      if (n.eq.1) go to 400
      fmi(kp1,1)=2.0*fni(kp1,1)
      fmi(kp1,2)=fni(kp1,2)+fni(kp1,3)
      fmi(kp1,3)=2.0*fni(kp1,4)
c     write (9,3100) n,a1,(xi(i),i=1,50),(xj(i),i=1,50),(v(i),i=1,240)
  400 continue
      if (itpow.eq.4)
     1  write (9,40) k,(fmi(kp1,i),i=1,3), (fni(kp1,i),i=1,4)
   40 format (/i5,3f10.6,5x,4f10.6)
      return
      end
      function s3j (fj1, fj2, fj3, fm1, fm2)
c
      implicit real*8 (a-h,o-z)
c
      s3j=0.0
      fm3=-fm1-fm2
      if (abs(fm3).gt.fj3) return
      a=fj1+fj2-fj3
      b=fj1-fj2+fj3
      c=-fj1+fj2+fj3
      d=fj1-fm1
      e=fj2+fm2
      x=fctrl(a)*fctrl(b)*fctrl(c)
      if (x.le.0.0) return
      x=x/fctrl(fj1+fj2+fj3+1.0)
      x=x*fctrl(d)*fctrl(fj1+fm1)*fctrl(fj2-fm2)*fctrl(e)
     1  *fctrl(fj3-fm3)*fctrl(fj3+fm3)
      if (x.le.0.0) return
      x=sqrt(x)
      b=d-b
      c=e-c
      zero=0.0
      fk=max(zero,b,c)
      k1=fk
      k2=min(a,d,e)
      if (k2.lt.k1) return
      f=(-1.0)**k1
      y=0.0
      do 150 k=k1,k2
      y=y+f/(fctrl(fk)*fctrl(fk-b)*fctrl(fk-c)*fctrl(a-fk)
     1  *fctrl(d-fk)*fctrl(e-fk))
      fk=fk+1.0
  150 f=-f
      s3j=x*y
      k=fj1-fj2-fm3
      s3j=s3j*((-1.0)**k)
      return
      end

