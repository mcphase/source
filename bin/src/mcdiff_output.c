// mcdiff - output of results
void print_sps(const char * filename,int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H)
{// fout = fopen_errchk (filename, "w");
int nofcomponents=3;
for (int i=1;i<=natmagnetic;++i){if((*jjjpars[i]).gJ==0){nofcomponents=6;}}
spincf spins(1,1,1,natmagnetic,nofcomponents);
FILE * fout;int i;
time_t curtime;
 struct tm * loctime;
   fout = fopen_errchk ("./results/mcdiff.sps", "w");
  fprintf(fout, "#{output file of program %s ",MCDIFFVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.sps-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcdiff - program to calculate neutron and magnetic Xray diffraction\n");
   fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
   fprintf(fout,"#**********************************************************\n");
   // printout the lattice and atomic positions
  fprintf(fout,"#\n# Lattice Constants (A)\n");
  fprintf(fout,"#! a=%8.5f b=%8.5f c=%8.5f alpha=%8.5f beta=%8.5f gamma=%8.5f\n",a,b,c,alpha,beta,gamma);
  fprintf(fout,"#! r1a=%8.5f r2a=%8.5f r3a=%8.5f\n",nr1*r1s[1],nr2*r2s[1],nr3*r3s[1]);
  fprintf(fout,"#! r1b=%8.5f r2b=%8.5f r3b=%8.5f   primitive lattice vectors [a][b][c]\n",nr1*r1s[2],nr2*r2s[2],nr3*r3s[2]);
  fprintf(fout,"#! r1c=%8.5f r2c=%8.5f r3c=%8.5f\n",nr1*r1s[3],nr2*r2s[3],nr3*r3s[3]);
  fprintf(fout,"#! nofatoms=%i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n",natmagnetic,spins.nofcomponents);
  fprintf(fout,"#*********************************************************************\n");

 for (i=1;i<=natmagnetic;++i)
 { Vector abc(1,3);
   abc=(*jjjpars[i]).xyz(1)*nr1*r1s+(*jjjpars[i]).xyz(2)*nr2*r2s+(*jjjpars[i]).xyz(3)*nr3*r3s;
  if((*jjjpars[i]).gJ!=0){ spins.m(1,1,1)(nofcomponents*(i-1)+1)=(*jjjpars[i]).mom(1);
   spins.m(1,1,1)(nofcomponents*(i-1)+2)=(*jjjpars[i]).mom(2);
   spins.m(1,1,1)(nofcomponents*(i-1)+3)=(*jjjpars[i]).mom(3);}
  else
 { spins.m(1,1,1)(nofcomponents*(i-1)+1)=(*jjjpars[i]).mom(4);
   spins.m(1,1,1)(nofcomponents*(i-1)+2)=(*jjjpars[i]).mom(5);
   spins.m(1,1,1)(nofcomponents*(i-1)+3)=(*jjjpars[i]).mom(6);
   spins.m(1,1,1)(nofcomponents*(i-1)+4)=(*jjjpars[i]).mom(7);
   spins.m(1,1,1)(nofcomponents*(i-1)+5)=(*jjjpars[i]).mom(8);
   spins.m(1,1,1)(nofcomponents*(i-1)+6)=(*jjjpars[i]).mom(9);}
   fprintf(fout,"#! da=%8.5f [a] db=%8.5f [b] dc=%8.5f [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g cffilename=%s\n",
   abc(1),abc(2),abc(3), (*jjjpars[i]).paranz, (*jjjpars[i]).diagonalexchange, (*jjjpars[i]).gJ, (*jjjpars[i]).cffilename);
 }
   fprintf (fout, "#!show_abc_unitcell=1.0\n");
   fprintf (fout, "#!show_primitive_crystal_unitcell=1.0\n");
   fprintf (fout, "#!show_magnetic_unitcell=1.0\n");
   fprintf (fout, "#!show_atoms=1.0\n");
   fprintf (fout, "#!spins_scale_moment=1.0\n");
   fprintf (fout, "#!scale_view_1=1.0 scale_view_2=1.0 scale_view_3=1.0\n");
   fprintf (fout, "#0 0 T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms nofmoment-components\n");
   fprintf (fout, "    #<Ma(1)> <Ma(2)> .... Momentconfiguration  \n");
   fprintf (fout, "    #<Mb(1)> <Mb(2)> .... UNITS:   [muB]\n");
   fprintf (fout, "    #<Mc(1)> <Mc(2)> ....}\n");
   fprintf (fout, " 0 0 %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i \n",
            T,Norm(H),H[1],H[2],H[3],spins.n()*spins.nofatoms,spins.nofatoms,spins.nofcomponents);
  spins.print(fout);
 fclose(fout);
}

void print_mf(const char * filename,mfcf & mfields,int natmagnetic,float a,float b,float c,float alpha,float beta,float gamma,int nr1,int nr2,int nr3,Vector r1s,Vector r2s,Vector r3s,jjjpar ** jjjpars,double T,Vector H)
{
FILE * fout;int i;
time_t curtime;
 struct tm * loctime;
   fout = fopen_errchk ("./results/mcdiff.mf", "w");
  fprintf(fout, "#{output file of program %s ",MCDIFFVERSION);
   curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
   fprintf(fout,"#!<--mcphas.mcphas.mf-->\n");
   fprintf(fout,"#*********************************************************\n");
   fprintf(fout,"# mcdiff - program to calculate neutron and magnetic Xray diffraction\n");
   fprintf(fout,"# reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
   fprintf(fout,"#**********************************************************\n");
   // printout the lattice and atomic positions
  fprintf(fout,"#\n# Lattice Constants (A)\n");
  fprintf(fout,"#! a=%8.5f b=%8.5f c=%8.5f alpha=%8.5f beta=%8.5f gamma=%8.5f\n",a,b,c,alpha,beta,gamma);
  fprintf(fout,"#! r1a=%8.5f r2a=%8.5f r3a=%8.5f\n",nr1*r1s[1],nr2*r2s[1],nr3*r3s[1]);
  fprintf(fout,"#! r1b=%8.5f r2b=%8.5f r3b=%8.5f   primitive lattice vectors [a][b][c]\n",nr1*r1s[2],nr2*r2s[2],nr3*r3s[2]);
  fprintf(fout,"#! r1c=%8.5f r2c=%8.5f r3c=%8.5f\n",nr1*r1s[3],nr2*r2s[3],nr3*r3s[3]);
  fprintf(fout,"#! nofatoms=%i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n",natmagnetic,mfields.nofcomponents);
  fprintf(fout,"#*********************************************************************\n");

 for (i=1;i<=natmagnetic;++i)
 { Vector abc(1,3);
   abc=(*jjjpars[i]).xyz(1)*nr1*r1s+(*jjjpars[i]).xyz(2)*nr2*r2s+(*jjjpars[i]).xyz(3)*nr3*r3s;
   fprintf(fout,"#! da=%8.5f [a] db=%8.5f [b] dc=%8.5f [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g cffilename=%s\n",
   abc(1),abc(2),abc(3), (*jjjpars[i]).paranz, (*jjjpars[i]).diagonalexchange, (*jjjpars[i]).gJ, (*jjjpars[i]).cffilename);
 }
   fprintf (fout, "#!show_abc_unitcell=1.0\n");
   fprintf (fout, "#!show_primitive_crystal_unitcell=1.0\n");
   fprintf (fout, "#!show_magnetic_unitcell=1.0\n");
   fprintf (fout, "#!show_atoms=1.0\n");
   fprintf (fout, "#!spins_scale_moment=1.0\n");
   fprintf (fout, "#!scale_view_1=1.0 scale_view_2=1.0 scale_view_3=1.0\n");
   fprintf (fout, "#0 0 T[K] |H| H[T] Ha[T] Hb[T] Hc[T] nofspins nofatoms nofmoment-components\n");
   fprintf (fout, "    #mfa(1) mfa(2) .... selfconsistent Mean field configuration \n");
   fprintf (fout, "    #mfb(1) mfb(2) .... UNITS: mf(i)=gJ*mu_B*heff(i)[meV] \n");
   fprintf (fout, "    #mfc(1) mfc(2) ....         (i.e. divide by gJ and mu_B=0.05788meV/T to get effective field[T]}\n");
   fprintf (fout, " 0 0 %4.4g %4.4g %4.4g  %4.4g %4.4g %i %i %i \n",
            T,Norm(H),H[1],H[2],H[3],mfields.n()*mfields.nofatoms,mfields.nofatoms,mfields.nofcomponents);
  mfields.print(fout);
 fclose(fout);
}

//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
void printeln(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
              Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int * J,int m,
              Vector * hkl,float * ikern,float * intmag,float * intmagdip,float * D,float * theta,float * out10,
              float * out11,complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,
              complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
              float a,float b,float c,int * colcode,Vector & P)
{// ausgabe auf file filename
 FILE * fout;char l[MAXNOFCHARINLINE];
 int i,j,chinr=0,ortho=1;
 double isave[]={0,0,0,0,0,0,0,0,0,0,0,0,0};
 double alpha,beta,gamma;
   extract(unitcell, "alpha", alpha); extract(unitcell, "beta", beta); extract(unitcell, "gamma", gamma);
   if(alpha!=90||beta!=90||gamma!=90){ortho=0;}
 time_t curtime;
 struct tm * loctime;
  fout = fopen_errchk (filename, "w");
 fprintf(fout, "#{output file of program mcdiff %s input file: %s %s ",filename,infile,MCDIFFVERSION);
 curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
 fprintf(fout,"#!<--mcdiff.mcdiff.out-->\n");
 fprintf(fout,"#! displaylines=false\n");
 fprintf(fout,"#***********************************************************************\n");
 fprintf(fout,"#*\n");
 fprintf(fout,"#* mcdiff - program to calculate neutron and magnetic xray diffraction\n");
 fprintf(fout,"#*\n");
 fprintf(fout,"#* reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
 fprintf(fout,"#***********************************************************************\n");

 fprintf(fout,"# unit cell:%s",unitcell);
 fprintf(fout,"#                   / %6.3f A \\     / %6.3f A \\     / %6.3f A \\ \n", r1(1), r2(1), r3(1));
 fprintf(fout,"#                r1=| %6.3f A |  r2=| %6.3f A |  r3=| %6.3f A |\n", r1(2), r2(2), r3(2));
 fprintf(fout,"#                   \\ %6.3f A /     \\ %6.3f A /     \\ %6.3f A /\n", r1(3), r2(3), r3(3));
 fprintf(fout, "#! Wavelength=%g A   number of atoms: %i\n",lambda, n);
 fprintf(fout, "#! T= %g K Ha= %g T Hb= %g T Hc= %g T\n",T,H(1),H(2),H(3));
 fprintf(fout, "#! Overall temperature factor B=%g A^2: Intensity is proportional to exp(-2*B*(sin(theta)/lambda)^2)\n",ovalltemp);

 if(lorenz == 0){sprintf(l,"100 no lorentz factor calculated");}
 if(lorenz == 1){sprintf(l,"1 / sin^2(2theta)   neutron powder flat sample");}
 if(lorenz == 2){sprintf(l,"1 / sin(2theta) / sin(theta)    neutron powder cyl. sample");}
 if(lorenz == 3){sprintf(l,"1 / sin(2theta)     neutron single crystal");}
 if(lorenz == 4){sprintf(l,"d^3  neutron TOF powder cyl sample... log scaled d-pattern");}
 if(lorenz == 5){sprintf(l,"d^4  neutron TOF powder cyl sample... d-pattern");}
 fprintf(fout, "# Lorentz Factor: %s\n#\n",l);
 if(code<2)
 {fprintf(fout, "# Lorentz Factor not considered for resonant magnetic xray scattering - F1 and F2 transition intensities calculated\n");
  fprintf(fout, "# according to fRMXS as given in equation (2) of Longfield et al. PRB 66 054417 (2002) and maximized with respect to azimuth.\n#\n");
 }
 fprintf(fout, "# List of atomic positions dr1 dr2 dr3, moments m scattering lengths sl,\n");
 fprintf(fout, "# Debye Waller factor (sqr(Intensity)~|sf| ~sum_i ()i exp(-2 DWFi sin^2(theta) / lambda^2)=EXP (-Wi),\n# units DWF [A^2], relation to other notations 2*DWF=B=8 pi^2 <u^2>)\n");
 fprintf(fout, "#  and  Lande factors total angular momentum J (=0 if dipole approximation is used) <j0> and <j2> formfactor\n# coefficients\n");
 if (ortho==1)
 {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
 } else
 {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]mi[MuB]mj[MuB]mk[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  fprintf(fout, "#                         ...with j||b, k||(a x b) and i normal to k and j\n");
 }
 for (i = 1;i<=n;++i)
 {if((double)(i)/50==(double)(i/50))
  {
   if (ortho==1)
   {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
   } else
   {fprintf(fout, "#  dr1[r1]dr2[r2]dr3[r3]mi[MuB]mj[MuB]mk[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
    fprintf(fout, "#                         ...with j||b, k||(a x b) and i normal to k and j\n");
   }
  }
  fprintf(fout, "# %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f%+6.3fi %6.3f %6.3f ",(*jjjpars[i]).xyz(1),(*jjjpars[i]).xyz(2),(*jjjpars[i]).xyz(3),(*jjjpars[i]).mom(1),(*jjjpars[i]).mom(2),(*jjjpars[i]).mom(3),(*jjjpars[i]).SLR,(*jjjpars[i]).SLI,(*jjjpars[i]).DWF,(*jjjpars[i]).gJ);
  if(J[i]==0||J[i]==-3){fprintf(fout,"F(Q) beyond dip.approx.");}
  if(J[i]==-1){fprintf(fout,"formfactor F(Q)=j0-(1-2/gJ)j2 (note that in case of transition metals mcdiff sets gJ=2) FF coefficients:");}
  if(J[i]==-2){fprintf(fout,"FL(Q)=(j0+j2)/2 and FS(Q)=j0 formfactors separate for spin and orb. moments, FF coefficients:");}
  if((*jjjpars[i]).Np(1)!=0){
  fprintf(fout," - formfactor calculated directly from radial wave function parameters in %s",(*jjjpars[i]).cffilename);
  }
  else
  {
  for (j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",(*jjjpars[i]).magFFj0(j));}
  for (j = 1;j<=7;++j)  {fprintf(fout,"%6.3f ",(*jjjpars[i]).magFFj2(j));}
  }
  fprintf(fout, "\n");
 }
 fprintf(fout, "#}\n");
 fprintf(fout, "#    REFLECTION LIST\n");
 if(colcode[10]>4||colcode[11]>4){fprintf(fout, "#    Projection Vector P= %g a + %g b + %g c \n",P(1),P(2),P(3));}
 hkl[0] = 0; D[0] = 100; theta[0] = 0; ikern[0] = 0; intmag[0] = 0;intmagdip[0] = 0; out10[0] = 0; out11[0] = 0;mx[0]=0;my[0]=0;mz[0]=0;mx2[0]=0;my2[0]=0;mz2[0]=0;mxmy[0]=0;mxmz[0]=0;mymz[0]=0;
 double rpvalue[]={0,0,0,0,0,0,0,0,0,0,0,0,0};
 double chisquared[]={0,0,0,0,0,0,0,0,0,0,0,0,0};
 double total=0;
 int imin=1;
 if(code==0)imin=0;
 for(i = imin;i<=m;++i)
 {if(code<2)
   {
    if((double)(i-imin)/50==(double)((i-imin)/50))
    {if (ortho==1)
     {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t)   Itot(2t) %s %s Imag_dip(2t) F1:max-Isigpi azim Ipisig azim Ipipig azim F2:max-Isigpi azim Ipisig azim Ipipig azim  |^ma_q| |^mb_q| |^mc_q| |^ma^2_q||^mb^2_q||^mc^2_q||(^ma*^mb)_q||(^ma*^mc)_q||(^mb*^mc)_q|}\n",colheader[colcode[10]],colheader[colcode[11]]);}
     else   // magnetic xrayscattering for nonortholattices currently not implemented
     {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t)   Itot(2t) %s %s Imag_dip(2t) \n",colheader[colcode[10]],colheader[colcode[11]]);}
    }
    // calculate alpha_i delta_i for reflection hkl[i](1..3)  [currently ok only for ortholattices !!!]
    double alpha1,alpha2,alpha3,delta1,delta2,delta3,sqr1,sqr2;
    alpha1=acos(-0.999999*hkl[i](1)*D[i]/a);   // the following lines should be extended to non ortho lattices !!!
    alpha2=acos(-0.999999*hkl[i](2)*D[i]/b);   // mind: in this section still the old convention is used: a||x,b||y,c||z ... this should be changed for nonortholattices
    alpha3=acos(-0.999999*hkl[i](3)*D[i]/c);
    delta1=acos(-1.0);
    sqr1=sqrt((1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)*(1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)+(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b)*(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b)+(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c)*(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c));
    sqr2=sqrt((1.0/b-hkl[i](2)*hkl[i](2)*D[i]*D[i]/b/b/b)*(1.0/b-hkl[i](2)*hkl[i](2)*D[i]*D[i]/b/b/b)+(hkl[i](2)*hkl[i](1)*D[i]*D[i]/b/b/a)*(hkl[i](2)*hkl[i](1)*D[i]*D[i]/b/b/a)+(hkl[i](2)*hkl[i](3)*D[i]*D[i]/b/b/c)*(hkl[i](2)*hkl[i](3)*D[i]*D[i]/b/b/c));
    delta2=acos(0.99999*hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b/b/sqr1/sqr2);
    // mind that delta2 is larger than pi if l is positive
    if(hkl[i](3)>0){delta2*=-1.0;}

    sqr1=sqrt((1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)*(1.0/a-hkl[i](1)*hkl[i](1)*D[i]*D[i]/a/a/a)+(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c)*(hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c)+(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b)*(hkl[i](1)*hkl[i](2)*D[i]*D[i]/a/a/b));
    sqr2=sqrt((1.0/c-hkl[i](3)*hkl[i](3)*D[i]*D[i]/c/c/c)*(1.0/c-hkl[i](3)*hkl[i](3)*D[i]*D[i]/c/c/c)+(hkl[i](3)*hkl[i](1)*D[i]*D[i]/c/c/a)*(hkl[i](3)*hkl[i](1)*D[i]*D[i]/c/c/a)+(hkl[i](3)*hkl[i](2)*D[i]*D[i]/c/c/b)*(hkl[i](3)*hkl[i](2)*D[i]*D[i]/c/c/b));
    delta3=acos(0.99999*hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c/c/sqr1/sqr2);
    // mind that delta3 is larger than pi if k is negative
    if(hkl[i](2)<0){delta3*=-1.0;}
    //printf("%g %g %g\n",delta1,delta2,delta3);
    // maximize IspF1 IppF1 IpsF1  IspF2 IppF2 IpsF2  and remember corresponding azimuth
    double IspF1=0,IppF1=0,IpsF1=0, IspF2=0, IppF2=0, IpsF2=0 ;
    double IspF1a=0,IppF1a=0,IpsF1a=0, IspF2a=0, IppF2a=0, IpsF2a=0;
    double azimuth;
    //printf("%g %g %g %g %g\n",hkl[i](1),hkl[i](2),hkl[i](3),delta3,hkl[i](1)*hkl[i](3)*D[i]*D[i]/a/a/c/c/sqr1/sqr2);

    for(azimuth=0.0;azimuth<=2*PI;azimuth+=PI/90)
     {complex <double> z1,z2,z3,z1z2,z2z3,z12,z32;
      double f1ps,f1sp,f1pp,f2ps,f2sp,f2pp;
      double st,ct,s2t;

      Matrix ang(1,3,1,3);
      ang(1,1)=sin(alpha1)*cos(azimuth+delta1);
      ang(1,2)=sin(alpha2)*cos(azimuth+delta2);
      ang(1,3)=sin(alpha3)*cos(azimuth+delta3);
      ang(2,1)=sin(alpha1)*sin(azimuth+delta1);
      ang(2,2)=sin(alpha2)*sin(azimuth+delta2);
      ang(2,3)=sin(alpha3)*sin(azimuth+delta3);
      ang(3,1)=cos(alpha1);
      ang(3,2)=cos(alpha2);
      ang(3,3)=cos(alpha3);

      z1=mx[i]*ang(1,1)+my[i]*ang(1,2)+mz[i]*ang(1,3);
      z2=mx[i]*ang(2,1)+my[i]*ang(2,2)+mz[i]*ang(2,3);
      z3=mx[i]*ang(3,1)+my[i]*ang(3,2)+mz[i]*ang(3,3);

      z1z2=mx2[i]*ang(1,1)*ang(2,1);
      z1z2+=my2[i]*ang(1,2)*ang(2,2);
      z1z2+=mz2[i]*ang(1,3)*ang(2,3);
      z1z2+=mxmy[i]*(ang(1,1)*ang(2,2)+ang(1,2)*ang(2,1));
      z1z2+=mxmz[i]*(ang(1,1)*ang(2,3)+ang(1,3)*ang(2,1));
      z1z2+=mymz[i]*(ang(1,2)*ang(2,3)+ang(1,3)*ang(2,2));

      z2z3=mx2[i]*ang(2,1)*ang(3,1);
      z2z3+=my2[i]*ang(2,2)*ang(3,2);
      z2z3+=mz2[i]*ang(2,3)*ang(3,3);
      z2z3+=mxmy[i]*(ang(2,1)*ang(3,2)+ang(2,2)*ang(3,1));
      z2z3+=mxmz[i]*(ang(2,1)*ang(3,3)+ang(2,3)*ang(3,1));
      z2z3+=mymz[i]*(ang(2,2)*ang(3,3)+ang(2,3)*ang(3,2));

      z12=mx2[i]*ang(1,1)*ang(1,1);
      z12+=my2[i]*ang(1,2)*ang(1,2);
      z12+=mz2[i]*ang(1,3)*ang(1,3);
      z12+=mxmy[i]*(ang(1,1)*ang(1,2)+ang(1,2)*ang(1,1));
      z12+=mxmz[i]*(ang(1,1)*ang(1,3)+ang(1,3)*ang(1,1));
      z12+=mymz[i]*(ang(1,2)*ang(1,3)+ang(1,3)*ang(1,2));

      z32=mx2[i]*ang(3,1)*ang(3,1);
      z32+=my2[i]*ang(3,2)*ang(3,2);
      z32+=mz2[i]*ang(3,3)*ang(3,3);
      z32+=mxmy[i]*(ang(3,1)*ang(3,2)+ang(3,2)*ang(3,1));
      z32+=mxmz[i]*(ang(3,1)*ang(3,3)+ang(3,3)*ang(3,1));
      z32+=mymz[i]*(ang(3,2)*ang(3,3)+ang(3,3)*ang(3,2));
      st=sin(theta[i]*PI/180);
      ct=cos(theta[i]*PI/180);
      s2t=sin(2.0*theta[i]*PI/180);

      f1ps=abs(z1*ct+z3*st); if(f1ps*f1ps>IpsF1){IpsF1=f1ps*f1ps;IpsF1a=azimuth*180/PI;}
      f1sp=abs(z3*st-z1*ct); if(f1sp*f1sp>IspF1){IspF1=f1sp*f1sp;IspF1a=azimuth*180/PI;}
      f1pp=-abs(z2*s2t); if(f1pp*f1pp>IppF1){IppF1=f1pp*f1pp;IppF1a=azimuth*180/PI;}

      f2ps=abs(-z1z2*st+z2z3*ct); if(f2ps*f2ps>IpsF2){IpsF2=f2ps*f2ps;IpsF2a=azimuth*180/PI;}
      f2sp=abs(z1z2*st+z2z3*ct); if(f2sp*f2sp>IspF2){IspF2=f2sp*f2sp;IspF2a=azimuth*180/PI;}
      f2pp=abs(-ct*ct*(z12*st*st/ct/ct+z32)); if(f2pp*f2pp>IppF2){IppF2=f2pp*f2pp;IppF2a=azimuth*180/PI;}
      if(code==1&&ortho==1)
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %5.4E %8.4f %5.4E %5.4E %5.4E        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
       hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],myround(ikern[i]), myround(intmag[i]),myround( (ikern[i]+intmag[i])),myround(out10[i]),
       myround(out11[i]),myround(intmagdip[i]),
       f1ps*f1ps,azimuth*180/PI,
       f1sp*f1sp,azimuth*180/PI,
       f1pp*f1pp,azimuth*180/PI,
       f2ps*f2ps,azimuth*180/PI,
       f2sp*f2sp,azimuth*180/PI,
       f2pp*f2pp,azimuth*180/PI,
       abs(mx[i]),abs(my[i]),abs(mz[i]),abs(mx2[i]),abs(my2[i]),abs(mz2[i]),abs(mxmy[i]),abs(mxmz[i]),abs(mymz[i]));}
     }

    if(IspF1+IpsF1+IppF1+IspF2+IpsF2+IppF2+ikern[i]+intmag[i]>0.0001)
      {if(ortho==1)
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %5.4E %8.4f %5.4E %5.4E %5.4E        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
        hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],
        myround(ikern[i]), myround(intmag[i]), myround((ikern[i]+intmag[i])),myround(out10[i]),myround(out11[i]),
        myround(intmagdip[i]),
        IspF1,IspF1a,IpsF1,IpsF1a,IppF1,IppF1a,IspF2,IspF2a,IpsF2,IpsF2a,IppF2,IppF2a,
        abs(mx[i]),abs(my[i]),abs(mz[i]),abs(mx2[i]),abs(my2[i]),abs(mz2[i]),abs(mxmy[i]),abs(mxmz[i]),abs(mymz[i]));}
       else
       {fprintf(fout, "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %5.4E %8.4f %5.4E %5.4E %5.4E\n",
        hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],
        myround(ikern[i]), myround(intmag[i]), myround((ikern[i]+intmag[i])),
        myround(out10[i]),myround(out11[i]),myround(intmagdip[i]));}
      }
    if(code==1&&ortho==1){fprintf(fout,"#\n");}
   }
   if(code==2)//calculate rpvalue and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t)   Itot(2t) %s %s Imag_dip(2t) Iobs\n",colheader[colcode[10]],colheader[colcode[11]]);}
      fprintf(fout,   "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %5.4E %8.4f %5.4E %5.4E %5.4E %8.4f\n",
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],
      myround(ikern[i]), myround(intmag[i]),myround((ikern[i]+intmag[i])),
      myround(out10[i]),myround(out11[i]),myround(intmagdip[i]),real(mx[i]));
     if(real(mx[i])>=0){total+=abs(mx[i]);
                        rpvalue[9]+=abs(isave[9]+ikern[i]+intmag[i]-abs(mx[i])); isave[9]=0;
                        rpvalue[10]+=abs(isave[10]+out10[i]-abs(mx[i])); isave[10]=0;
                        rpvalue[11]+=abs(isave[11]+out11[i]-abs(mx[i])); isave[11]=0;
                        rpvalue[12]+=abs(isave[12]+ikern[i]+intmagdip[i]-abs(mx[i]));isave[12]=0;
                      }
     else {isave[9]+=ikern[i]+intmag[i];isave[10]+=out11[i];isave[11]+=out11[i];isave[12]+=ikern[i]+intmagdip[i];
          }
   }
   if(code==3)//calculate also rpvalue and chisquared and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {fprintf(fout, "#{h     k      l      d[A]    |Q|[1/A] 2theta  Inuc(2t) Imag(2t)   Itot(2t) %s %s Imag_dip(2t) Iobs        error\n",colheader[colcode[10]],colheader[colcode[11]]);}
      fprintf(fout,    "%6.3f %6.3f %6.3f %7.4f %7.4f %7.3f %8.4f %5.4E %8.4f %5.4E %5.4E %5.4E %5.4E %5.4E\n",
      hkl[i](1), hkl[i](2), hkl[i](3),D[i],2 * PI / D[i],2 * theta[i],
      myround(ikern[i]), myround(intmag[i]), myround((ikern[i]+intmag[i])),
      myround(out10[i]),myround(out11[i]),myround(intmagdip[i]),
      real(mx[i]),abs(my[i]));
     if(real(mx[i])>=0){total+=abs(mx[i]);
      chisquared[9]+=(isave[9]+ikern[i]+intmag[i]-abs(mx[i]))*(isave[9]+ikern[i]+intmag[i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      chisquared[10]+=(isave[10]+out10[i]-abs(mx[i]))*(isave[10]+out10[i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      chisquared[11]+=(isave[11]+out11[i]-abs(mx[i]))*(isave[11]+out11[i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      chisquared[12]+=(isave[12]+ikern[i]+intmagdip[i]-abs(mx[i]))*(isave[12]+ikern[i]+intmagdip[i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      ++chinr;
      rpvalue[9]+=abs(isave[9]+ikern[i]+intmag[i]-abs(mx[i])); isave[9]=0;
      rpvalue[10]+=abs(isave[10]+out10[i]-abs(mx[i])); isave[10]=0;
      rpvalue[11]+=abs(isave[11]+out11[i]-abs(mx[i])); isave[11]=0;
      rpvalue[12]+=abs(isave[12]+ikern[i]+intmagdip[i]-abs(mx[i])); isave[12]=0;
                      }
     else {isave[9]+=ikern[i]+intmag[i];isave[10]+=out11[i];isave[11]+=out11[i];isave[12]+=ikern[i]+intmagdip[i];
          }
   }

 }
if (code>=2){rpvalue[9]*=100.0/total;fprintf(fout,"#!rpvalue=%6.2f\n",rpvalue[9]);
             rpvalue[10]*=100.0/total;fprintf(fout,"#!rpvaluecol10=%6.2f\n",rpvalue[10]);
             rpvalue[11]*=100.0/total;fprintf(fout,"#!rpvaluecol11=%6.2f\n",rpvalue[11]);
             rpvalue[12]*=100.0/total;fprintf(fout,"#!rpvaluedip=%6.2f\n",rpvalue[12]);}
if (code==3){chisquared[9]*=1.0/(double)chinr;fprintf(fout,"#!chisquared=%6.4f\n",chisquared[9]);
             chisquared[10]*=1.0/(double)chinr;fprintf(fout,"#!chisquaredcol10=%6.4f\n",chisquared[10]);
             chisquared[11]*=1.0/(double)chinr;fprintf(fout,"#!chisquaredcol11=%6.4f\n",chisquared[11]);
             chisquared[12]*=1.0/(double)chinr;fprintf(fout,"#!chisquareddip=%6.4f\n",chisquared[12]);}

fclose(fout);
return;}

