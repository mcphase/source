// mcdiff - output of results
void print_sps(inimcdiff & ini)
{int nofcomponents=3;
spincf spins(1,1,1,ini.natmagnetic,nofcomponents);
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
  fprintf(fout,"#! a=%8.5f b=%8.5f c=%8.5f alpha=%8.5f beta=%8.5f gamma=%8.5f\n",ini.a,ini.b,ini.c,ini.alpha,ini.beta,ini.gamma);
  fprintf(fout,"#! r1a=%8.5f r2a=%8.5f r3a=%8.5f\n",ini.nr1*ini.r1s[1],ini.nr2*ini.r2s[1],ini.nr3*ini.r3s[1]);
  fprintf(fout,"#! r1b=%8.5f r2b=%8.5f r3b=%8.5f   primitive lattice vectors [a][b][c]\n",ini.nr1*ini.r1s[2],ini.nr2*ini.r2s[2],ini.nr3*ini.r3s[2]);
  fprintf(fout,"#! r1c=%8.5f r2c=%8.5f r3c=%8.5f   (strained using strain tensor)\n",ini.nr1*ini.r1s[3],ini.nr2*ini.r2s[3],ini.nr3*ini.r3s[3]);
  fprintf(fout,"#! nofatoms=%i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n",ini.natmagnetic,spins.nofcomponents);
  fprintf(fout,"#*********************************************************************\n");

 for (i=1;i<=ini.natmagnetic;++i)
 { Vector abc(1,3);
   abc=(*ini.jjjpars[i]).xyz(1)*ini.nr1*ini.r1s+(*ini.jjjpars[i]).xyz(2)*ini.nr2*ini.r2s+(*ini.jjjpars[i]).xyz(3)*ini.nr3*ini.r3s;
   spins.m(1,1,1)(nofcomponents*(i-1)+1)=(*ini.jjjpars[i]).mom(1);
   spins.m(1,1,1)(nofcomponents*(i-1)+2)=(*ini.jjjpars[i]).mom(2);
   spins.m(1,1,1)(nofcomponents*(i-1)+3)=(*ini.jjjpars[i]).mom(3);
   fprintf(fout,"#! da=%8.5f [a] db=%8.5f [b] dc=%8.5f [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g sipffilename=%s\n",
   abc(1),abc(2),abc(3), (*ini.jjjpars[i]).paranz, (*ini.jjjpars[i]).diagonalexchange, (*ini.jjjpars[i]).gJ, (*ini.jjjpars[i]).sipffilename);
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
            ini.T,Norm(ini.H),ini.H[1],ini.H[2],ini.H[3],spins.n()*spins.nofatoms,spins.nofatoms,spins.nofcomponents);
  spins.print(fout);
 fclose(fout);
}

void print_mf(inimcdiff & ini,mfcf & mfields)
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
  fprintf(fout,"#! a=%8.5f b=%8.5f c=%8.5f alpha=%8.5f beta=%8.5f gamma=%8.5f\n",ini.a,ini.b,ini.c,ini.alpha,ini.beta,ini.gamma);
  fprintf(fout,"#! r1a=%8.5f r2a=%8.5f r3a=%8.5f\n",ini.nr1*ini.r1s[1],ini.nr2*ini.r2s[1],ini.nr3*ini.r3s[1]);
  fprintf(fout,"#! r1b=%8.5f r2b=%8.5f r3b=%8.5f   primitive lattice vectors [a][b][c]\n",ini.nr1*ini.r1s[2],ini.nr2*ini.r2s[2],ini.nr3*ini.r3s[2]);
  fprintf(fout,"#! r1c=%8.5f r2c=%8.5f r3c=%8.5f\n",ini.nr1*ini.r1s[3],ini.nr2*ini.r2s[3],ini.nr3*ini.r3s[3]);
  fprintf(fout,"#! nofatoms=%i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n",ini.natmagnetic,mfields.nofcomponents);
  fprintf(fout,"#*********************************************************************\n");

 for (i=1;i<=ini.natmagnetic;++i)
 { Vector abc(1,3);
   abc=(*ini.jjjpars[i]).xyz(1)*ini.nr1*ini.r1s+(*ini.jjjpars[i]).xyz(2)*ini.nr2*ini.r2s+(*ini.jjjpars[i]).xyz(3)*ini.nr3*ini.r3s;
   fprintf(fout,"#! da=%8.5f [a] db=%8.5f [b] dc=%8.5f [c] nofneighbours=%i diagonalexchange=%i gJ=%4.6g sipffilename=%s\n",
   abc(1),abc(2),abc(3), (*ini.jjjpars[i]).paranz, (*ini.jjjpars[i]).diagonalexchange, (*ini.jjjpars[i]).gJ, (*ini.jjjpars[i]).sipffilename);
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
            ini.T,Norm(ini.H),ini.H[1],ini.H[2],ini.H[3],mfields.n()*mfields.nofatoms,mfields.nofatoms,mfields.nofcomponents);
  mfields.print(fout);
 fclose(fout);
}

//*******************************************************************************************************
//*******************************************************************************************************
//*******************************************************************************************************
//void printheader(jjjpar ** jjjpars,int code,const char * filename,const char* infile,char * unitcell,double T,
//              Vector & H,float lambda,float ovalltemp,int lorenz,Vector r1,Vector r2,Vector r3,int n,int m,
//              float a,float b,float c,Vector & P, Vector & Pxyz)
void printheader(inimcdiff & ini,int code,int m)
{// output of the header to filename .. length is governed by abs(ini.colcod[0]);
 FILE * fout;char l[MAXNOFCHARINLINE];
 int i,ortho=1;
 double alpha,beta,gamma;
   extract(ini.unitcellstr, "alpha", alpha); extract(ini.unitcellstr, "beta", beta); extract(ini.unitcellstr, "gamma", gamma);
   if(alpha!=90||beta!=90||gamma!=90){ortho=0;}
 time_t curtime;
 struct tm * loctime;
  fout = fopen_errchk (ini.outfilename, "w");
 fprintf(fout, "#{output file of program mcdiff %s input file: %s %s ",ini.outfilename,ini.infile,MCDIFFVERSION);
 curtime=time(NULL);loctime=localtime(&curtime);fputs (asctime(loctime),fout);
 fprintf(fout,"#!<--mcdiff.mcdiff.out-->\n");
if(abs(ini.colcod[0])>0){
 fprintf(fout,"#! displaylines=false\n");
 fprintf(fout,"#***********************************************************************\n");
 fprintf(fout,"#*\n");
 fprintf(fout,"#* mcdiff - program to calculate neutron and magnetic xray diffraction\n");
 fprintf(fout,"#*\n");
 fprintf(fout,"#* reference: M. Rotter and A. Boothroyd PRB 79 (2009) 140405R\n");
 fprintf(fout,"#***********************************************************************\n");

 fprintf(fout,"# lattice parameters:%s",ini.unitcellstr);
 fprintf(fout,"# prim. unit cell   / %10.7f A \\     / %10.7f A \\     / %10.7f A \\ \n", ini.r1(1), ini.r2(1), ini.r3(1));
 fprintf(fout,"#                b1=| %10.7f A |  b2=| %10.7f A |  b3=| %10.7f A |\n", ini.r1(2), ini.r2(2), ini.r3(2));
 fprintf(fout,"#                   \\ %10.7f A /     \\ %10.7f A /     \\ %10.7f A /\n", ini.r1(3), ini.r2(3), ini.r3(3));
 fprintf(fout, "#! Wavelength=%g A   number of atoms: %i\n",ini.lambda, ini.n);
 fprintf(fout, "#! T= %g K Ha= %g T Hb= %g T Hc= %g T\n",ini.T,ini.H(1),ini.H(2),ini.H(3));
 fprintf(fout, "#! Overall temperature factor B=%g A^2: Intensity is proportional to exp(-2*B*(sin(theta)/lambda)^2)\n",ini.ovalltemp);

 if(ini.lorenz == 0){snprintf(l,MAXNOFCHARINLINE,"1.0 no lorentz factor calculated");}
 if(ini.lorenz == 1){snprintf(l,MAXNOFCHARINLINE,"1 / sin^2(2theta)   neutron powder flat sample");}
 if(ini.lorenz == 2){snprintf(l,MAXNOFCHARINLINE,"1 / sin(2theta) / sin(theta)    neutron powder cyl. sample");}
 if(ini.lorenz == 3){snprintf(l,MAXNOFCHARINLINE,"1 / sin(2theta)     neutron single crystal");}
 if(ini.lorenz == 4){snprintf(l,MAXNOFCHARINLINE,"d^3  neutron TOF powder cyl sample... log scaled d-pattern");}
 if(ini.lorenz == 5){snprintf(l,MAXNOFCHARINLINE,"d^4  neutron TOF powder cyl sample... d-pattern");}
 fprintf(fout, "# Lorentz Factor: %s\n#\n",l);
 if(code<2)
 {fprintf(fout, "# Lorentz Factor not considered for resonant magnetic xray scattering - F1 and F2 transition intensities calculated\n");
  fprintf(fout, "# according to fRMXS as given in equation (2) of Longfield et al. PRB 66 054417 (2002) and maximized with respect to azimuth.\n#\n");
 }
 fprintf(fout, "#! nofatoms=%i atoms in unit cell\n",ini.n);
 if(abs(ini.colcod[0])>1){
 fprintf(fout, "# it follows a list of atomic positions db1 db2 db3, moments m scattering lengths sl,\n");
 fprintf(fout, "# Debye Waller factor (sqr(Intensity)~|sf| ~sum_i ()i exp(-2 DWFi sin^2(theta) / lambda^2)=EXP (-Wi),\n# units DWF [A^2], relation to other notations 2*DWF=B=8 pi^2 <u^2>)\n");
 fprintf(fout, "#  and  Lande factors gJ for total angular momentum J  <j0> and <j2> formfactor\n# coefficients\n");
 if (ortho==1)
 {fprintf(fout, "#  db1[b1]db2[b2]db3[b3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
 } else
 {fprintf(fout, "#  db1[b1]db2[b2]db3[b3]mi[MuB]mj[MuB]mk[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
  fprintf(fout, "#                         ...with j||b, k||(a x b) and i normal to k and j\n");
 }
 for (i = 1;i<=ini.n;++i)
 {if((double)(i)/50==(double)(i/50))
  {
   if (ortho==1)
   {fprintf(fout, "#  db1[b1]db2[b2]db3[b3]ma[MuB]mb[MuB]mc[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
   } else
   {fprintf(fout, "#  db1[b1]db2[b2]db3[b3]mi[MuB]mj[MuB]mk[MuB]sl[10^-12cm]  DWF[A^2] gJ     <j0>:A a      B      b      C      c      D      <j2>A  a      B      b      C      c      D\n");
    fprintf(fout, "#                         ...with j||b, k||(a x b) and i normal to k and j\n");
   }
  }
  fprintf(fout, "# %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f%+6.3fi %6.3f %6.3f ",myround((*ini.jjjpars[i]).xyz(1)),myround((*ini.jjjpars[i]).xyz(2)),myround((*ini.jjjpars[i]).xyz(3)),myround((*ini.jjjpars[i]).mom(1)),myround((*ini.jjjpars[i]).mom(2)),myround((*ini.jjjpars[i]).mom(3)),myround((*ini.jjjpars[i]).SLR),myround((*ini.jjjpars[i]).SLI),myround((*ini.jjjpars[i]).DWF),myround((*ini.jjjpars[i]).gJ));
  (*ini.jjjpars[i]).FFinfo(fout);
 }
 fprintf(fout, "#}\n");
 } // abs(ini.colcod[0]) >1
 fprintf(fout, "#Intensities integrated with respect to 2theta are in barn/atom and are proportional to\n" \
               "# |SF|^2 with SF=sum_(atom j=1-nofatoms in unit cell) sl_j exp(iQr_j-DWF_j Q^2 / 16 pi^2) / nofatoms} \n");
} // abs(ini.colcod[0]) >0
 fclose(fout);
}

void printreflist(inimcdiff & ini,int code,int m,
              Vector * hkl,float * totint,float * D,float ** out,
              complex <double>*mx,complex <double>*my,complex <double>*mz,complex <double>*mxmy,
              complex <double>*mxmz,complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2
              )
{FILE * fout;
 int i,chinr=0;
 double isave[]={0,0,0,0,0,0,0,0,0,0,0,0,0};
   fout = fopen_errchk (ini.outfilename, "a");

fprintf(fout, "#    REFLECTION LIST\n");
 fprintf(fout, "#    Polarization Vector P= %g a + %g b + %g c (|P|=%g)\n",ini.P(1),ini.P(2),ini.P(3),Norm(ini.Pxyz));
 hkl[0] = 0; D[0] = 10000;  totint[0] = 0; 
 if(ini.colcod[0]<0){ mx[0]=0;my[0]=0;mz[0]=0;mx2[0]=0;my2[0]=0;mz2[0]=0;mxmy[0]=0;mxmz[0]=0;mymz[0]=0;}
 double rpvalue[]={0,0,0,0,0,0,0,0,0,0,0,0,0};
 double chisquared[]={0,0,0,0,0,0,0,0,0,0,0,0,0};
 double total=0;
 int imin=1;
 if(code==0)imin=0;
 for(i = imin;i<=m;++i)
 {if(code<2)
   {
    if((double)(i-imin)/50==(double)((i-imin)/50))
    {ini.print_usrdefcolhead(fout);
     if (ini.colcod[0]<0)
     {fprintf(fout, " F1:max-Isigpi azim Ipisig azim Ipipig azim F2:max-Isigpi azim Ipisig azim Ipipig azim " 
                    " |^ma_q| |^mb_q| |^mc_q| |^ma^2_q||^mb^2_q||^mc^2_q||(^ma*^mb)_q||(^ma*^mc)_q||(^mb*^mc)_q|\n");
     }
     else   // no magnetic xrayscattering 
     {fprintf(fout, "\n");
     }
    }
       if(ini.colcod[0]<0)
    {
   // calculate alpha_i delta_i for reflection hkl[i](1..3) 
    // mx my mz should be components along xyz of magnetic moment(Q) with
    //   y||b, z||(a x b) and x perpendicular to y and z
    double alpha1,alpha2,alpha3,delta1,delta2,delta3;
    Vector u3(1,3);u3=ini.rtoijk_rez*hkl[i];u3/=Norm(u3); // u3 in xyz=ijk coordinates
    // alpha_i are the projections of u3(Scattering Vector) onto xyz coordinates
    alpha1=acos(-0.999999*u3(1));
    alpha2=acos(-0.999999*u3(2));
    alpha3=acos(-0.999999*u3(3));
    
   //  project xyz unit vectors in plane perpendicular to q:  xperp=x-(u3*x) yperp=y-(u3*y) zperp=z-(u3*z)
    Vector xperp(1,3),yperp(1,3),zperp(1,3),u1(1,3);
     xperp=0;xperp(1)=1-u3(1);xperp(2)= -u3(1);xperp(3)= -u3(1);
     yperp=0;yperp(1)= -u3(2);yperp(2)=1-u3(2);yperp(3)= -u3(2);
     zperp=0;zperp(1)= -u3(3);zperp(2)= -u3(3);zperp(3)=1-u3(3);
    //In the chosen experimental geometry
//$\Psi=0$ when $\mbf x$ points to the x-ray source (and for the special case
// of $\mbf Q|| \mbf x$ when $\mbf y$ points to the x-ray source).
 // choose u1||xperp (unless xperp==0), else u1||yperp 
     if(Norm(xperp)>0){u1=xperp/Norm(xperp);}else{u1=yperp/Norm(yperp);}
     delta1=acos(-0.999999*(u1*xperp/Norm(xperp)));
     delta2=acos(-0.999999*(u1*yperp/Norm(yperp)));
     delta3=acos(-0.999999*(u1*zperp/Norm(zperp)));
  
    //printf("%g %g %g\n",alpha1,alpha2,alpha3);
    //printf("%g %g %g\n",delta1,delta2,delta3);exit(1);
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
      st=ini.lambda*0.5/D[i];
      ct=sqrt(1-st*st);
      s2t=2*st*ct;

      f1ps=abs(z1*ct+z3*st); if(f1ps*f1ps>IpsF1){IpsF1=f1ps*f1ps;IpsF1a=azimuth*180/PI;}
      f1sp=abs(z3*st-z1*ct); if(f1sp*f1sp>IspF1){IspF1=f1sp*f1sp;IspF1a=azimuth*180/PI;}
      f1pp=-abs(z2*s2t); if(f1pp*f1pp>IppF1){IppF1=f1pp*f1pp;IppF1a=azimuth*180/PI;}

      f2ps=abs(-z1z2*st+z2z3*ct); if(f2ps*f2ps>IpsF2){IpsF2=f2ps*f2ps;IpsF2a=azimuth*180/PI;}
      f2sp=abs(z1z2*st+z2z3*ct); if(f2sp*f2sp>IspF2){IspF2=f2sp*f2sp;IspF2a=azimuth*180/PI;}
      f2pp=abs(-ct*ct*(z12*st*st/ct/ct+z32)); if(f2pp*f2pp>IppF2){IppF2=f2pp*f2pp;IppF2a=azimuth*180/PI;}
      if(code==1)
       {ini.print_usrdefcols(fout,out,i);fprintf(fout, "        %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
       f1ps*f1ps,azimuth*180/PI,
       f1sp*f1sp,azimuth*180/PI,
       f1pp*f1pp,azimuth*180/PI,
       f2ps*f2ps,azimuth*180/PI,
       f2sp*f2sp,azimuth*180/PI,
       f2pp*f2pp,azimuth*180/PI,
       abs(mx[i]),abs(my[i]),abs(mz[i]),abs(mx2[i]),abs(my2[i]),abs(mz2[i]),abs(mxmy[i]),abs(mxmz[i]),abs(mymz[i]));}
     }
        if(IspF1+IpsF1+IppF1+IspF2+IpsF2+IppF2+totint[i]>SMALLINTENSITY)
        {ini.print_usrdefcols(fout,out,i);
         fprintf(fout, "       %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f   %8.6f %3.0f %8.6f %3.0f %8.6f %3.0f  %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
        myround(SMALLINTENSITY,IspF1),myround(SMALLINTENSITY,IspF1a),myround(SMALLINTENSITY,IpsF1),myround(SMALLINTENSITY,IpsF1a),myround(SMALLINTENSITY,IppF1),myround(SMALLINTENSITY,IppF1a),
        myround(SMALLINTENSITY,IspF2),myround(SMALLINTENSITY,IspF2a),myround(SMALLINTENSITY,IpsF2),myround(SMALLINTENSITY,IpsF2a),myround(SMALLINTENSITY,IppF2),myround(SMALLINTENSITY,IppF2a),
        myround(SMALLINTENSITY,abs(mx[i])),myround(SMALLINTENSITY,abs(my[i])),myround(SMALLINTENSITY,abs(mz[i])),myround(SMALLINTENSITY,abs(mx2[i])),myround(SMALLINTENSITY,abs(my2[i])),myround(SMALLINTENSITY,abs(mz2[i])),myround(SMALLINTENSITY,abs(mxmy[i])),
        myround(SMALLINTENSITY,abs(mxmz[i])),myround(SMALLINTENSITY,abs(mymz[i])));
        }
       }
       else
       {if(totint[i]>SMALLINTENSITY){
        ini.print_usrdefcols(fout,out,i); fprintf(fout, "\n");                               }
       }
      
    if(code==1&&ini.colcod[0]<0){fprintf(fout,"#\n");}
   }
   if(code==2)//calculate rpvalue and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {ini.print_usrdefcolhead(fout);fprintf(fout," Iobs\n");} 
      ini.print_usrdefcols(fout,out,i);fprintf(fout, "%5.4E\n",real(mx[i]));
     if(real(mx[i])>=0){total+=abs(mx[i]);
                        for(int j=1;j<= ini.nofoutputcolumns;++j)
                        {rpvalue[j]+=abs(isave[j]+out[j][i]-abs(mx[i])); 
                         isave[j]=0;}                        
                      }
     else {for(int j=1;j<= ini.nofoutputcolumns;++j)isave[j]+=out[j][i];
          }
   }
   if(code==3)//calculate also rpvalue and chisquared and output neutrons only
   {if((double)(i-imin)/50==(double)((i-imin)/50))
    {ini.print_usrdefcolhead(fout); fprintf(fout," Iobs      error\n");} 
     ini.print_usrdefcols(fout,out,i);  fprintf(fout, "%5.4E %5.4E\n",real(mx[i]),abs(my[i]));
     if(real(mx[i])>=0){total+=abs(mx[i]);
       for(int j=1;j<= ini.nofoutputcolumns;++j)
      {chisquared[j]+=(isave[j]+out[j][i]-abs(mx[i]))*(isave[j]+out[j][i]-abs(mx[i]))/abs(my[i])/abs(my[i]);
      rpvalue[j]+=abs(isave[j]+out[j][i]-abs(mx[i])); isave[j]=0;
      }
      ++chinr;
                            }
     else {for(int j=1;j<= ini.nofoutputcolumns;++j)isave[j]+=out[j][i];
          }
   }

 }
if (code>=2){for(int j=1;j<= ini.nofoutputcolumns;++j){rpvalue[j]*=100.0/total;
                                                  fprintf(fout,"#!rpvaluecol%i=%6.2f\n",j,rpvalue[j]);
                                                 }
            }
if (code==3){for(int j=1;j<= ini.nofoutputcolumns;++j){
             chisquared[j]*=1.0/(double)chinr;fprintf(fout,"#!chisquaredcol%i=%6.4f\n",j,chisquared[j]);
                                                 }            
            }

fclose(fout);
return;}

