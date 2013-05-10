/************************************************************************/
void physpropclc(Vector H,double T,spincf & sps,mfcf & mf,physproperties & physprops,inipar & ini,par & inputpars)
{ int i,j,k,l,n,m1;div_t result;float mmax; char text[100];FILE * fin_coq;
 //save fe and u
 // calculate nettomoment from spinstructure
    Vector mom(1,3),d1(1,inputpars.nofcomponents);physprops.m=0;
    for (l=1;l<=inputpars.nofatoms;++l){
    // go through magnetic unit cell and sum up the contribution of every atom
    for(i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k){
      for(m1=1;m1<=inputpars.nofcomponents;++m1){d1[m1]=mf.mf(i,j,k)[inputpars.nofcomponents*(l-1)+m1];}                  
     (*inputpars.jjj[l]).mcalc(mom,T, d1,H,(*inputpars.jjj[l]).Icalc_parstorage);
     physprops.m+=mom;
    }}}}
    physprops.m/=(double)sps.n()*(double)sps.nofatoms;


// thermal expansion - magnetostricton correlation-functions
// according to each neighbour given in mcphas.j a correlation
// function is calculated - up to ini.nofspincorrs neighbours
 int nmax,i1,j1,k1,l1,i2,j2,k2,is;Vector xyz(1,3),d(1,3),d_rint(1,3);
 nmax=ini.nofspincorrs;
 for (l=1;l<=inputpars.nofatoms;++l){if(nmax>(*inputpars.jjj[l]).paranz){nmax=(*inputpars.jjj[l]).paranz;}}

 for(n=1;n<=nmax;++n){physprops.jj[n]=0;for(l=1;l<=inputpars.nofatoms;++l){
    //calculate spincorrelation function of neighbour n of sublattice l
    // 1. transform dn(n) to primitive lattice
     xyz=(*inputpars.jjj[l]).dn[n]-(*inputpars.jjj[l]).xyz;
     d=inputpars.rez*(const Vector&)xyz;
     for (i=1;i<=3;++i)d_rint(i)=rint(d(i)); //round relative position to integer numbers (to do
                                             // something sensible if not integer, i.e. if sublattice
					     // of neighbour has not been identified by par.cpp)

     // go through magnetic unit cell and sum up the contribution of every atom
      for(i=1;i<=sps.na();++i){for(j=1;j<=sps.nb();++j){for(k=1;k<=sps.nc();++k){
         // now we are at atom ijk of sublattice l: the spin is  sps.m(i,j,k)(1..3+3*(l-1))
	 // which i1 j1 k1 l1 is the neighbour in distance d ?
	l1=(*inputpars.jjj[l]).sublattice[n];  // ... yes, this is the sublattice of the neighbour

	i1=i+(int)(d_rint(1));
	j1=j+(int)(d_rint(2));
	k1=k+(int)(d_rint(3));
        while (i1<=0) i1+=sps.na();result=div(i1,sps.na());i1=result.rem;
        while (j1<=0) j1+=sps.nb();result=div(j1,sps.nb());j1=result.rem;
        while (k1<=0) k1+=sps.nc();result=div(k1,sps.nc());k1=result.rem;
        result=div(i1,sps.na());i1=result.rem; if(i1==0)i1=sps.na();
        result=div(j1,sps.nb());j1=result.rem; if(j1==0)j1=sps.nb();
        result=div(k1,sps.nc());k1=result.rem; if(k1==0)k1=sps.nc();

	// sum up correlation function
           for(i2=1;i2<=inputpars.nofcomponents;++i2)
               {physprops.jj[n](i2+inputpars.nofcomponents*inputpars.nofcomponents*(l-1))+=
	        (sps.m(i,j,k)(i2+inputpars.nofcomponents*(l-1)))*
		(sps.m(i1,j1,k1)(i2+inputpars.nofcomponents*(l1-1))); //<JaJa>,<JbJb> ...
               }
               k2=inputpars.nofcomponents;
           for(i2=1;i2<=inputpars.nofcomponents-1;++i2)
              {for(j2=i2+1;j2<=inputpars.nofcomponents;++j2)
                        {++k2;physprops.jj[n](k2+inputpars.nofcomponents*inputpars.nofcomponents*(l-1))+=
			 (sps.m(i,j,k)(i2+inputpars.nofcomponents*(l-1)))*
			 (sps.m(i1,j1,k1)(j2+inputpars.nofcomponents*(l1-1))); //<JaJb>
			++k2;physprops.jj[n](k2+inputpars.nofcomponents*inputpars.nofcomponents*(l-1))+=
			 (sps.m(i,j,k)(j2+inputpars.nofcomponents*(l-1)))*
			 (sps.m(i1,j1,k1)(i2+inputpars.nofcomponents*(l1-1))); //<JbJa>
			}
	      }
      }}}
 }physprops.jj[n]/=sps.n(); // divide by number of basis sets in magnetic unit cell
 }

// neutron intensities (structure factor)
      // check if maxnofhklis was modified by user
  if (ini.maxnofhkls!=physprops.maxnofhkls){physprops.update_maxnofhkls(ini.maxnofhkls);}
  mmax=0; int qh,qk,ql;j=0;
  ComplexVector a(1,3),aFT(1,3),qeuklid(1,3);
    complex<double> piq(0,2*3.1415926535),g,gFT;
  ComplexVector * mq;
  Vector ri(1,3);
      Vector abc(1,6); abc(1)=inputpars.a; abc(2)=inputpars.b; abc(3)=inputpars.c;
                       abc(4)=inputpars.alpha; abc(5)=inputpars.beta; abc(6)=inputpars.gamma;
  double QQ;
  mq = new ComplexVector [mf.in(mf.na(),mf.nb(),mf.nc())+2];
  for(i=0;i<=mf.in(mf.na(),mf.nb(),mf.nc())+1;++i)
  {mq[i]=ComplexVector(1,3*mf.nofatoms);}
  float in;
  //sps.FT(mq); //Fourier trafo of momentum configuration
  for(qh=0;qh<mf.na();++qh){for(qk=0;qk<mf.nb();++qk){for(ql=0;ql<mf.nc();++ql)
     {mq[mf.in(qh,qk,ql)]=0;
     for(i=0;i<mf.na();++i){for(j=0;j<mf.nb();++j){for(k=0;k<mf.nc();++k)
      { g=exp(piq*((double)qh*i/mf.na()+(double)qk*j/mf.nb()+(double)ql*k/mf.nc()));
	for (l=1;l<=mf.nofatoms;++l) {
                                    for(m1=1;m1<=inputpars.nofcomponents;++m1){d1[m1]=mf.mf(i+1,j+1,k+1)[inputpars.nofcomponents*(l-1)+m1];}                  
                                    (*inputpars.jjj[l]).mcalc(mom,T, d1,H,(*inputpars.jjj[l]).Icalc_parstorage);
                                    mq[mf.in(qh,qk,ql)](3*(l-1)+1)+=g*mom(1);
                                    mq[mf.in(qh,qk,ql)](3*(l-1)+2)+=g*mom(2);
                                    mq[mf.in(qh,qk,ql)](3*(l-1)+3)+=g*mom(3);
                                   }
      }}}
//      fprintf(stdout,"%g %g %g %g %g \n",qh,qk,ql,imag(mq[in(qh,qk,ql)](1)),real(mq[in(qh,qk,ql)](1)));
     }}}
 for (qk=1;qk<=mf.nb();++qk)
 {for (ql=1;ql<=mf.nc();++ql) mq[mf.in(mf.na(),qk,ql)]=mq[mf.in(0,qk,ql)];}
 for (qh=1;qh<=mf.na();++qh)
 {for (ql=1;ql<=mf.nc();++ql) mq[mf.in(qh,mf.nb(),ql)]=mq[mf.in(qh,0,ql)];}
 for (qh=1;qh<=mf.na();++qh)
 {for (qk=1;qk<=mf.nb();++qk) mq[mf.in(qh,qk,mf.nc())]=mq[mf.in(qh,qk,0)];}
 mq[mf.in(mf.na(),mf.nb(),mf.nc())]=mq[mf.in(0,0,0)];

     // mq[n]=mq[0];

  j=0;for(qh=0;qh<=mf.na();++qh){for(qk=0;qk<=mf.nb();++qk){for(ql=0;ql<=mf.nc();++ql)
   {
      Vector q(1,3),hkl(1,3),Q(1,3);
      // this is q- Vector
      q(1)=(double)qh/mf.na();q(2)=(double)qk/mf.nb();q(3)=(double)ql/mf.nc();

      // try different Q vectors corresponding to q !!
     int i1,j1,k1;
//     ri(1)=inputpars.a*inputpars.r(1,1);
//     ri(2)=inputpars.b*inputpars.r(2,1);
//     ri(3)=inputpars.c*inputpars.r(3,1);
//     i1max=(int)(ini.maxQ*abs(ri)/2/3.1415926535)+1;
//     ri(1)=inputpars.a*inputpars.r(1,2);
//     ri(2)=inputpars.b*inputpars.r(2,2);
//     ri(3)=inputpars.c*inputpars.r(3,2);
//     j1max=(int)(ini.maxQ*abs(ri)/2/3.1415926535)+1;
//     ri(1)=inputpars.a*inputpars.r(1,3);
//     ri(2)=inputpars.b*inputpars.r(2,3);
//     ri(3)=inputpars.c*inputpars.r(3,3);
//     k1max=(int)(ini.maxQ*abs(ri)/2/3.1415926535)+1;

// inserted 10.5.10 to make comaptible with nonortholattices
     Matrix abc_in_ijk(1,3,1,3),p(1,3,1,3),pstar(1,3,1,3);
     get_abc_in_ijk(abc_in_ijk,abc);
     p=abc_in_ijk*inputpars.r; // p is the primitive crystal unit cell in ijk coordinates
     pstar=2*PI*p.Inverse().Transpose();
     Vector nmin(1,3),nmax(1,3);
     nlimits_calc(nmin, nmax, ini.maxQ, pstar);
     // problem: we want to find all lattice vectors Rn=ni*ai which are within a
     // sphere of radius r from the origin (ai = column vectors of matrix a)
     // this routine returns the maximum and minimum values of ni i=1,2,3
     // by probing the corners of a cube



              for (i1=(int)nmin(1);i1<=nmax(1);++i1){//for(i2=-1;i2<=1&&i2-i1!=1;i2+=2){
              for (j1=(int)nmin(2);j1<=nmax(2);++j1){//for(j2=-1;j2<=1&&j2-j1!=1;j2+=2){
              for (k1=(int)nmin(3);k1<=nmax(3);++k1){//for(k2=-1;k2<=1&&k2-k1!=1;k2+=2){
       Q(1)=q(1)+i1;Q(2)=q(2)+j1;Q(3)=q(3)+k1;
        //project back to big lattice
       hkl=inputpars.rez.Transpose()*Q;

      // qeuklid is Q in ijk coordinate system !
      hkl2ijk(ri,hkl,abc);qeuklid(1)=ri(1);qeuklid(2)=ri(2);qeuklid(3)=ri(3);
      QQ=Norm(qeuklid);

     a=0;aFT=0;
     for(l=1;l<=inputpars.nofatoms;++l)
      {//multiply mq by lattice positions exp(iqr_i) and sum into a
       ri=inputpars.rez*(const Vector&)(*inputpars.jjj[l]).xyz; // ri ... atom position with respect to primitive lattice
       gFT=exp(piq*(Q*ri));
       g=gFT*(*inputpars.jjj[l]).debyewallerfactor(QQ)*(*inputpars.jjj[l]).F(QQ); // and formfactor + debey waller factor
        for(l1=1;l1<=3;++l1){
         a(l1)+=g*mq[mf.in(qh,qk,ql)](3*(l-1)+l1);
         aFT(l1)+=gFT*mq[mf.in(qh,qk,ql)](3*(l-1)+l1);
                            }
      }
      if (QQ<ini.maxQ||(i1==0&&j1==0&&k1==0&&ini.maxQ==0))
          {
	   //calculate intensity using polarization factor (a^*.a-|a.q/|q||^2)

	   if (Q==0.0) // tests showed that a*a is square of norm(a) !!! (therefore no conjugate !!!) changed 12.4.03 !!
	    {in=abs((complex<double>)(a*a));}
           else
	    {in=abs((complex<double>)(a*a)-(a.Conjugate()*qeuklid)*(a*qeuklid)/(qeuklid*qeuklid));}

           in/=sps.n()*sps.n()*sps.nofatoms*sps.nofatoms;
           //fprintf(stdout,"%i %i %i %g %g %g\n",qh,qk,ql,in,real(a(1)),imag(a(1)));
           //not implemented: lorentzfactor

           if (in>0.0001&&ini.maxnofhkls>0)
	    {if(j<ini.maxnofhkls)
        	{++j; physprops.hkli[j](1)=hkl(1);
	         physprops.hkli[j](2)=hkl(2);
                 physprops.hkli[j](3)=hkl(3);
	         physprops.hkli[j](4)=in;
		 physprops.hkli[j](5)=abs((complex<double>)aFT(1)/(double)sps.n()/(double)sps.nofatoms);
		 physprops.hkli[j](6)=abs((complex<double>)aFT(2)/(double)sps.n()/(double)sps.nofatoms);
		 physprops.hkli[j](7)=abs((complex<double>)aFT(3)/(double)sps.n()/(double)sps.nofatoms);
		 }
	    else
	        {int jmin,jj; //determine minimal intensity in list
		 jmin=1;for(jj=1;jj<=j;++jj){if(physprops.hkli[jj](4)<physprops.hkli[jmin](4)){jmin=jj;}}
		 if (in>physprops.hkli[jmin](4)) //if calc int is bigger than smallest value in list
		    { physprops.hkli[jmin](1)=hkl(1);
	              physprops.hkli[jmin](2)=hkl(2);
                      physprops.hkli[jmin](3)=hkl(3);
	              physprops.hkli[jmin](4)=in;
	              physprops.hkli[jmin](5)=abs((complex<double>)a(1)*(1.0/(double)sps.n()/sps.nofatoms));
	              physprops.hkli[jmin](6)=abs((complex<double>)a(2)*(1.0/(double)sps.n()/sps.nofatoms));
	              physprops.hkli[jmin](7)=abs((complex<double>)a(3)*(1.0/(double)sps.n()/sps.nofatoms));}
		}
            }
           }
      }}}//}}}

   }}}physprops.nofhkls=j;
//save spinarrangement
      physprops.sps=sps;

//save mfarrangement
      physprops.mf=mf;
   spincf * magmom;
    magmom=new spincf(sps.na(),sps.nb(),sps.nc(),inputpars.nofatoms,3);
                   for (l1=1;l1<=inputpars.nofatoms;++l1){
                    // go through magnetic unit cell and sum up the contribution of every atom
                  for(i1=1;i1<=sps.na();++i1){for(j1=1;j1<=sps.nb();++j1){for(k1=1;k1<=sps.nc();++k1){
                   for(m1=1;m1<=inputpars.nofcomponents;++m1){d1[m1]=mf.mf(i1,j1,k1)[inputpars.nofcomponents*(l1-1)+m1];}                  
                   (*inputpars.jjj[l1]).mcalc(mom,T,d1,H,(*inputpars.jjj[l1]).Icalc_parstorage);
                    for(m1=1;m1<=3;++m1){(*magmom).m(i1,j1,k1)(3*(l1-1)+m1)=mom(m1);}
                    }}}}
  
// display spinstructure
  float * x;x=new float[inputpars.nofatoms+1];float *y;y=new float[inputpars.nofatoms+1];float*z;z=new float[inputpars.nofatoms+1];
		 for (is=1;is<=inputpars.nofatoms;++is)
		   {x[is]=(*inputpars.jjj[is]).xyz[1];
 		    y[is]=(*inputpars.jjj[is]).xyz[2];
		    z[is]=(*inputpars.jjj[is]).xyz[3];}
    sprintf(text,"physpropclc:T=%gK, |H|=%gT, Ha=%gT, Hb=%gT, Hc=%gT  %i spins",myround(T),myround(Norm(H)),myround(physprops.H(1)),myround(physprops.H(2)),myround(physprops.H(3)),sps.n());
                    fin_coq = fopen_errchk ("./results/.spins3dab.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,4,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dac.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,5,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
                    fin_coq = fopen_errchk ("./results/.spins3dbc.eps", "w");
                     sps.eps3d(fin_coq,text,abc,inputpars.r,x,y,z,6,inputpars.gJ,(*magmom));
                    fclose (fin_coq);
   fin_coq = fopen_errchk ("./results/.spins.eps", "w");
    sps.eps(fin_coq,text);
   fclose (fin_coq);
delete[]x;delete []y; delete []z;
  //sps.display(text);
delete []mq;
delete magmom;
}
