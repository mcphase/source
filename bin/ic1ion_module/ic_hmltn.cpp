/* ic_hmltn.cpp
 *
 * Calculates the intermediate coupling Hamiltonian matrix.
 *
 * Functions:
 *    sMat<double> convH2H(sMat<double> H,int l,std::vector<std::vector<int> > c); // Converts from one basis to another
 *    bool ic_parseheader(const char *filename, icpars &pars);                     // Determines the paramters from saved file
 *    sMat<double> ic_Hcso(icpars &pars);                                          // Calculates the Hamiltonian H=H_c+H_so
 *    sMat<double> ic_hmltn(sMat<double> &H_cfi, icpars &pars);                    // Calculates the IC Hamiltonian matrix
 *    std::vector<double> ic_mag(sMat<double> &Hic, sMat<double> &iHic,            // Calculates the magnetisation
 *           sMat<double> &Jmat, sMat<double> &iJmat, std::vector<double> &T, 
 *           double H_mag=1., int nev=30);
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include <fstream>
#include <sstream>
#include <string>

// --------------------------------------------------------------------------------------------------------------- //
// Routine to convert a matrix from one basis to another, e.g from |vUSL> to |vUSLJmJ>
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> convH2H(sMat<double> Hin, int lnOut, std::vector< std::vector<int> > cv)
{
   int id,i,j,m,n,ldiag;
   sMat<double> Hout(lnOut,lnOut);
   std::vector< std::vector<int> > nz = Hin.find();
   int sz = (int)nz.size();

   for(id=0; id<sz; id++)
   {
      m = cv[nz[id][0]][1]-cv[nz[id][0]][0]+1; n = cv[nz[id][1]][1]-cv[nz[id][1]][0]+1;
      ldiag = (m<n) ? m : n;
      j = cv[nz[id][1]][0];
      for(i=cv[nz[id][0]][0]; i<(cv[nz[id][0]][0]+ldiag); i++)
         Hout(i,j++) = Hin(nz[id][0],nz[id][1]);
   }

   return Hout;
}
void convH2H(complexdouble *Hin, complexdouble *Hout, int lnIn, int lnOut, std::vector< std::vector<int> > cv)
{
   int i,j,ii,jj,m,n,ldiag; memset(Hout,0,lnOut*lnOut*sizeof(complexdouble));

   for(i=0; i<lnIn; i++)
      for(j=0; j<lnIn; j++)
      {
         if(fabs(Hin[i+j*lnIn].r)>DBL_EPSILON || fabs(Hin[i+j*lnIn].i)>DBL_EPSILON)
         {
            m = cv[i][1]-cv[i][0]+1; n = cv[j][1]-cv[j][0]+1; ldiag = (m<n) ? m : n; jj = cv[j][0];
            for(ii=cv[i][0]; ii<(cv[i][0]+ldiag); ii++)
            {
               Hout[ii+jj*lnOut].r = Hin[i+j*lnIn].r; Hout[ii+jj*lnOut].i = Hin[i+j*lnIn].i; jj++;
            }
         }
      }
}

// --------------------------------------------------------------------------------------------------------------- //
// Reads the header of a file and determines if the parameters saved match that of pars
// --------------------------------------------------------------------------------------------------------------- //
bool ic_parseheader(const char *filename, icpars &pars)
{
   std::string ss,sp,e_units,b_units,norm,type;
   std::stringstream st;
   size_t i,j,ii;
   int n,k,q; orbital l;
   std::vector<double> F(4,0.), alpha(3,0.); double xi=0.;
   cfpars B;

   std::fstream FILEIN; FILEIN.open(filename, std::fstream::in);
   if(FILEIN.fail()==true) return false;

   getline(FILEIN,ss);             // Gets the first line % Date and module info
   getline(FILEIN,ss);             // Free ion configuration
   if(ss.find("Configuration")==std::string::npos) getline(FILEIN,ss); // Ion Name :  (if exists)
   i = ss.find_last_of("^");
   n = ss[i+1]-48; sp = ss.substr(i-1,1); l = Lin(sp);
   if(n!=pars.n || l!=pars.l) { FILEIN.close(); return false; }
   getline(FILEIN,ss);             // Free ion parameters
   i = ss.find("("); j = ss.find(")",i); e_units = ss.substr(i+1,j-i-1);
   i = ss.find("=",j); j = ss.find("F",i); st.str(ss.substr(i+1,j-i-1)); st >> F[1];
   i = ss.find("=",j); if(l==2) j = ss.find("x",i); else j = ss.find("F",i); st.str(ss.substr(i+1,j-i-1)); st >> F[2];
   if(l==3) { i = ss.find("=",j); j = ss.find("x",i); sp = ss.substr(i+1,j-i-1); st >> F[3]; }
   i = ss.find("=",j); j = ss.find("a",i); st.str(ss.substr(i+1,j-i-1)); st >> xi;
   i = ss.find("=",j); j = ss.find("b",i); st.str(ss.substr(i+1,j-i-1)); st >> alpha[0];
   if(l==3) { i = ss.find("=",j); j = ss.find("g",i); st.str(ss.substr(i+1,j-i-1)); st >> alpha[1];
              i = ss.find("=",j); st.str(ss.substr(i+1)); st >> alpha[2]; }
   else { i = ss.find("=",j); st.str(ss.substr(i+1)); st >> alpha[1]; }
   conv_e_units(pars,e_units);
   if(fabs(xi-pars.xi)>1e-5) { FILEIN.close(); return false; }
   if(fabs(F[0]-pars.F[0])>1e-5 || fabs(F[1]-pars.F[1])>1e-5 || fabs(F[2]-pars.F[2])>1e-5 || fabs(F[3]-pars.F[3])>1e-5)
      { FILEIN.close(); return false; }
   if(fabs(alpha[0]-pars.alpha[0])>1e-5 || fabs(alpha[1]-pars.alpha[1])>1e-5 || fabs(alpha[2]-pars.alpha[2])>1e-5) 
      { FILEIN.close(); return false; }
   getline(FILEIN,ss);             // Crystal Field Normalisation
   i = ss.find(":"); norm = ss.substr(i+2);
   getline(FILEIN,ss);             // Crystal Field parameters
   FILEIN.close();
   i = ss.find("("); j = ss.find(")",i); b_units = ss.substr(i,j-i);
   B.conv_e_units(b_units); B.conv_B_norm(norm);
   i = ss.find_first_of("ABDLVW",j); ii = ss.find_first_of("246");
   type = ss.substr(i,ii-i); B.conv(type); i = ss.find(type,j); 
   while(i!=std::string::npos)
   {
      i+=type.length(); ii = ss.find("=",i); j = ss.find(",",ii);
      k = ss[i]-48; i++; if(ss[i]==61||ss[i]==45) q = -(ss[i+1]-48); else q = ss[i]-48; // 61==m 45==- 48==0
      sp = ss.substr(ii+1,j-ii-1); //val = atof(sp.c_str());
      B.assign(type,k,q,atof(sp.c_str())); i = ss.find(type,j);
   }

   return B==pars.B;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculate the Hamiltonian matrix: H = H_c + H_so for the Coulomb and spin-orbit terms only
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> ic_Hcso(icpars &pars)
{
   sMat<double> emat,H_so;
   std::vector<double> E;
   int n = pars.n; orbital e_l = pars.l;
   int nn = n; if(n>(2*e_l+1)) n = 4*e_l+2-n; 

   switch(e_l)
   {
      case F: 
         E = racah_FtoF_k(pars._F); E = racah_FtoE(E);
         emat = racah_emat(nn,E[0],E[1],E[2],E[3]) + racah_ci(nn,pars._alpha[0],pars._alpha[1],pars._alpha[2]);
         H_so = racah_so(n,pars._xi);
         break;
      case D:
         emat = racah_emat(nn,pars._F[0],pars._F[1],pars._F[2]) + racah_ci(nn,pars._alpha[0],pars._alpha[1]);
         H_so = racah_so(n,pars._xi,D);
         break;
      case P: 
         emat = racah_emat(nn,pars._F[0],pars._F[1]) + racah_ci(nn,pars._alpha[0]); H_so = racah_so(n,pars._xi,D);
         break;
      default:
         std::cerr << "ic_hmltn(): l!=2 or l!=3, only d- and f- electrons have been implemented so far.\n";
   }
 
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
   int i,icv,J2min,J2max;
   std::vector< std::vector<int> > cvSI2SO;
   std::vector<int> cvEl;

   cvSI2SO.reserve(num_states); icv=0;
   // Goes through all the states and gets the conversion matrix to use with convH2H, which is just a running index of blocks
   //    of the same L,S values (different J's) for |vUSL> to |vUSLJ> (cvSI2SO), and same L,S,J for |vUSL> to |vUSLJmJ> (cvSI2CF)
   for(i=0; i<num_states; i++)
   {
      J2min = abs(abs(conf.states[i].L)*2 - conf.states[i].S2);
      J2max =     abs(conf.states[i].L)*2 + conf.states[i].S2;
      cvEl.push_back(icv); icv += (J2max-J2min)/2+1; cvEl.push_back(icv-1); cvSI2SO.push_back(cvEl); cvEl.clear();
   }

   if(nn>(2*e_l+1)) H_so = convH2H(emat,icv,cvSI2SO) - H_so; else H_so += convH2H(emat,icv,cvSI2SO);

   return H_so;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculate the intermediate coupling Hamiltonian matrix
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> ic_hmltn(sMat<double> &H_cfi, icpars &pars)
{
 //char rmat[255],imat[255]; strcpy(rmat,"results/mcphas.icmatr"); strcpy(imat,"results/mcphas.icmati");
   sMat<double> H_cf;
 //if(ic_parseheader(rmat,pars) && ic_parseheader(imat,pars)) { H_cf = mm_gin(rmat); H_cfi = mm_gin(imat); return H_cf; }

   sMat<double> Upq,Umq,emat,H_so;
   std::vector<double> E;
   double p = 1./( pow(-1.,(double)abs(pars.l))*(2.*pars.l+1.) );
   double icfact[] = {0,0,p/threej(2*pars.l,4,2*pars.l,0,0,0),0,p/threej(2*pars.l,8,2*pars.l,0,0,0),0,p/threej(2*pars.l,12,2*pars.l,0,0,0)};
   int n = pars.n; orbital e_l = pars.l;
   int nn = n; //if(n>(2*e_l+1)) n = 4*e_l+2-n; 

   switch(e_l)
   {
      case F: 
         E = racah_FtoF_k(pars._F); E = racah_FtoE(E);
         emat = racah_emat(nn,E[0],E[1],E[2],E[3]) + racah_ci(nn,pars._alpha[0],pars._alpha[1],pars._alpha[2]);
         H_so = racah_so(n,pars._xi);
         break;
      case D:
         emat = racah_emat(nn,pars._F[0],pars._F[1],pars._F[2]) + racah_ci(nn,pars._alpha[0],pars._alpha[1]);
         H_so = racah_so(n,pars._xi,D);
         break;
      case P: 
         emat = racah_emat(nn,pars._F[0],pars._F[1]) + racah_ci(nn,pars._alpha[0]); H_so = racah_so(n,pars._xi,P);
         break;
      default:
         std::cerr << "ic_hmltn(): l!=2 or l!=3, only d- and f- electrons have been implemented so far.\n";
   }
 
   int k,iq,q;
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
   int i,j,icv,icv1,icv2,sum2Jmin_JmaxP1;
   int J2min,J2max;
   std::vector< std::vector<int> > cvSI2SO,cvSI2CF,cvSO2CF;
   std::vector<int> cvEl;

   char nstr[6]; char basename[255]; char filename[255]; strcpy(basename,"results/mms/");
   if(pars.save_matrices) {
   #ifndef _WINDOWS
   struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
      if(mkdir("results/mms",0777)!=0) std::cerr << "ic_hmltn(): Can't create mms dir, " << strerror(errno) << "\n";
   #else
   DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
      if (!CreateDirectory("results\\mms", NULL)) std::cerr << "ic_hmltn(): Cannot create mms directory\n";
   #endif
   nstr[0] = (e_l==F?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
   } else { strcpy(basename,"nodir/"); }

   cvSI2SO.reserve(num_states); cvSI2CF.reserve(num_states); cvSO2CF.reserve(num_states*5); icv=0; icv1=0; icv2=0;
   // Goes through all the states and gets the conversion matrix to use with convH2H, which is just a running index of blocks
   //    of the same L,S values (different J's) for |vUSL> to |vUSLJ> (cvSI2SO), and same L,S,J for |vUSL> to |vUSLJmJ> (cvSI2CF)
   for(i=0; i<num_states; i++)
   {
      J2min = abs(abs(conf.states[i].L)*2 - conf.states[i].S2);
      J2max =     abs(conf.states[i].L)*2 + conf.states[i].S2;
      sum2Jmin_JmaxP1 = 0; for(j=J2min; j<=J2max; j+=2) sum2Jmin_JmaxP1 += j+1;
      cvEl.push_back(icv1); icv1 += (J2max-J2min)/2+1; cvEl.push_back(icv1-1); cvSI2SO.push_back(cvEl); cvEl.clear();
      cvEl.push_back(icv2); icv2 += sum2Jmin_JmaxP1;   cvEl.push_back(icv2-1); cvSI2CF.push_back(cvEl); cvEl.clear();
      for(j=J2min; j<=J2max; j+=2)
      {
         cvEl.push_back(icv); icv += j+1; cvEl.push_back(icv-1); cvSO2CF.push_back(cvEl); cvEl.clear();
      }
   }

   sMat<double> temat = convH2H(emat,icv1,cvSI2SO);
   if(nn>(2*e_l+1)) H_so = temat - H_so; else H_so += temat; temat.clear();

   H_cf.zero(icv,icv); H_cfi.zero(icv,icv);

// Calculates the crystal field Hamiltonian
//
//         ---   k    [  k        q  k    ]     ---  k  k     ---   k [  k       q  k  ]
// V   = i >    B     | O   - (-1)  O     |  +  >   B  O   +  >    B  | O  + (-1)  O   |
//  cf     ---   -|q| [  |q|         -|q| ]     ---  0  0     ---   q [  q          -q ]
//        k,q<0                                  k           k,q>0  
//
#define NSTR(K,Q) nstr[1] = K+48; nstr[2] = Q+48; nstr[3] = 0
#define MSTR(K,Q) nstr[1] = K+48; nstr[2] = 109;  nstr[3] = Q+48; nstr[4] = 0
   for(k=2; k<=6; k+=2)
      for(iq=0; iq<(2*k+1); iq++)           //  racah_ukq is slightly faster as it involves less operations. However,
      {                                     //     it needs more memory as it stores zero-value entries as well.
         q = iq-k;                          //  For n=5, it needs around 1G, or it hits cache (slow!), for n=6,7 it
         if(pars.B(k,q)!=0)                 //     needs around 3G. 
         {                                  //  fast_ukq does slightly more operations, uses less memory as it only
            if(q==0)                        //     calculates for nonzero reduced matrix elements. It needs around
            {                               //     300Mb for n=5 and about 1.8G for n=6,7
               NSTR(k,0); 
               strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Upq = mm_gin(filename);
               if(Upq.isempty())
               {
                //if(n<6)
                     Upq = racah_ukq(n,k,0,e_l);
                //else
                //   Upq = fast_ukq(n,k,0);
                  rmzeros(Upq); mm_gout(Upq,filename);
               }
             /*if(nn>(2*e_l+1)) H_cf -= Upq * (pars.B(k,q)/icfact[k]); else*/ H_cf += Upq * (pars.B(k,q)/icfact[k]);
            }
            else
            {
               NSTR(k,abs(q)); 
               strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Upq = mm_gin(filename);
               if(Upq.isempty())
               {
                /*if(n<6)*/Upq = racah_ukq(n,k,abs(q),e_l);
                //else Upq = fast_ukq(n,k,abs(q));
                  rmzeros(Upq); mm_gout(Upq,filename);
               }
               MSTR(k,abs(q)); 
               strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
               Umq = mm_gin(filename);
               if(Umq.isempty())
               {
                /*if(n<6)*/Umq = racah_ukq(n,k,-abs(q),e_l);
                //else Umq = fast_ukq(n,k,-abs(q));
                  rmzeros(Umq); mm_gout(Umq,filename);
               }
               if(q<0)
                /*if(nn>(2*e_l+1)) H_cfi-= (Upq - Umq*pow(-1.,q)) * (pars.B(k,q)/icfact[k]); 
                  else   H_cfi+= (Upq - Umq*pow(-1.,q)) * (pars.B(k,q)/icfact[k]); changed MR 15.12.09 */
                         H_cfi-= (Upq - Umq*pow(-1.,q)) * (pars.B(k,q)/icfact[k]);
               else
                /*if(nn>(2*e_l+1)) H_cf -= (Upq + Umq*pow(-1.,q)) * (pars.B(k,q)/icfact[k]); 
                  else*/ H_cf += (Upq + Umq*pow(-1.,q)) * (pars.B(k,q)/icfact[k]); 
            }
         }
      }

   H_cf += convH2H(H_so,icv2,cvSO2CF); rmzeros(H_cf); rmzeros(H_cfi);

 //ic_printheader(rmat,pars); mm_gout(H_cf,rmat);
 //ic_printheader(imat,pars); mm_gout(H_cfi,imat);

   return H_cf;
}
// --------------------------------------------------------------------------------------------------------------- //
// Function to calculate the magnetisation given an IC Hamiltonian matrix and applied magnetic field
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> ic_mag(sMat<double> &Hic, sMat<double> &iHic, sMat<double> &Jmat, sMat<double> &iJmat, std::vector<double> &T, double H_mag, int nev) 
{                                                                                       //H_mag=1., nev=30 by default
   std::vector<double> mag(T.size(),0.);

   if(Hic.nr() != Hic.nc()) { std::cerr << "ic_mag: Hic not square!\n"; return mag; }

   int i,j,Esz,ind_j;//,ind_i,id_j;
   double dt,dt2,me,*vt=0;
   std::vector<int> ind,i_ndegen;
   std::vector<double> E,me_ndegen,Z(T.size(),0.);
   sMat<double> V,Mt,iV,iMt;
   complexdouble *zt=0,zme;

   double *fJmat=0; complexdouble *zJmat=0;
   if(iJmat.isempty() && iHic.isempty()) fJmat = Jmat.f_array(); else zJmat = zmat2f(Jmat,iJmat);

   // Physical constants. Taken from NIST Reference on Constants, Units, and Uncertainty,
   //     http://physics.nist.gov/cuu/Constants/
   const double mu_B  = 0.46686437;       // cm^{-1} / Tesla - Bohr magneton
   const double k_B   = 1.3806505e-23;    // J/K - Boltzmann constant
   const double h     = 6.62606896e-34;   // Js - Planck's constant
   const double c     = 299792458;        // m/s - speed of light in vacuum

   int Hsz = Hic.nr();

   // Defines smallness criteria
 //double small = 1e-6 * (Q_e/1000);                                  // small = 1e-6 * (Q_e/1000);

   // Calculates the total Hamiltonian as a function of field (last index)
   sMat<double> Hmltn = Hic + Jmat*(-H_mag*mu_B);                     // Hmltn = Hic + (-H(ind_H)*mu_B).*Jmat;
   sMat<double> iHmltn = iHic + iJmat*(-H_mag*mu_B);

   // Calculates the eigenvectors V and eigenvalues (enegies) E
   // Where:            ---
   //        | V  >  =  >    a  |j, j    >
   //           i       ---i  i      z,i
   //
   int lda = Hmltn.nc(), n = Hmltn.nr(), info = 0;
   char jobz = 'V', uplo = 'U';
   int lwork = 4*n, incx = 1;
   double *eigval = (double*)malloc(n*sizeof(double));
   double *m=0, *work=0, *mz=0;
   complexdouble *z=0, *zwork=0, *zm=0;
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   double alpha = 1, beta = 0;
   if(iHmltn.isempty())
   {
      char range = 'A'; double vl,vu; int il=1,iu=n/10,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = (int*)malloc(2*n*sizeof(int));
      m = (double*)malloc((ldz*n)*sizeof(double));
      lwork = 26*n; int liwork=10*n;
      mz = Hmltn.f_array(); work = (double*)malloc(lwork*sizeof(double));
      int *iwork = (int*)malloc(liwork*sizeof(int));
      F77NAME(dsyevr)(&jobz, &range, &uplo, &n, mz, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, m, &ldz, isuppz, work, &lwork, iwork, &liwork, &info);
      free(isuppz); free(work); free(iwork); free(mz);
   }
   else
   {
      char range = 'A'; double vl,vu; int il,iu,numfnd,ldz=n;
      double abstol = 0.00001; int *isuppz = (int*)malloc(2*n*sizeof(int));
      zm = zmat2f(Hmltn,iHmltn); zwork = (complexdouble*)malloc(lwork*sizeof(complexdouble));
      z = (complexdouble*)malloc((ldz*n)*sizeof(complexdouble));
      int lrwork = 24*n,liwork=10*n;
      double *rwork = (double*)malloc(lrwork*sizeof(double));
      int *iwork = (int*)malloc(liwork*sizeof(int));
      F77NAME(zheevr)(&jobz, &range, &uplo, &n, zm, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, eigval, 
              z, &ldz, isuppz, zwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);
      free(isuppz); free(zwork); free(rwork); free(iwork); free(zm);
   }

   double elem; int ii; ind.reserve(n); for(i=0;i<n;i++) ind.push_back(i); i=1; j=2; 
   while(i<n)
   {
      if(eigval[i-1]<=eigval[i]) { i=j; j++; }
      else { elem = eigval[i-1]; eigval[i-1] = eigval[i]; eigval[i] = elem; ii = ind[i-1]; ind[i-1] = ind[i]; ind[i] = ii; i--; if(i==0) i=1; }
   }
                           
   E.reserve(nev); ind.reserve(nev); dt = eigval[0];
   for(i=0; i<nev; i++) 
   {
      E.push_back(eigval[i]);
      if(eigval[i]<dt) dt = eigval[i];
   }
   Esz = E.size(); free(eigval); 

   // Sets energy levels relative to lowest level.
   for(i=0; i<Esz; i++) E[i] -= dt;                                   // E = E - min(E);

   // Converts energy levels from cm^{-1} to J - NB. E = hf = hc/lambda in metres!
   for(i=0; i<Esz; i++) E[i] *= h*c*100.;                             // E = E .* (h*c*100);

   for(ind_j=0; ind_j<Esz; ind_j++)                                   // for ind_j = 1:length(E)
   {  // Calculates the matrix elements <Vi|J.H|Vi>
      if(iHmltn.isempty())
      {
         vt = (double*)malloc(Hsz*sizeof(double)); 
         F77NAME(dsymv)(&uplo, &Hsz, &alpha, fJmat, &Hsz, &m[ind[ind_j]*Hsz], &incx, &beta, vt, &incx);
#ifdef _G77 
         F77NAME(ddot)(me, &Hsz, &m[ind[ind_j]*Hsz], &incx, vt, &incx);
#else	 
         me = F77NAME(ddot)(&Hsz, &m[ind[ind_j]*Hsz], &incx, vt, &incx);
#endif
      }
      else
      {
         zt = (complexdouble*)malloc(Hsz*sizeof(complexdouble));
         F77NAME(zhemv)(&uplo, &Hsz, &zalpha, zJmat, &Hsz, &z[ind[ind_j]*Hsz], &incx, &zbeta, zt, &incx);
#ifdef _G77 
         F77NAME(zdotc)(&zme, &Hsz, &z[ind[ind_j]*Hsz], &incx, zt, &incx);
#else	 
         zme = F77NAME(zdotc)(&Hsz, &z[ind[ind_j]*Hsz], &incx, zt, &incx);
#endif
         me = zme.r;                                                  //   me = V(:,ind_j)' * Jmat * V(:,ind_j);
      }

      if(iHic.isempty() && iJmat.isempty()) free(vt); else free(zt);

      // Calculates the elements of the expectation of J.H: <i|J.H|i> exp(-Ei(H)/kT) and the partition function: exp(-Ei(H)/kT)
      for(i=0; i<(int)T.size(); i++)                                  //   JH(:,ind_j) = (me - sum(me_ndegen)) .* exp(-beta .* E(ind_j));
      {                                                               //   Z(:,ind_j) = exp(-beta .* E(ind_j));
         dt2 = exp(-E[ind_j]/(k_B*T[i])); mag[i]+=me*dt2; Z[i]+=dt2;
      }
   }                                                                  // end

   if(iHic.isempty() && iJmat.isempty()) { free(m); free(fJmat); } else { free(z); free(zJmat); }

   // Calculates the expectation <<J.H>> = sum(<i|J.H|i>exp(-Ei/kT)) / sum(exp(-Ei/kT))
   for(i=0; i<(int)T.size(); i++) mag[i] /= Z[i];
                                                                      // for ind_T = 1:size(T,2)
   return mag;                                                        //   exp_JH(ind_H,ind_T) = sum(JH(ind_T,:)) / sum(Z(ind_T,:));
                                                                      // end

   // Calculates the magnetisation M(ind_H,ind_T) per unit volume per atom;
   // mag = exp_JH;  % in u_B/atom
   // %M = g * mu_B * exp_JH .* N_A; % J/T/m^3/mol = A/m/mol
   // To get magnetisation in emu/g (cgs units): Get magnetisation in J/T/m^3/kg by
   // %  M = g*exp_JH * mu_B * N_A/(molar mass). Then multiply by 4pi*10^{-7}.
   // To get susceptibility in emu/mol (cgs units): Take M in J/T/m^3/mol, divide by
   // %  field in Tesla, and divide result by 4pi
   // NB. H = B/mu_0 = B/(4pi*10^{-7}) where B is in Tesla. 
   //   and chi = M/H not M/B => chi = M/(field in Tesla)*(4pi*10^{-7}) in m^3/mol.
   // NB. 1T = 1Netwon/A/m 
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: 
//    g++ states.cpp cfp.cpp njsyms.cpp so_cf.cpp coulomb.cpp ic_hmltn.cpp && ./a.out
/* --------------------------------------------------------------------------------------------------------------- //
#include <ctime>
int main(int argc, char *argv[])
{
   clock_t start,mid,end; start = clock();

   int n;
   if(argc>1) { n = atoi(argv[1]); } else { n = 2; }

   std::vector<double> F(4,1.);  // F = [1 1 1 1];
   std::vector<double> B2(5,0.);  
   std::vector<double> B4(9,0.);
   std::vector<double> B6(13,0.);
   std::vector< std::vector<double> > B;
   sMat<double> Hic,Hici;

   B4[4] = 1.; B4[8] = sqrt(5./14.); B6[6] = 1.; B6[10] = -sqrt(7./2.);
   B.push_back(B2); B.push_back(B4); B.push_back(B6);
   
 //Hic = racah_ukq(n,2,0); rmzeros(Hic); mm_gout(Hic,"n2U20.mm"); 
 //sMat<double> Hic2 = mm_gin("n2U20.mm"); mm_gout(Hic2,"n2U20a.mm"); // Test with $ diff n2U20*mm

   Hic = ic_hmltn(Hici,n,F,1.,B);      // ic_hmltn with racah_ukq : real    19m14.121s : fast_ukq : 11m33.990s
 //Hic = racah_ukq(n,4,0);             //   n=7                     user    11m18.798s              11m9.790s
 //Hic = racah_emat(n,1.,1.,1.,1.);    //                            sys     0m17.925s              0m5.996s
 //Hic = racah_so(n,1.);               // fast_ukq doesn't hit the cache (much - about 10%), unlike racah_ukq
 //std::cout << "x = " <<  Hic.display_full(); 
   // To compare with matlab: ./ic_hmltn > out.m && mv out.m ~/work/matlab/saficf/racah
   // x=[]; out; y=full(ic_hmltn('f2',[1 1 1 1],1,{zeros(1,5) [0 0 0 0 1 0 0 0 sqrt(5/14)] ...
   //       [0 0 0 0 0 0 1 0 0 0 -sqrt(7/2) 0 0]})); y(abs(y)<1e-4)=0; sum(sum(abs(x-y)))
 //std::cout << "x = sparse(" << Hic.nr() << "," << Hic.nc() << ");\n" << Hic;
 //std::cout << "ic_hmltn(" << n << ",[1 1 1 1],1,{B4=1,B6=1}) =\n" <<  Hic;
 //sMat<double> Hicconj = Hic; Hicconj.transpose();
 //sMat<double> comp = Hic - Hicconj;
 //std::cout << "sum(sum(abs( Hic-Hic.transpose() ))) = " << vsum(msum(mabs(comp))) << "\n"; // Checks H_IC is hermitian
 
 //rmzeros(Hic); mm_sout(Hic,"Hic.mm"); Hic2 = mm_sin("Hic.mm"); mm_gout(Hic2,"Hic2.mm"); // Test with $ diff Hic*mm

   // Test the zeeman term
 //Hic += racah_mumat(n,0);
   mid = clock();

   // Using ARPACK++ to solve the IC Hamiltonian for a few eigenvalues
 //char whichp[] = "LM"; int nev = 18;
 //int ncvp = (2*nev+1)>(Hsz-1) ? (2*nev+1) : Hsz-1; 
 //ARSymStdEig<double, sMat<double> > dprob(Hsz, nev, &Hmltn, &sMat<double>::MultMv, whichp, ncvp);
 //dprob.FindEigenvectors();
 //int i;
 //for(i=0; i<nev; i++)
 //   std::cout << "E[" << i << "] = " << dprob.Eigenvalue(i) << "\n";
 
   // Testing ic_mag()
 //int i;
 //std::vector<double> T(50,0.); for(i=0; i<50; i++) T[i] = i+1;
 //sMat<double> Jmat = racah_mumat(n,0);
 //std::vector<double> mag = ic_mag(Hic,Jmat,T);
 //for(i=0; i<50; i++) std::cout << T[i] << "\t" << mag[i] << "\n"; 
 
   // Full diagonalisation with internal routine
   eigVE<double> dg = eig(Hic); 
   std::cout << "E = " << dispvect(dg.E) << "\nx=[];" << dg.V << "\nV = x;\n";

   end = clock(); std::cerr << "Time to calculate Hic = " << (double)(mid-start)/CLOCKS_PER_SEC << "s. Time to diagonalise = " << (double)(end-mid)/CLOCKS_PER_SEC << "s.\n";

   return 0;
}*/
