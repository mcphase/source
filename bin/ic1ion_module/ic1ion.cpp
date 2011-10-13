/* ic1ion.cpp
 *
 * Main routine for the stand-alone program, handles input/output, and calculates single-ion magnetisation, 
 * energy levels and wavefunctions.
 *
 * Functions:
 *   void ic_printheader(const char *filename, icpars &pars);                  // Prints header to file
 *   void ic_showoutput(const char *filename, icpars &pars);                   // Prints calculated spectra to file
 *   void ic_showoutput(const char*file, eigVE<double>&d, fconf&f, icpars&p);  // As above but with sMat<> format
 *   void ic_cmag(const char *filename, icpars &pars);                         // Calcs. magnetisation using icmfmat::
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include <fstream>
#include <iomanip>
#include <ctime>


#ifdef _INTEGRAL
#include "../include/vector.h"
extern "C" void mcalc(Vector &J, double *T, Vector &gjmbH, double *gJ, Vector &ABC, char **sipffilename, double *lnZ, double *U, ComplexMatrix &est);
extern "C" int estates(ComplexMatrix &est, Vector &gjmbheff, double *gJ, double *T, Vector &ABC, char **sipffilename);
extern "C" int du1calc(int &tn, double &T, Vector &gjmbH, double &g_J, Vector &ABC, char **sipffilename, ComplexMatrix &mat, float &delta, ComplexMatrix &est);
extern "C" int mq(ComplexVector &Mq, double &th, double &ph, double &J0, double &J2, double &J4, double &J6, ComplexMatrix &est);
extern "C" int dv1calc(int &tn, double &th, double &ph, double &J0, double &J2, double &J4, double &J6, ComplexMatrix &est, double &T, ComplexMatrix &mat);
extern "C" void mcalc_parameter_storage_matrix_init(ComplexMatrix *est, Vector &gjmbheff, double *g_J, double *T, Vector &ABC, char **sipffilename);
#endif

// --------------------------------------------------------------------------------------------------------------- //
// Converts eigenvectors to different basis 
// --------------------------------------------------------------------------------------------------------------- //
void ic_conv_basis(icpars &pars, iceig &VE, fconf &conf)
{
   #define CFS conf.states
   #define CJS confJ.states
   std::string basis; basis.assign(pars.basis); strtolower(basis);
   if(basis.find("msml")!=std::string::npos)
   {
      // Enumerate the states in the |LmL,SmS> basis
      fconf confLS(pars.n,pars.l), confJ(pars.n,1,pars.l); int L2,S2,ns=0; char Jlabel[12];
      for(int ii=0; ii<(int)confLS.states.size(); ii++)
      {
         L2 = 2*abs(confLS.states[ii].L); S2 = confLS.states[ii].S2;
         for(int mL2=-L2; mL2<=L2; mL2+=2)
         {
            for(int mS2=-S2; mS2<=S2; mS2+=2)
            {
               conf.states[ns] = confLS.states[ii]; conf.states[ns].J2 = mL2; conf.states[ns].mJ2 = mS2; 
               conf.states[ns].id.assign(confLS.states[ii].id);
               if(mL2%2==0) sprintf(Jlabel,",mL=%i",mL2/2); else sprintf(Jlabel,",mL=%i/2",mL2); conf.states[ns].id.append(Jlabel);
               if(mS2%2==0) sprintf(Jlabel,",mS=%i",mS2/2); else sprintf(Jlabel,",mS=%i/2",mS2); conf.states[ns].id.append(Jlabel);
//std::cout << "State: " << ns+1 << "\t";
//std::cout << "U=" << CFS[ns].U << ",v=" << CFS[ns].v << ",L=" << CFS[ns].L << ",S2=" << CFS[ns].S2 << ",mL=" << CFS[ns].J2/2. << ",mS=" << CFS[ns].mJ2/2. << "  \tid=" << CFS[ns].id << "\t";
//std::cout << "U=" << CJS[ns].U << ",v=" << CJS[ns].v << ",L=" << CJS[ns].L << ",S2=" << CJS[ns].S2 << ",J="  << CJS[ns].J2/2. << ",mJ=" << CFS[ns].mJ2/2. << "\tid=" << CJS[ns].id << "\n";
               ns++;
            }
         }
      }
      std::cout << "ic_conv_basis(): Converting from |LSmJ> to |LmL,SmS> basis.\n";
      std::cout << "ic_conv_basis(): States Check. Number of LS states: " << ns << ", Number of mJ states " << (int)confJ.states.size() << "\n";     
      char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
      if(pars.save_matrices) {
      #ifndef _WINDOWS
      struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
         if(mkdir("results/mms",0777)!=0) std::cerr << "ic_conv_basis(): Can't create mms dir, " << strerror(errno) << "\n";
      #else
      DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
         if (!CreateDirectory("results\\mms", NULL)) std::cerr << "ic_conv_basis(): Cannot create mms directory\n";
      #endif 
      nstr[0] = (pars.l==F?102:100); if(pars.n<10) { nstr[1] = pars.n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = pars.n+38; nstr[3] = 0; }
      strcat(basename,nstr); strcpy(filename,basename); strcat(filename,"_JmJ2mSmL.mm");
      } else { strcpy(filename,"nodir/nofile"); }
      sMat<double> convmat; convmat = mm_gin(filename); int J2,mL2,mS2,mJ2;
      if(convmat.isempty())    // Conversion matrix not previously saved... Needs to be calculated
      {
         convmat.zero(ns,ns);
         for(int i=0; i<ns; i++)
            for(int j=0; j<ns; j++)
            {
               if(CFS[i].L==CJS[j].L && CFS[i].S2==CJS[j].S2 && CFS[i].U==CJS[j].U && CFS[i].v==CJS[j].v)
               {
                  L2=2*abs(CFS[i].L); S2=CFS[i].S2; J2=CJS[j].J2; mL2=CFS[i].J2; mS2=CFS[i].mJ2; mJ2=CJS[j].mJ2;
                  convmat(i,j) = sqrt(J2+1.) * threej(L2,S2,J2,mL2,mS2,mJ2); if((L2-S2+mJ2)%4==2) convmat(i,j)=-convmat(i,j);
               }
            }
	 rmzeros(convmat); mm_gout(convmat,filename);
      }
      char /*transpose = 'C',*/ notranspose = 'N';
      if(VE.iscomplex())
      {
         sMat<double> zeros(ns,ns); complexdouble *zmt,zalpha,zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0;
         complexdouble *zConv = zmat2f(convmat,zeros); zmt = new complexdouble[ns*ns];
         F77NAME(zgemm)(&notranspose,&notranspose,&ns,&ns,&ns,&zalpha,zConv,&ns,VE.zV(0),&ns,&zbeta,zmt,&ns);
       //F77NAME(zgemm)(&notranspose,&notranspose,&ns,&ns,&ns,&zalpha,VE.zV(0),&ns,zConv,&ns,&zbeta,zmt,&ns);
         memcpy(VE.zV(0),zmt,ns*ns*sizeof(complexdouble)); delete[]zmt; free(zConv);
      }
      else
      {
         double *dConv = convmat.f_array(); double alpha=1., beta=0.; char notranspose = 'N'; double *dmt = new double[ns*ns];
         F77NAME(dgemm)(&notranspose, &notranspose, &ns, &ns, &ns, &alpha, dConv, &ns, VE.V(0), &ns, &beta, dmt, &ns);
       //F77NAME(dgemm)(&notranspose, &notranspose, &ns, &ns, &ns, &alpha, VE.V(0), &ns, dConv, &ns, &beta, dmt, &ns);
         memcpy(VE.V(0),dmt,ns*ns*sizeof(double)); delete[]dmt; free(dConv);
      }
      // Checks that the eigenvalues are orthonormal
      char transa='C', transb='N'; double summm=0., alpha=1., beta=0.; complexdouble zalpha,zbeta; zalpha.r=1; zalpha.i=0; zbeta.r=0; zbeta.i=0; int incx=1;
      if(VE.iscomplex())
      {
         complexdouble *zmm = (complexdouble*)malloc(ns*ns*sizeof(complexdouble)); 
         complexdouble *vet = (complexdouble*)malloc(ns*ns*sizeof(complexdouble)); memcpy(vet,VE.zV(0),ns*ns*sizeof(complexdouble));
         F77NAME(zgemm)(&transa, &transb, &ns, &ns, &ns, &zalpha, vet, &ns, VE.zV(0), &ns, &zbeta, zmm, &ns);
         for(int ii=0; ii<ns; ii++) { zmm[ii*ns+ii].r-=1.; summm += F77NAME(dzasum)(&ns, &zmm[ii*ns], &incx); if(VE.E(ii+1)==0) break; }
         std::cout << "ic_conv_basis(): Orthonomality Test. Sum(V^TV-I) = " << summm << "\n";
         free(zmm); free(vet);
      }
      else
      {
         double *dmm = (double*)malloc(ns*ns*sizeof(double)); 
         double *vet = (double*)malloc(ns*ns*sizeof(double)); memcpy(vet,VE.V(0),ns*ns*sizeof(double));
         F77NAME(dgemm)(&transa, &transb, &ns, &ns, &ns, &alpha, vet, &ns, VE.V(0), &ns, &beta, dmm, &ns);
         for(int ii=0; ii<ns; ii++) { dmm[ii*ns+ii]-=1.; summm += F77NAME(dasum)(&ns, &dmm[ii*ns], &incx); if(VE.E(ii+1)==0) break; }
         std::cout << "ic_conv_basis(): Orthonomality Test. Sum(V^TV-I) = " << summm << "\n";
         free(dmm); free(vet);
      }
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Prints out a header to a specified file
// --------------------------------------------------------------------------------------------------------------- //
void ic_printheader(const char *outfile, icpars &pars)
{
   time_t curtime = time(NULL);
   std::fstream FILEOUT; FILEOUT.open(outfile, std::fstream::out);
   std::string Lstr = Lout(pars.l); strtolower(Lstr);
   FILEOUT << "#{ ic1ionmodule version " << IC1IONMODULE_VERSION << " " << ctime(&curtime);
   if(!pars.ionname.empty()) FILEOUT << "# Ion name: " << pars.ionname << "\n";
   FILEOUT << "# Free ion configuration: " << Lstr << "^" << pars.n << "\n";
   FILEOUT << "# Free ion parameters (" << pars.e_units << "): F^2=" << pars.F[1] << " F^4=" << pars.F[2];
   if(pars.l==F) FILEOUT << " F^6=" << pars.F[3]; FILEOUT << " zeta=" << pars.xi << " alpha=" << pars.alpha[0] << " beta=" << pars.alpha[1]; 
   if(pars.l==D) FILEOUT << "\n"; else FILEOUT << " gamma=" << pars.alpha[2] << "\n";
   FILEOUT << "# Crystal Field parameters normalisation: " << pars.B.norm() << "\n";
   std::string norm=pars.B.norm(); strtolower(norm); if(norm.find("stev")!=std::string::npos)
   {
      std::string op; if(pars.B.op_equiv==Lt) op.assign("L"); else op.assign("J");
      FILEOUT << "# Stevens Factors: <" << op << "||alpha||" << op << ">=" << pars.B.alpha();
      FILEOUT <<                  ", <" << op << "||beta||" << op << ">=" << pars.B.beta();
      if(pars.l==F) FILEOUT << ", <" << op << "||gamma||" << op << ">=" << pars.B.gamma(); FILEOUT << "\n";
   }
   FILEOUT << "# Crystal Field parameters (" << pars.B.units() << "): " << pars.B.cfparsout(", ") << "\n";
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      FILEOUT << "# With magnetic field: Bx=" << pars.Bx << ", By=" << pars.By << ", Bz=" << pars.Bz << " Tesla.\n";
   }
   FILEOUT.close();
}
 
// --------------------------------------------------------------------------------------------------------------- //
// Outputs the energy and wavefunctions to a specified file
// --------------------------------------------------------------------------------------------------------------- //
void ic_showoutput(const char *filename,                        // Output file name - default "results/mcphas.icr"
                   icpars &pars,                                // Input parameters
                   iceig &VE,                                   // Eigenstates class
                   int iconf,
                   std::vector<int> ikeepJ)
{
   fconf conf(pars.n,iconf,pars.l);

   unsigned int iE,iV,i=1,j=2,num_states=conf.states.size(),ns;
   std::vector<int> isE; isE.reserve(num_states);
   int ii; 
   double elem,conv=1.; complexdouble elc;

   if(pars.e_units.find("meV")!=std::string::npos) conv = 1./MEV2CM; 
   else if(pars.e_units.find("K")!=std::string::npos) conv = CM2K;

   if(iconf==1) ic_conv_basis(pars,VE,conf);

   ic_printheader(filename,pars);
   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out | std::fstream::app); // Opens file for appending
   FILEOUT << "# Energy offset, E0=" << VE.E(0)*conv << pars.e_units << "\n";
   if(!VE.iscomplex()) FILEOUT << "# Energy(" << pars.e_units << ")\tWavefunctions(^{2S+1}L_J,mJ) }\n"; else
   FILEOUT << "# Energy(" << pars.e_units << ")\tAmplitude\t|Amplitude|^2\tWavefunctions(^{2S+1}L_J,mJ) }\n";

   bool istrunc = !ikeepJ.empty(); if(istrunc) num_states=ikeepJ.size();
   std::vector<int> isV(num_states,0), _isV; _isV.reserve(num_states);
   if(istrunc) _isV = ikeepJ; else for(ii=0; ii<(int)num_states; ii++) _isV.push_back(ii);

   double *V=0; complexdouble *zV=0; if(VE.iscomplex()) zV = new complexdouble[num_states]; else V = new double[num_states];
   if((unsigned int)pars.num_eigv > num_states) { ns = num_states; } else { ns = (unsigned int)pars.num_eigv; }
   for(iE=0; iE<num_states; iE++)
   {
      if(VE.E(iE)==-DBL_MAX) { FILEOUT << "# Higher levels truncated\n"; break; }      // For case of truncated matrices.
      FILEOUT << (VE.E(iE)-VE.E(0))*conv << "\t\t";
      isV = _isV; i=1; j=2;
      if(VE.iscomplex())
      {
         memcpy(zV,VE.zV(iE),num_states*sizeof(complexdouble));
         while(i<num_states)
         {
            if((zV[i-1].r*zV[i-1].r+zV[i-1].i*zV[i-1].i) >= (zV[i].r*zV[i].r+zV[i].i*zV[i].i)) { i=j; j++; }
            else { elc = zV[i-1]; zV[i-1] = zV[i]; zV[i] = elc; ii=isV[i-1]; isV[i-1]=isV[i]; isV[i]=ii; i--; if(i==0) i=1; }
         }
         elc = zV[0];
         FILEOUT << " (" << elc.r; if(elc.i>0) FILEOUT << "+"; else FILEOUT << "-";
         FILEOUT << "i" << fabs(elc.i) << ")\t\t" << (elc.r*elc.r+elc.i*elc.i) << "\t";
         FILEOUT << "|" << conf.states[isV[0]].id << ">";
         for(iV=1; iV<ns; iV++)
         {
            elc = zV[iV]; FILEOUT << "\n\t\t+";
            FILEOUT << "(" << elc.r; if(elc.i>0) FILEOUT << "+"; else FILEOUT << "-";
	    FILEOUT << "i" << fabs(elc.i) << ")\t\t" << (elc.r*elc.r+elc.i*elc.i) << "\t";
            FILEOUT << "|" << conf.states[isV[iV]].id << ">";
         }
         FILEOUT << "\n";
      }
      else
      {
         memcpy(V,VE.V(iE),num_states*sizeof(double));
         while(i<num_states)
         {
            if(fabs(V[i-1])>=fabs(V[i])) { i=j; j++; }
            else { elem = V[i-1]; V[i-1] = V[i]; V[i] = elem; ii=isV[i-1]; isV[i-1]=isV[i]; isV[i]=ii; i--; if(i==0) i=1; }
         }
         FILEOUT << V[0] << "|" << conf.states[isV[0]].id << ">\t";
         for(iV=1; iV<ns; iV++)
         {
            if(iV%4==0) FILEOUT << "\n\t\t";
            if(V[iV]>0) FILEOUT << "+";
            FILEOUT << V[iV] << "|" << conf.states[isV[iV]].id << ">\t";
         }
         FILEOUT << "\n";
      }
   }
   if(VE.iscomplex()) delete[]zV; else delete[]V;
   FILEOUT.close();
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the magnetisation
// --------------------------------------------------------------------------------------------------------------- //
void ic_cmag(const char *filename, icpars &pars)
{
   if(!(pars.calcphys & PHYSPROP_MAGBIT)) return;
   icmfmat mfmat(pars.n,pars.l,6,pars.save_matrices);

   ic_printheader(filename,pars);
   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out | std::fstream::app); // Opens file for appending

   if(pars.xT==0. && pars.yT==0.) { std::cerr << "ic_showphysprop(): Either x- or y-axis must be temperature.\n"; return; }

   double Tmin, Tstep, Tmax, Hmin, Hstep, Hmax;
   if(pars.xT!=0.) { Tmin = pars.xMin; Tstep = pars.xStep; Tmax = pars.xMax; Hmin = pars.yMin; Hstep = pars.yStep; Hmax = pars.yMax; }
              else { Tmin = pars.yMin; Tstep = pars.yStep; Tmax = pars.yMax; Hmin = pars.xMin; Hstep = pars.xStep; Hmax = pars.xMax; }

   if(pars.xT!=0.) FILEOUT << "# T(K)\tHa(T)\tHb(T)\tHc(T)\t"; else FILEOUT << "# Ha(T)\tHb(T)\tHc(T)\tT(K)\t";
   if(pars.mag_units==0) FILEOUT << "Magnetisation(uB/atom)\tMa\tMb\tMc\tM_parallel\n";
   else if(pars.mag_units==1) FILEOUT << "Magnetisation(emu/mol)\tMa\tMb\tMc\tM_parallel\n";
   else if(pars.mag_units==2) FILEOUT << "Magnetisation(Am^2/mol)\tMa\tMb\tMc\tM_parallel\n";
   int i,j,k; double Hm=0.,dt,Z;
   int nT = (int)ceil((Tmax-Tmin)/Tstep), nH = (int)ceil((Hmax-Hmin)/Hstep) + 1;
   std::vector<double> T(nT,0.); for(i=0; i<nT; i++) T[i] = Tmin+i*Tstep;
   std::vector<double> mag(nT,0.),ma(nT,0.),mb(nT,0.),mc(nT,0.);

   if((pars.xHa!=0. && pars.yHa!=0.) || (pars.xHb!=0. && pars.yHb!=0.) || (pars.xHc!=0. && pars.yHc!=0.))
   {
      std::cerr << "ic_showphysprop(): You may only specify magnetic field on either x- and y-axis - not both.\n"; return;
   }
   double xnorm = sqrt(pars.xHa*pars.xHa+pars.xHb*pars.xHb+pars.xHc*pars.xHc); if(xnorm==0) xnorm=1.;
   double ynorm = sqrt(pars.yHa*pars.yHa+pars.yHb*pars.yHb+pars.yHc*pars.yHc); if(ynorm==0) ynorm=1.;
   std::vector<double> gjmbH(6,0.), gjmbHmeV(6,0.); 
   if(pars.xT==0.) gjmbH[1]=pars.xHa/xnorm; else gjmbH[1]=pars.yHa/ynorm; gjmbH[0]=GS*gjmbH[1];
   if(pars.xT==0.) gjmbH[3]=pars.xHb/xnorm; else gjmbH[3]=pars.yHb/ynorm; gjmbH[2]=GS*gjmbH[3];
   if(pars.xT==0.) gjmbH[5]=pars.xHc/xnorm; else gjmbH[5]=pars.yHc/ynorm; gjmbH[4]=GS*gjmbH[5];

   iceig VE;
   std::vector<double> ex; std::vector< std::vector<double> > matel, exj;
   sMat<double> J,iJ,H,iH; H = ic_hmltn(iH,pars); H/=MEV2CM; iH/=MEV2CM;

   double convfact=1; if(pars.mag_units==1) convfact = NAMUB*1e3; else if(pars.mag_units==2) convfact = NAMUB;  // 1==cgs, 2==SI

   for(i=0; i<nH; i++)
   {
      for(j=0; j<6; j++) gjmbHmeV[j] = gjmbH[j]*(-MUB*(Hmin+i*Hstep)); 
      mfmat.Jmat(J,iJ,gjmbHmeV,pars.save_matrices); J+=H; iJ+=iH; 
      if(pars.partial) VE.lcalc(pars,J,iJ); 
      #ifndef NO_ARPACK
         else if(pars.arnoldi) VE.acalc(pars,J,iJ); else VE.calc(J,iJ); 
      #endif
      ex = mfmat.expJ(VE,Tmax,matel,pars.save_matrices); ex.assign(matel[0].size(),0.); for(j=0; j<6; j+=2) exj.push_back(ex);
      for(j=0; j<nT; j++) 
      {
         ma[j] = 0.; mb[j] = 0.; mc[j] = 0.;  Z = 0.;
         for(k=0; k<(int)matel[0].size(); k++) 
         { 
            if(j==0) for(int ii=0; ii<6; ii+=2) exj[ii/2][k] = GS*matel[ii][k]+matel[ii+1][k]; 
            dt = exp(-(VE.E(k)-VE.E(0))/(KB*T[j])); Z+=dt;
            ma[j] += exj[0][k] * dt; mb[j] += exj[1][k] * dt; mc[j] += exj[2][k] * dt;
         }
         ma[j]/=Z; mb[j]/=Z; mc[j]/=Z; mag[j] = sqrt(ma[j]*ma[j]+mb[j]*mb[j]+mc[j]*mc[j]);
      }
      if(pars.xT!=0.)
      {
         Hm = (Hmin+i*Hstep)/ynorm;
         for(j=0; j<nT; j++) 
            FILEOUT << T[j] << "\t" << Hm*pars.yHa << "\t" << Hm*pars.yHb << "\t" << Hm*pars.yHc << "\t" << mag[j]*convfact
                    << "   \t" << ma[j]*convfact << "   \t" << mb[j]*convfact << "   \t" << mc[j]*convfact << "   \t"
                    << (ma[j]*gjmbH[1] + mb[j]*gjmbH[3] + mc[j]*gjmbH[5])*convfact << "\n";
      }
      else
      {
         Hm = (Hmin+i*Hstep)/xnorm;
         for(j=0; j<nT; j++) 
            FILEOUT << Hm*pars.xHa << "\t" << Hm*pars.xHb << "\t" << Hm*pars.xHc << "\t" << T[j] << "\t" << mag[j]*convfact
                    << "   \t" << ma[j]*convfact << "   \t" << mb[j]*convfact << "   \t" << mc[j]*convfact << "   \t" 
                    << (ma[j]*gjmbH[1] + mb[j]*gjmbH[3] + mc[j]*gjmbH[5])*convfact << "\n";
      }
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Reads in parameters from mcphas.ic or a file specified on the command line.
// Calculates the IC Hamilton matrix; solve its eigensystem and outputs the energies and wavefunctions calculated.
// Also calculates the magnetisation and heat capacity if desired (specify iccalcmag=1,2 and/or iccalcCp=1)
// --------------------------------------------------------------------------------------------------------------- //
int main(int argc, char *argv[])
{
   char infile[255], outfile[255], physfile[255];
   icpars pars;
   std::string norm,units;

   bool headonly=false, hcsoonly=false, truncate=false; int ai[]={0,0,0,0}, ib=1, narg=(argc>4?4:argc); double elim=2e3;
   for(int ia=1; ia<narg; ia++) {
      if(strncmp(argv[ia],"-h",2)==0) headonly=true; else
      if(strncmp(argv[ia],"-t",2)==0) { truncate=true; elim = atoi(argv[++ia]); } else
      if(strncmp(argv[ia],"-s",2)==0) hcsoonly=true; else ai[ib++]=ia; }

   if(ai[1]!=0) strcpy(infile,  argv[ai[1]]); else strcpy(infile,"mcphas.ic");
   if(ai[2]!=0) strcpy(outfile, argv[ai[2]]); else strcpy(outfile,"results/ic1ion.out");
   if(ai[3]!=0) strcpy(physfile,argv[ai[3]]); else strcpy(physfile,"results/ic1ion.mag");

   // Gets input parameters and what to calculate from input single-ion parameters file
   ic_parseinput(infile,pars);

   if(headonly) { std::cout << "ic1ion: Outputing only header\n"; ic_printheader(outfile,pars); return(0); }
   if(hcsoonly) { iceig VE; sMat<double> Hcso = ic_Hcso(pars); VE.calc(Hcso); strcpy(outfile,"results/mcphas.icpJ"); ic_showoutput(outfile,pars,VE,0); return(0); }
   if(truncate) { iceig VE = spectre_eig(pars, elim); ic_showoutput(outfile,pars,VE,1); return(0); }

   // For quick and dirty timing routines...
   clock_t start,end; start = clock();

   // Calculates the intermediate coupling Hamilton matrix
   sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars);
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      std::vector<double> gjmbH(6,0.);
      if(fabs(pars.Bx)>DBL_EPSILON) { gjmbH[1]=-MUBc*pars.Bx; gjmbH[0]=GS*gjmbH[1]; }
      if(fabs(pars.By)>DBL_EPSILON) { gjmbH[3]=-MUBc*pars.By; gjmbH[2]=GS*gjmbH[3]; }
      if(fabs(pars.Bz)>DBL_EPSILON) { gjmbH[5]=-MUBc*pars.Bz; gjmbH[4]=GS*gjmbH[5]; }
      sMat<double> J,iJ; icmfmat mfmat(pars.n,pars.l,6,pars.save_matrices); mfmat.Jmat(J,iJ,gjmbH,pars.save_matrices); Hic+=J; iHic+=iJ;
   }

// std::cout << std::setprecision(16) << "Hic=" << Hic.display_full() << "; Hic=Hic./" << MEV2CM << ";\n";
// std::cout << std::setprecision(16) << "iHic=" << iHic.display_full() << "; iHic=iHic./" << MEV2CM << ";\n";

   end = clock(); std::cerr << "Time to calculate Hic = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";

   // Fully diagonalise IC Hamilton matrix and saves results to <outfile>
   iceig VE; 
   if(pars.partial_standalone)      if(iHic.isempty()) VE.lcalc(pars,Hic); else VE.lcalc(pars,Hic,iHic);
   #ifndef NO_ARPACK
   else if(pars.arnoldi_standalone) if(iHic.isempty()) VE.acalc(pars,Hic); else VE.acalc(pars,Hic,iHic);
   #endif
   else                             if(iHic.isempty()) VE.calc(Hic); else VE.calc(Hic,iHic);
   start = clock(); std::cerr << "Time to diagonalise = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
// iceig VE; VE.lcalc(pars,Hic,iHic); int i; for(i=0; i<6; i++) std::cout << (VE.E(i)-VE.E(0))/MEV2CM << " "; std::cout << "\n";
   ic_showoutput(outfile,pars,VE);
// sMat<double> Hcso = ic_Hcso(pars); VE.calc(Hcso); strcpy(outfile,"results/mcphas.icpJ"); ic_showoutput(outfile,pars,VE,0);
   end = clock(); std::cerr << "Time to save file = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";

   // If required calculates some single-ion physical properties and saves to <physfile>
   ic_cmag(physfile,pars);
   start = clock(); std::cerr << "Time to calculate and save magnetisation is = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";

// iceig VE(Hic,iHic);
/* fconf conf(pars.n,0,pars.l);
   int j,icv,imax=0; double vel,vmax; std::vector< std::vector<int> > cvSO2CF; std::vector<int> cvEl; icv=0;
   eigVE<double> VEcso = eig(Hcso); //std::cout << "VEcso=" << VEcso.V.display_full() << "\n";
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      cvEl.push_back(icv); icv += conf.states[imax].J2+1; cvEl.push_back(icv-1); cvSO2CF.push_back(cvEl); cvEl.clear();
   }
   sMat<double> Vcso = convH2H(VEcso.V,Hic.nr(),cvSO2CF); //std::cout << "Vcso=" << Vcso.display_full() << "\n";
   iceig VE2(Hic.nr(),En,Vcso.f_array());
   strcpy(outfile,"results/mcphas.icpmJ"); ic_showoutput(outfile,pars,VE2);
   delete []En; */
   
#ifdef _INTEGRAL
 //clock_t start,end; end = clock();
   Vector J(1,6,0.), gmbH(1,6,.0578838263), ABC; gmbH[1]*=2; gmbH[3]*=2; gmbH[5]*=2; 
   double T=2.0,lnZ=0.,U=0.,gJ=0.;
   char *filearray[1]; 
   filearray[0] = infile;
   ComplexMatrix est; mcalc_parameter_storage_matrix_init(&est,gmbH,&gJ,&T,ABC,filearray);
 //ComplexMatrix est; int Hsz=getdim(pars.n,pars.l); est = ComplexMatrix(0,Hsz,0,Hsz);
   end = clock();

   mcalc(J,&T,gmbH,&gJ,ABC,filearray,&lnZ,&U,est);
   start = clock(); std::cerr << "Time to do mcalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "lnZ = " << lnZ << ", U = " << U << "\n";
   std::cerr << "J[1] = " << J[1] << ", J[2] = " << J[2] << ", J[3] = " << J[3] << ", J[4] = " << J[4] << ", J[5] = " << J[5] << ", J[6] = " << J[6] << "\n";

//for (int it=0; it<1000; it++) {
// end = clock();
// mcalc(J,&T,gmbH,&gJ,ABC,filearray,&lnZ,&U,est);
// start = clock(); std::cerr << "Time to do mcalc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
// std::cerr << "lnZ = " << lnZ << ", U = " << U << "\n";
// std::cerr << "J[1] = " << J[1] << ", J[2] = " << J[2] << ", J[3] = " << J[3] << ", J[4] = " << J[4] << ", J[5] = " << J[5] << ", J[6] = " << J[6] << "\n";
// if(it%100==0) { std::cerr << it << " "; } } std::cerr << "\n";
 //gmbH[1] = 0.; gmbH[2] = 0.; gmbH[3] = 0.; gmbH[4] = 0.; gmbH[5] = 0.; gmbH[6] = 0.; 
   estates(est,gmbH,&gJ,&T,ABC,filearray);
   end = clock(); std::cerr << "Time to do estates() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   
   int imq, tn = 2; float delta=0.; ComplexMatrix mat6(1,6,1,6);
   imq = du1calc(tn,T,gmbH,gJ,ABC,filearray,mat6,delta,est);
   start = clock(); std::cerr << "Time to calculate du1calc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";

   ComplexVector Mq;
   double th=PI/4, ph=PI/4, J0=1., J2=1., J4=1., J6=1.;
 //double th=0., ph=0., J0=1., J2=1., J4=0., J6=0.;
   imq = mq(Mq,th,ph,J0,J2,J4,J6,est);
   end = clock(); std::cerr << "Time to calculate mq() = " << (double)(end-start)/CLOCKS_PER_SEC << "s.";
   std::cerr << " Mq = [" << Mq[1].real() << "+" << Mq[1].imag() << "i "
                          << Mq[2].real() << "+" << Mq[2].imag() << "i "
                          << Mq[3].real() << "+" << Mq[3].imag() << "i]\n";

   ComplexMatrix mat(1,6,1,6);
   imq = dv1calc(tn,th,ph,J0,J2,J4,J6,est,T,mat);
   start = clock(); std::cerr << "Time to calculate dv1calc() = " << (double)(start-end)/CLOCKS_PER_SEC << "s.\n";
#endif

   return 0;
}
