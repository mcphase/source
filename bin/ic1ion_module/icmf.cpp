/* icmf.cpp
 *
 * Holds classes which manipulate the CF and free ion parameters, and calculates the moment.
 *
 * Classes (in icpars.hpp)
 *   iceig   - Contains and calculates eigenvectors/values by various methods (full,partial,arnoldi)
 *   icmfmat - Contains the L and S matrices and calculates the magnetic moment
 *
 * Functions:
 *   void strtolower(std::string &instring);                           // Converts a string to lower case
 *   void conv_e_units(icpars &flags, std::string &newunits);          // Converts 1-ion pars to diff. energy units
 *   std::vector<double> stev_thetak(int n, orbital l);                // Calculates the Stevens factors.
 *   std::vector<double> rk_int(std::string &ionname);                 // Looks up value of radial integrals
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */
#include "ic1ion.hpp"
#include <cctype>                  // For std::tolower
#include <fstream>


// --------------------------------------------------------------------------------------------------------------- //
// Member function for complexdouble struct
// --------------------------------------------------------------------------------------------------------------- //
complexdouble complexdouble::operator=(const double v) { complexdouble t; t.r=v; t.i=0.; return t; }

// --------------------------------------------------------------------------------------------------------------- //
// Constructors for class iceig::
// --------------------------------------------------------------------------------------------------------------- //
iceig::iceig(int Hsz, bool isreal)
{
   _Hsz = Hsz; _E = new double[_Hsz]; 
   if(isreal) { _V = new double[_Hsz*_Hsz]; _zV = 0; } else { _zV = new complexdouble[_Hsz*_Hsz]; _V = 0; }
}
iceig::iceig(sMat<double>&H)
{
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; _zV = 0;
   int info = ic_diag(H,_V,_E); 
   if(info!=0) { std::cerr << "iceig(H) - Error diagonalising, info==" << info << "\n"; }
}
iceig::iceig(sMat<double>&H, sMat<double>&iH)
{
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; _V = 0;
   int info = ic_diag(H,iH,_zV,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; }
}
iceig::iceig(int Hsz, double *E, double *V)
{
   _Hsz = Hsz; _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; _zV = 0;
   memcpy(_E,E,Hsz*sizeof(double)); memcpy(_V,V,Hsz*Hsz*sizeof(double));
}
iceig::iceig(int Hsz, double *E, complexdouble *zV)
{
   _Hsz = Hsz; _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; _V = 0;
   memcpy(_E,E,Hsz*sizeof(double)); memcpy(_zV,zV,Hsz*Hsz*sizeof(complexdouble));
}
iceig::iceig(int Hsz, double *E, complexdouble *zV, int step)
{
   _Hsz = Hsz; _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; _V = 0;
   memcpy(_E,E,Hsz*sizeof(double)); 
   int i; for(i=0; i<Hsz; i++) memcpy(&_zV[i*Hsz],&zV[i*Hsz+i+step],Hsz*sizeof(complexdouble));
}
// --------------------------------------------------------------------------------------------------------------- //
// Copy constructors
// --------------------------------------------------------------------------------------------------------------- //
iceig::iceig(const iceig &p) { *this = p; }
iceig &iceig::operator = (const iceig &p) 
{ 
   if(_Hsz==p._Hsz)
   {
      if(_E==0) {_E = new double[_Hsz];} memcpy(_E,p._E,_Hsz*sizeof(double));
      if(p._V!=0) { if(_V==0) _V = new double[_Hsz*_Hsz]; memcpy(_V,p._V,_Hsz*_Hsz*sizeof(double)); }
      if(p._zV!=0) { if(_zV==0) _zV = new complexdouble[_Hsz*_Hsz]; memcpy(_zV,p._zV,_Hsz*_Hsz*sizeof(complexdouble)); }
      if(p._V==0 && _V!=0) { delete[]_V; _V=0; } if(p._zV==0 && _zV!=0) { delete[]_zV; _zV=0; } 
   }
   else
   {
      if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
      _Hsz = p._Hsz; _E = new double[_Hsz]; memcpy(_E,p._E,_Hsz*sizeof(double));
      if(p._V!=0) { _V = new double[_Hsz*_Hsz]; memcpy(_V,p._V,_Hsz*_Hsz*sizeof(double)); }
      if(p._zV!=0) { _zV = new complexdouble[_Hsz*_Hsz]; memcpy(_zV,p._zV,_Hsz*_Hsz*sizeof(complexdouble)); }
   }
   return *this; 
}
// --------------------------------------------------------------------------------------------------------------- //
// Destructor
// --------------------------------------------------------------------------------------------------------------- //
iceig::~iceig()
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
}
// --------------------------------------------------------------------------------------------------------------- //
// Methods for iceig:: class
// --------------------------------------------------------------------------------------------------------------- //
void iceig::calc(sMat<double>&H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz];
   int info = ic_diag(H,_V,_E); 
   if(info!=0) { std::cerr << "iceig(H) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
}
void iceig::calc(sMat<double>&H, sMat<double>&iH)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   int info = ic_diag(H,iH,_zV,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::calc(int Hsz, complexdouble *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = Hsz; _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   int info = ic_diag(Hsz,H,_zV,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::calc(int Hsz, double *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = Hsz; _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz];
   int info = ic_diag(H,Hsz,Hsz,_V,_E); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
}
void iceig::lcalc(icpars &pars, sMat<double>&H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_V,0,_Hsz*_Hsz*sizeof(double));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm
   }
   if(nev>_Hsz) nev=_Hsz;
   int info = ic_leig(H,_V,_E,nev); 
   if(info!=0) { std::cerr << "iceig(H) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
}
void iceig::lcalc(icpars &pars, sMat<double>&H, sMat<double>&iH)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   if(nev>_Hsz) nev=_Hsz;
   int info = ic_leig(H,iH,_zV,_E,nev); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
void iceig::lcalc(icpars &pars, complexdouble *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = getdim(pars.n,pars.l); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz];
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   if(nev>_Hsz) nev=_Hsz;
   int info = ic_leig(_Hsz,H,_zV,_E,nev); 
   if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
}
#ifndef NO_ARPACK
void iceig::acalc(icpars &pars, sMat<double>&H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _V = new double[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_V,0,_Hsz*_Hsz*sizeof(double));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm
   }
   if(nev>=_Hsz)   // We want all eigenvalues - better not to use the Arnoldi method
   {
      int info = ic_diag(H,_V,_E); 
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
   }
   else
   {
      double *dH = H.f_array(); int info = ic_arpackeig(_Hsz,dH,_V,_E,nev); free(dH);
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_V; _V=0; }
   }
}
void iceig::acalc(icpars &pars, sMat<double>&H, sMat<double>&iH)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = H.nc(); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm
   }
   if(nev>=_Hsz)   // We want all eigenvalues - better not to use the Arnoldi method
   {
      int info = ic_diag(H,iH,_zV,_E); 
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
   }
   else
   {
      complexdouble *zH = zmat2f(H,iH); int info = ic_arpackeig(_Hsz,zH,_zV,_E,nev); free(zH);
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; }
   }
}
void iceig::acalc(icpars &pars, complexdouble *H)
{
   if(_E!=0) { delete[]_E; _E=0; } if(_V!=0) { delete[]_V; _V=0; } if(_zV!=0) { delete[]_zV; _zV=0; }
   _Hsz = getdim(pars.n,pars.l); _E = new double[_Hsz]; _zV = new complexdouble[_Hsz*_Hsz]; 
   memset(_E,0,_Hsz*sizeof(double)); memset(_zV,0,_Hsz*_Hsz*sizeof(complexdouble));
   sMat<double> Hcso = ic_Hcso(pars); rmzeros(Hcso); eigVE<double> VEcso = eig(Hcso); 
   fconf conf(pars.n,0,pars.l); int i,j,imax=0,nev=0; double vel,vmax;
   for(i=0; i<Hcso.nr(); i++) 
   { 
      vmax = 0; for(j=0; j<Hcso.nr(); j++) { vel = fabs(VEcso.V(j,i)); if(vel>vmax) { vmax=vel; imax=j; } }
      nev += conf.states[imax].J2+1; if(exp(-(VEcso.E[i]-VEcso.E[0])/(208.510704))<DBL_EPSILON) break;   // 208.5==300K in 1/cm 
   }
   if(nev>=_Hsz)   // We want all eigenvalues - better not to use the Arnoldi method
   {
      int info = ic_diag(_Hsz,H,_zV,_E); 
      if(info!=0) { std::cerr << "iceig(H,iH) - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; exit(EXIT_FAILURE); }
   }
   else
   {
      int info = ic_arpackeig(_Hsz,H,_zV,_E,nev); 
      if(info!=0) { std::cerr << "iceig::acalc() - Error diagonalising, info==" << info << "\n"; delete[]_E; _E=0; delete[]_zV; _zV=0; exit(EXIT_FAILURE); }
   }
}
#endif
std::string iceig::strout() 
{
   int i,j;
   std::stringstream ss;
   for(i=0; i<_Hsz; i++) {ss << _E[i];} ss << "\n";
   for(i=0; i<_Hsz; i++) { for(j=0; j<_Hsz; j++) ss << _V[j*_Hsz+i]; ss << "\n"; }
   return ss.str();
}
// --------------------------------------------------------------------------------------------------------------- //
// Constructor for class icmfmat::
// --------------------------------------------------------------------------------------------------------------- //
icmfmat::icmfmat()
{ 
   sMat<double> t; J.assign(6,t); T.assign(8,NULL);
   iflag.assign(6,0); iflag[1]=1; iflag[4]=1;
   _n = 1; _l = S; _num_op = 1; _save_matrices=false;_xyz=0;
   #ifdef JIJCONV
   jijconv.assign(1,0);
   #endif
}
icmfmat::icmfmat(int n, orbital l, int num_op, bool save_matrices, int xyz)
{
   _n = n; _l = l; _num_op = num_op; _xyz = xyz;_save_matrices=save_matrices;
   sMat<double> t; J.assign(num_op>6?num_op:6,t); T.assign(num_op>6?num_op+2:8,NULL);
   iflag.assign(num_op>6?num_op:6,0); 
   if(xyz==0)for(int m=0;m<(num_op>6&&!save_matrices?num_op:6);++m)op_generate(m);
}

icmfmat::icmfmat(const icmfmat & pp)
{_n=pp._n;_l=pp._l;_num_op=pp._num_op;_xyz=pp._xyz;_save_matrices=pp._save_matrices;
 J=pp.J;T=pp.T;iflag=pp.iflag;
}

icmfmat::~icmfmat()
{
for(int i=0;i<(int)T.size();++i)if(T[i]!=NULL)delete []T[i];
}

complexdouble * icmfmat::balcar_Mq(int xyz, int K, int Q, int n, orbital l)
{bool imag;
 int Hsz = getdim(n,l);sMat<double> zeroes; zeroes.zero(Hsz,Hsz);
 if(_n!=n||_l!=l){std::cerr << "icmfmat::balcar_Mq() - n or l do not match!\n"; exit(EXIT_FAILURE); }
 sMat<double> retval=balcar_Mq(xyz,  K, Q, imag);
 if(imag)return zmat2f(zeroes,retval);else return zmat2f(retval,zeroes);
}

sMat<double>  icmfmat::balcar_Mq(int xyz, int K, int Q, bool & imag )
{ imag=false;
#define NSTRB(K,Q) nstr[3] = K+48; nstr[4] = Q+48; nstr[5] = 0
#define MSTRB(K,Q) nstr[3] = K+48; nstr[4] = 109;  nstr[5] = Q+48; nstr[6] = 0
#define NPOS std::string::npos
   int Hsz = getdim(_n,_l);
   sMat<double> retval(Hsz,Hsz),qpp,qmp,qpm,qmm;
   // Using the reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero... see balcar* in lovesey.cpp
   if((K%2==1) || (K>2*_l)) {  retval.zero(Hsz,Hsz); return retval; }
      
   char nstr[7]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   if(_save_matrices) {
  nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 77;   // ASCII codes: 77="M", 83=="S", 100=="d", 102=="f", 109=="m", 112="p"
   } else { strcpy(basename,"nodir/"); }

   if(xyz==1||xyz==2)
   {
      nstr[1]=83;
      if(Q!=0)
      {
         NSTRB(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MSq(1,K,abs(Q),_n,_l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTRB(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MSq(1,K,-abs(Q),_n,_l); rmzeros(qmp); mm_gout(qmp,filename); }
         NSTRB(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MSq(-1,K,abs(Q),_n,_l); rmzeros(qpm); mm_gout(qpm,filename); }
         MSTRB(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmm = mm_gin(filename); if(qmm.isempty()) { qmm = balcar_MSq(-1,K,-abs(Q),_n,_l); rmzeros(qmm); mm_gout(qmm,filename); }
	 // sum to coeff of Zlm (neglecting a 1/sqrt(2) factor) 
         if(Q<0) { if(Q%2==0) { qmp -= qpp; qmm -= qpm; } else { qmp += qpp; qmm += qpm; } }
         else    { if(Q%2==0) { qmp += qpp; qmm += qpm; } else { qmp -= qpp; qmm -= qpm; } }
         // add spherical components and multiply by addition factor 1/sqrt(2)which was neglected in the line above
         if(xyz==1) { if(Q<0) {imag=true;retval = (qmm-qmp)/(-2.);}    else retval = (qmm-qmp)/2.; }
         if(xyz==2) { if(Q<0) retval = (qmm+qmp)/2.; else {imag=true;retval = (qmm+qmp)/2.;} }// changed MR 25.5.2010 // Q<0 signs changed MR 28.5.2010
      }
      else
      {
         NSTRB(K,0); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MSq(1,K,0,_n,_l); rmzeros(qpp); mm_gout(qpp,filename); }
         NSTRB(K,0); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MSq(-1,K,0,_n,_l); rmzeros(qpm); mm_gout(qpm,filename); }      
         if(xyz==1) {  retval = (qpp-qpm)/(-sqrt(2.)); }
         if(xyz==2) {  imag=true;retval = (qpp+qpm)/sqrt(2.); }// changed MR 25.5.2010
      }
   }
   else if(xyz==3)
   {
      nstr[1]=83;
      if(Q!=0)
      {
         NSTRB(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MSq(0,K,abs(Q),_n,_l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTRB(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MSq(0,K,-abs(Q),_n,_l); rmzeros(qmp); mm_gout(qmp,filename); }
         if(Q<0) {  if(Q%2==0) qmp -= qpp; else qmp += qpp;  imag=true;retval = qmp/(-sqrt(2.));}
         else    {  if(Q%2==0) qmp += qpp; else qmp -= qpp;  retval = qmp/sqrt(2.);}// changed by MR 28.5.2010
      }
      else
      {
         NSTRB(K,0); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         retval = mm_gin(filename); if(retval.isempty()) { retval = balcar_MSq(0,K,0,_n,_l); rmzeros(retval); mm_gout(retval,filename); }
      }
   }
   else if(xyz==-1||xyz==-2)
   {
      nstr[1]=76;
      if(Q!=0)
      {
         NSTRB(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MLq(1,K,abs(Q),_n,_l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTRB(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MLq(1,K,-abs(Q),_n,_l); rmzeros(qmp); mm_gout(qmp,filename); }
         NSTRB(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MLq(-1,K,abs(Q),_n,_l); rmzeros(qpm); mm_gout(qpm,filename); }
         MSTRB(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmm = mm_gin(filename); if(qmm.isempty()) { qmm = balcar_MLq(-1,K,-abs(Q),_n,_l); rmzeros(qmm); mm_gout(qmm,filename); }
	 // sum to coeff of Zlm (neglecting a 1/sqrt(2) factor) 
         if(Q<0) { if(Q%2==0) { qmp -= qpp; qmm -= qpm; } else { qmp += qpp; qmm += qpm; } }
         else    { if(Q%2==0) { qmp += qpp; qmm += qpm; } else { qmp -= qpp; qmm -= qpm; } }
         // add spherical components and multiply by addition factor 1/sqrt(2)which was neglected in the line above
         if(xyz==-1) { if(Q<0) {imag=true;retval = (qmm-qmp)/(-2.);}    else retval = (qmm-qmp)/2.; }
         if(xyz==-2) { if(Q<0) retval = (qmm+qmp)/2.; else {imag=true;retval = (qmm+qmp)/2.;} }// changed MR 25.5.2010 // Q<0 signs changed MR 28.5.2010
      }
      else
      {
         NSTRB(K,0); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MLq(1,K,0,_n,_l); rmzeros(qpp); mm_gout(qpp,filename); }
         NSTRB(K,0); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MLq(-1,K,0,_n,_l); rmzeros(qpm); mm_gout(qpm,filename); }
         if(xyz==-1) {  retval = (qpp-qpm)/(-sqrt(2.)); }
         if(xyz==-2) {  imag=true;retval = (qpp+qpm)/sqrt(2.); }// changed MR 25.5.2010
      }
   }
   else if(xyz==-3)
   {
      nstr[1]=76;
      if(Q!=0)
      {
         NSTRB(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MLq(0,K,abs(Q),_n,_l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTRB(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MLq(0,K,-abs(Q),_n,_l); rmzeros(qmp); mm_gout(qmp,filename); }
         if(Q<0) {  if(Q%2==0) qmp -= qpp; else qmp += qpp; imag=true;retval = qmp/(-sqrt(2.));}
         else    {  if(Q%2==0) qmp += qpp; else qmp -= qpp; retval = qmp/sqrt(2.);}// changed by MR 28.5.2010
      }
      else
      {
         NSTRB(K,0); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         retval = mm_gin(filename); if(retval.isempty()) { retval = balcar_MLq(0,K,0,_n,_l); rmzeros(retval); mm_gout(retval,filename); }
      }
   }
   return retval;
}

sMat<double> icmfmat::op_generate(int i)
{if(_num_op<=i) {_num_op = i+1; iflag.resize(_num_op,0); // extend operator storage if more operators are required
                 sMat<double> t; J.resize(_num_op,t);T.resize(_num_op+2,NULL);
                                 }
 if(J[i].isempty()) // only do something if Operator matrix is empty
 {if(_xyz!=0)
 {// Indices fpr balcar Mq 4-8 are k=2 quadrupoles; 9-15:k=3; 16-24:k=4; 25-35:k=5; 36-48:k=6
   int k[] = {0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
   int q[] = {0,-1,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
   if(i>48){std::cerr << "Error icmfmat::op_generate - generating Operator balcar_Mq for " << i << " > 48 not possible\n";exit(EXIT_FAILURE);}
   bool imag; J[i]=balcar_Mq(_xyz,k[i], q[i], imag );iflag[i]=imag;}
 else
 {if(i>50){std::cerr << "Error icmfmat::op_generate - generating Operator " << i << " > 50 not possible\n";exit(EXIT_FAILURE);}
  char nstr[6];  char basename[255]; strcpy(basename,"results/mms/");
      if(_save_matrices) {
      #ifndef _WINDOWS
      struct stat status; stat("results/mms",&status); if(!S_ISDIR(status.st_mode))
         if(mkdir("results/mms",0777)!=0) std::cerr << "icmfmat::op_generate(): Can't create mms dir, " << strerror(errno) << "\n";
      #else
      DWORD drAttr = GetFileAttributes("results\\mms"); if(drAttr==0xffffffff || !(drAttr&FILE_ATTRIBUTE_DIRECTORY)) 
         if (!CreateDirectory("results\\mms", NULL)) std::cerr << "icmfmat::op_generate(): Cannot create mms directory\n";
      #endif
      nstr[0] = (_l==F?102:100); if(_n<10) { nstr[1] = _n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = _n+38; nstr[3] = 0; }
      strcat(basename,nstr); strcat(basename,"_"); 
      } else { strcpy(basename,"nodir/"); }
 if(i<6){
    // Determines the filename strings for where the moment operator matrices are stored if previously calculated
    char Lfilestr[255], Sfilestr[255]; 
   if(i==0||i==1||i==3||i==4)
  {
      iflag[1]=1; iflag[4]=1;
   nstr[0] = 76;   // 76 is ASCII for "L", 85=="U", 100=="d" and 102=="f"
   nstr[1] = 49; // 49=="1"
  
   // Calculates the L and S operator matrix for each direction
   sMat<double> Sp1, Sm1, Lp1, Lm1; 
   nstr[2]=120; nstr[3]=0; strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");   
   nstr[0]=83;             strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");
   J[3] = mm_gin(Lfilestr); J[0] = mm_gin(Sfilestr); 
   nstr[2]=121; nstr[3]=0; strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");   
   nstr[0]=76;             strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");
   J[4] = mm_gin(Lfilestr); J[1] = mm_gin(Sfilestr); 
   if(J[0].isempty() || J[1].isempty() || J[3].isempty() || J[4].isempty())
    { 
      racah_mumat(_n,1,Lp1,Sp1,_l); rmzeros(Sp1); rmzeros(Lp1);
      racah_mumat(_n,-1,Lm1,Sm1,_l); rmzeros(Sm1); rmzeros(Lm1);
      J[0] = (Sm1-Sp1)/sqrt(2); J[1] = (Sm1+Sp1)/sqrt(2); Sm1.clear(); Sp1.clear();   // Sx and Sy
      J[3] = (Lm1-Lp1)/sqrt(2); J[4] = (Lm1+Lp1)/sqrt(2); Lm1.clear(); Lp1.clear();   // Lx ans Ly
      mm_gout(J[1],Sfilestr); mm_gout(J[4],Lfilestr);
      nstr[2]=120; nstr[3]=0; strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");   
      nstr[0]=83;             strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");
      mm_gout(J[0],Sfilestr); mm_gout(J[3],Lfilestr);
    }
   
  }else
  {
   nstr[0]=76; nstr[2]=122; nstr[3]=0; strcpy(Lfilestr,basename); strcat(Lfilestr,nstr); strcat(Lfilestr,".mm");
   nstr[0]=83;                         strcpy(Sfilestr,basename); strcat(Sfilestr,nstr); strcat(Sfilestr,".mm");
   J[2] = mm_gin(Sfilestr); J[5] = mm_gin(Lfilestr);                               // Sz and Lz
   if(J[2].isempty() || J[5].isempty()) { 
      racah_mumat(_n,0,J[5],J[2],_l); rmzeros(J[2]); rmzeros(J[5]); mm_gout(J[2],Sfilestr); mm_gout(J[5],Lfilestr); 
                                        }
   
  }
/* // Checks the moment operator matrices against those given by Chan and Lam.
// int ii,jj=0; sMat<double> mu; double g_s = 2.0023193043622; // electronic g-factor
// chanlam_mumat(n,1,mu,l); for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    if(fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[0](ii,jj))>10*DBL_EPSILON) { std::cerr << "icmfmat: Magnetic moment operator x does not agree.\n"; break; }
//    if(ii==mu.nr() && jj==mu.nc()) std::cerr << "icmfmat: Magnetic moment operator x agrees.\n";
// chanlam_mumat(n,2,mu,l); for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    if(fabs(-mu(ii,jj)-J[4](ii,jj)-g_s*J[1](ii,jj))>10*DBL_EPSILON) { std::cerr << "icmfmat: Magnetic moment operator y does not agree.\n"; break; }
//    if(ii==mu.nr() && jj==mu.nc()) std::cerr << "icmfmat: Magnetic moment operator y agrees.\n";
// chanlam_mumat(n,3,mu,l); for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    if(fabs(-mu(ii,jj)-J[5](ii,jj)-g_s*J[2](ii,jj))>10*DBL_EPSILON) { std::cerr << "icmfmat: Magnetic moment operator z does not agree.\n"; break; }
//    if(ii==mu.nr() && jj==mu.nc()) std::cerr << "icmfmat: Magnetic moment operator z agrees.\n";
   double sumcheck;
   chanlam_mumat(n,1,mu,l); sumcheck = 0.; for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    std::cout << -mu(ii,jj) << "\t" << J[3](ii,jj)+g_s*J[0](ii,jj)  << "\t" << fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[0](ii,jj)) << "\n";
      sumcheck += fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[0](ii,jj)); std::cout << "Moment Matrix Check: sum(-mu_x(ChanLam) - (Lx+gSx)) = " << sumcheck << "\n";
   chanlam_mumat(n,2,mu,l); sumcheck = 0.; for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
//    std::cout << -mu(ii,jj) << "\t" << J[4](ii,jj)+g_s*J[1](ii,jj)  << "\t" << fabs(-mu(ii,jj)-J[3](ii,jj)-g_s*J[2](ii,jj)) << "\n";
      sumcheck += fabs(-mu(ii,jj)-J[4](ii,jj)-g_s*J[1](ii,jj)); std::cout << "Moment Matrix Check: sum(-mu_y(ChanLam) - (Ly+gSy)) = " << sumcheck << "\n";
   chanlam_mumat(n,3,mu,l); sumcheck = 0.; for(ii=0; ii<mu.nr(); ii++) for(jj=0; jj<mu.nc(); jj++) 
      sumcheck += fabs(-mu(ii,jj)-J[5](ii,jj)-g_s*J[2](ii,jj)); std::cout << "Moment Matrix Check: sum(-mu_z(ChanLam) - (Lz+gSz)) = " << sumcheck << "\n"; 
*/
        } // 50>=i>=6
    else{
    char filename[255];
            // for i>=6 calculates operator J[i] either loading Umq Upq from file or calculating the matrices
      nstr[0] = 85;   // 85 is ASCII for "U", 100=="d" and 102=="f"
#define NSTR(K,Q) nstr[1] = K+48; nstr[2] = Q+48; nstr[3] = 0
#define MSTR(K,Q) nstr[1] = K+48; nstr[2] = 109;  nstr[3] = Q+48; nstr[4] = 0
      // Indices 6-10 are k=2 quadrupoles; 11-17:k=3; 18-26:k=4; 27-37:k=5; 38-50:k=6
      int k[] = {1,1,1,1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
      int q[] = {0,0,0,0,0,0,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};
    //int im[]= {0,0,1,1,0,0, 1, 1,0,0,0, 1, 1, 1,0,0,0,0, 1, 1, 1, 1,0,0,0,0,0, 1, 1, 1, 1, 1,0,0,0,0,0,0, 1, 1, 1, 1, 1, 1,0,0,0,0,0,0,0};
      sMat<double> Upq,Umq; double redmat; int n = _n; //if(n>(2*_l+1)) n = 4*_l+2-n; 
            if(q[i]<0) iflag[i]=1;
            if(k[i]%2==1){J[i].zero(J[0].nr(),J[0].nc());} // continue;   // Using the  reduced matrix element with at (l k l; 0 0 0) 3-j symbol, odd k gives zero...
            else if(k[i]>_l*2){J[i].zero(J[0].nr(),J[0].nc());} //continue;
            else {
            redmat = pow(-1.,(double)abs(_l)) * (2*_l+1) * threej(2*_l,2*k[i],2*_l,0,0,0);// * wy2stev(i);
            NSTR(k[i],abs(q[i])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Upq = mm_gin(filename); if(Upq.isempty()) { Upq = racah_ukq(n,k[i],abs(q[i]),_l); rmzeros(Upq); mm_gout(Upq,filename); }
            if(q[i]==0) { J[i]= Upq * (redmat); }
            else {
            MSTR(k[i],abs(q[i])); strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
            Umq = mm_gin(filename); if(Umq.isempty()) { Umq = racah_ukq(n,k[i],-abs(q[i]),_l); rmzeros(Umq); mm_gout(Umq,filename); }
            if(q[i]<0) { 
               if((q[i]%2)==0)  J[i]= (Umq - Upq) * (redmat); else   J[i]= (Umq + Upq) * (redmat); }
            else {
               if((q[i]%2)==0)  J[i]= (Umq + Upq) * (redmat); else  J[i]= (Umq - Upq) * (redmat); } 
                 }
               }
      }
 
  }
 }
return J[i];
}

void icmfmat::op_free(int i)
{
if(_save_matrices&&i>5) {J[i].clear();}
}


// --------------------------------------------------------------------------------------------------------------- //
// Calculates the mean field matrix sum_i (H_i*J_i)
// --------------------------------------------------------------------------------------------------------------- //
void icmfmat::Jmat(sMat<double>&Jmat, sMat<double>&iJmat, std::vector<double>&gjmbH)
{  int i; Jmat.zero(J[0].nr(),J[0].nc()); iJmat.zero(J[0].nr(),J[0].nc()); 
   if(_num_op<(int)gjmbH.size()) {_num_op = (int)gjmbH.size(); iflag.resize(_num_op,0); 
                                   sMat<double> t; J.resize(_num_op,t);T.resize(_num_op+2,NULL);
                                 }
   for(i=0; i<((int)gjmbH.size()>6?6:_num_op); i++)
      if(fabs(gjmbH[i])>DBL_EPSILON*100) { if(iflag[i]==1) iJmat += J[i]*gjmbH[i]; else Jmat += J[i]*gjmbH[i]; }
   // Higher order than dipole operators needed
   if((int)gjmbH.size()>6)
   {
     for(i=6; i<(int)gjmbH.size(); i++)if(fabs(gjmbH[i])>DBL_EPSILON)
      { 
         op_generate(i);            
         if(iflag[i]==1) iJmat += J[i]*gjmbH[i]; else Jmat += J[i]*gjmbH[i];
         op_free(i);
      }
   }
}
// --------------------------------------------------------------------------------------------------------------- //
// Calculates the expectation values <V|J|V>exp(-beta*T) given a set of eigenstates
// --------------------------------------------------------------------------------------------------------------- //

std::vector<double>  icmfmat::expJ(iceig &VE, double T, std::vector< std::vector<double> > &matel, int num_op)
{ std::vector<double> E,eb, ex((num_op+2),0.), me; matel.clear();
   int iJ, ind_j, Esz, Hsz=VE.Hsz();op_generate(0); 
   if(Hsz!=J[0].nr()) { std::cerr << "icmfmat::expJ() - Hamiltonian matrix size not same as mean field operator!\n"; return E; }
   sMat<double> zeroes; zeroes.zero(J[0].nr(),J[0].nc());
   //complexdouble zalpha; zalpha.r=1; zalpha.i=0; 
   //complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   // Checks that the eigenvactors are orthonormal
/* char transa='C', transb='N'; double summm=0.;
   if(VE.iscomplex())
   {
      complexdouble *zmm = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble)); 
      complexdouble *vet = (complexdouble*)malloc(Hsz*Hsz*sizeof(complexdouble)); memcpy(vet,VE.zV(0),Hsz*Hsz*sizeof(complexdouble));
      F77NAME(zgemm)(&transa, &transb, &Hsz, &Hsz, &Hsz, &zalpha, vet, &Hsz, VE.zV(0), &Hsz, &zbeta, zmm, &Hsz);
      for(int ii=0; ii<Hsz; ii++) { zmm[ii*Hsz+ii].r-=1.; summm += F77NAME(dzasum)(&Hsz, &zmm[ii*Hsz], &incx); if(VE.E(ii+1)==0) break; }
      std::cout << "#ic1ion: Sum(V^TV-I) = " << summm << "\n";
      free(zmm); free(vet);
   }
   else
   {
      double *dmm = (double*)malloc(Hsz*Hsz*sizeof(double)); 
      double *vet = (double*)malloc(Hsz*Hsz*sizeof(double)); memcpy(vet,VE.V(0),Hsz*Hsz*sizeof(double));
      F77NAME(dgemm)(&transa, &transb, &Hsz, &Hsz, &Hsz, &alpha, vet, &Hsz, VE.V(0), &Hsz, &beta, dmm, &Hsz);
      for(int ii=0; ii<Hsz; ii++) { dmm[ii*Hsz+ii]-=1.; summm += F77NAME(dasum)(&Hsz, &dmm[ii*Hsz], &incx); if(VE.E(ii+1)==0) break; }
      std::cout << "#ic1ion: Sum(V^TV-I) = " << summm << "\n";
      free(dmm); free(vet);
   }
*/

   for(int ii=0; ii<Hsz; ii++) for(int jj=0; jj<Hsz; jj++) 
      if(fabs(VE.zV(ii,jj).r*VE.zV(ii,jj).r+VE.zV(ii,jj).i*VE.zV(ii,jj).i)<DBL_EPSILON*100000) 
      {
         VE.zV(ii,jj) = 0;
      }  
   
   if(_xyz!=0) 
   { 
      std::cout << "#Calculating the expectation of the moment density operator ";
      switch(_xyz)
      { case 1: std::cout << "Sx\n";break;
        case 2: std::cout << "Sy\n";break;
        case 3: std::cout << "Sz\n";break;
        case -1: std::cout << "Lx\n";break;
        case -2: std::cout << "Ly\n";break;
        case -3: std::cout << "Lz\n";break;
        default: break;
      }
      
   }
// Sets energy levels relative to lowest level, and determines the maximum energy level needed.
   for(Esz=0; Esz<J[0].nr(); Esz++) { E.push_back(VE.E(Esz)-VE.E(0)); if(exp(-E[Esz]/(KB*T))<DBL_EPSILON || VE.E(Esz+1)==0 || VE.E(Esz+1)==-DBL_MAX) break; }

   if (T<0){printf ("Temperature T=%g<0: please choose probability distribution for the -T=%i lowest energy states by hand\n",T,(int)(-T));
                         printf ("Number   Excitation Energy\n");
     for (ind_j=0;ind_j<(int)(-T);++ind_j) printf ("%i    %4.4g meV\n",ind_j+1,E[ind_j]);
     } // MR 10.9.2010

   // calculate Matrix Elements 
   for(iJ=0; iJ<num_op; iJ++)
   { op_generate(iJ);            
      
      me.assign(Esz,0.);
     if(!J[iJ].isempty()) // might be empty because redmat is zero ...
      {if(!VE.iscomplex() && _xyz==0&&iflag[iJ]==0)
       {      
            for(ind_j=0; ind_j<Esz; ind_j++)
            { // Calculates the matrix elements <Vi|J|Vi>
          me[ind_j]=J[iJ].MultvxMv(VE.V(ind_j));
            }
       }  
      else
      {   for(ind_j=0; ind_j<Esz; ind_j++)
         {  // Calculates the matrix elements <Vi|J|Vi>
           me[ind_j]=J[iJ].MultvxMv(VE.zV(ind_j),iflag[iJ]);
                   }
         
       }
      }
    matel.push_back(me); 

    op_free(iJ);
   }

// For first run calculate also the partition function and internal energy
// Rest of the runs only calculate the new matrix elements

  double U=0;
  double Z=0;eb.assign(Esz,0.);
  for(iJ=0; iJ<num_op; iJ++)
   {ex[iJ]=0;
    for(ind_j=0; ind_j<Esz; ind_j++)
    {
     if(iJ==0)
        {if (T<0)
         { Esz=(int)(-T); char instr[MAXNOFCHARINLINE];
            printf("eigenstate %i: %4.4g meV  - please enter probability w(%i):",ind_j+1,E[ind_j],ind_j+1);
            if(fgets(instr, MAXNOFCHARINLINE, stdin)==NULL) { printf("Error in input. Exiting\n"); exit(-1); }
            eb[ind_j]=strtod(instr,NULL);
         }
         else
         { eb[ind_j] = exp(-E[ind_j]/(KB*T));} 
         Z+=eb[ind_j]; 
         U+=(E[ind_j]+VE.E(0))*eb[ind_j];
        }  
        ex[iJ]+=matel[iJ][ind_j]*eb[ind_j];
    }       
    ex[iJ]/=Z; 
    if(fabs(ex[iJ])<DBL_EPSILON) ex[iJ]=0.; 
            
 } // iJ
 ex[iJ] = log(Z)-VE.E(0)/(KB*T); // set lnZ
 ex[iJ+1] = U/Z;// set U
    return ex;
}



// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix M_ab=<i|Ja|j><j|Jb|i>{exp(-beta_i*T)-exp(-beta_j*T)} for some state i,j
// --------------------------------------------------------------------------------------------------------------- //
int icmfmat::u1(complexdouble*u,int sz, double T,// * sqrt{exp(-beta_i*T)-exp(-beta_j*T)}
        int tn, float &delta,complexdouble *V,complexdouble * ev,int Hsz,int & n,int &nd)
{// check if printout should be done and make tn positive
   int pr=0; if (tn<0) { pr=1; tn*=-1; }
   double ninit=u[0].r;
   double pinit=u[0].i;
   int i,j=0,k=0; for(i=0; i<Hsz; ++i) { for(j=i; j<Hsz; ++j) { ++k; if(k==tn) break; } if(k==tn) break; }
   double maxE=delta;n=i;nd=j;
   if((delta=(ev[j].r-ev[i].r))<=maxE)
   {  double *en = new double[Hsz]; for(k=0; k<Hsz; k++) en[k] = ev[k].r;
      // catch eigenvectors VE from est
      iceig VE(Hsz,en,V,1);
      // Calculates the transition matrix elements:
      //    u1 = <i|Ja|j> * sqrt[(exp(-Ei/kT)-exp(-Ej/kT)) / Z ]   if delta > small
      //    u1 = <i|Ja-<Ja>|j> * sqrt[(exp(-Ei/kT)) / kTZ ]             if delta < small (quasielastic scattering)
      //    See file icpars.cpp, function mfmat::Mab() to see the actual code to calculate this.
      double  Z=0., therm; complexdouble zme; zme.r=0; zme.i=0.;
      std::vector<double> mij(sz,0.);
      int iJ;
     // Calculates the matrix elements: <i|Ja|j> and <j|Ja|i> for each of the  Ja's
     for(iJ=0; iJ<sz; iJ++)
     {u[iJ]=zme;
      op_generate(iJ); 
      if(!VE.iscomplex() && iflag[iJ]==0)
      { u[iJ].r=J[iJ].MultuxMv(VE.V(i),VE.V(j));
      } 
      else
      { u[iJ]=J[iJ].MultuxMv(VE.zV(i),VE.zV(j),iflag[iJ]);
       }
      op_free(iJ);
     }

   if(i==j&&T>0) {//subtract thermal expectation value from zij=zii
            std::vector< std::vector<double> > matel;
            std::vector<double> vJ = expJ(VE,T,matel,sz);
           for(iJ=0; iJ<sz; iJ++)u[iJ].r-=vJ[iJ];
            }
   if (T<0){T=-T;}

   delta = VE.E(j)-VE.E(i);
   if(delta<-0.000001)
   {
      std::cerr << "ERROR module ic1ion - du1calc: energy gain delta gets negative\n"; 
      exit(EXIT_FAILURE);
   }
   if(j==i)delta=-SMALL; // if transition within the same level: take negative delta !!- this is needed in routine intcalc

   // Calculates the partition function
   for(iJ=0; iJ<Hsz; iJ++) { therm = exp(-(VE.E(iJ)-VE.E(0))/(KB*T)); Z += therm; if(therm<DBL_EPSILON) break; }

   // do some printout if wishes and set correct occupation factor
   if (delta>SMALL)
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T)) - exp(-(VE.E(j)-VE.E(0))/(KB*T));
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         for(iJ=0;iJ<sz;++iJ)
         {switch(_xyz)
           {case  0: printf(" |<%i|I%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case  1: printf(" |<%i|Msx%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case  2: printf(" |<%i|Msy%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case  3: printf(" |<%i|Msz%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case -1: printf(" |<%i|Mlx%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case -2: printf(" |<%i|Mly%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case -3: printf(" |<%i|Mlz%i|%i>|^2=%6.3f\t ",i+1,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
         }}
         printf(" n%i-n%i=%6.3f\n",i,j,therm / Z);
      }
   }
   else
   {
      therm = exp(-(VE.E(i)-VE.E(0))/(KB*T))/(KB*T);    // quasielastic scattering has not wi-wj but wj*epsilon/kT
      if(pr==1)
      {
         printf("delta(%i->%i)=%6.3fmeV\n",i+1,j+1,delta);
         for(iJ=0;iJ<sz;++iJ)
         {switch(_xyz)
          { case  0: printf(" |<%i|I%i-<I%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case  1: printf(" |<%i|Msx%i-<Msx%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case  2: printf(" |<%i|Msy%i-<Msy%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case  3: printf(" |<%i|Msz%i-<Msz%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case -1: printf(" |<%i|Mlx%i-<Mlx%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case -2: printf(" |<%i|Mly%i-<Mly%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
           case -3: printf(" |<%i|Mlz%i-<Mlz%i>|%i>|^2=%6.3f\t ",i+1,iJ,iJ,j+1,u[iJ].r*u[iJ].r+u[iJ].i*u[iJ].i);break;
         }}printf(" n%i=%6.3f\n",i,(KB*T)*therm/Z);
      }
   }

   // multiply matrix Mab by occupation factor
   for(iJ=0; iJ<sz; iJ++)
      { u[iJ].r *= sqrt(therm/Z); u[iJ].i *= sqrt(therm/Z); }

     }
   // determine number of thermally reachable states
   if (ninit>Hsz)ninit=Hsz;
   double zsum=0,zi,x;
   int noft=0; 
   if(T>0)for(i=0; (i<ninit)&((((x=(ev[i].r-ev[0].r)/(KB*fabs(T)))<200)? zi=exp(-x):zi=0)>=(pinit*zsum)); ++i)
   {//fprintf(stderr,"i=%i zi=%g zsum=%g noft=%i Hsz=%i Ei-E1=%g\n",i,zi,zsum,noft,Hsz,ev[0][i].r-ev[0].r);
      noft += Hsz-i; 
      zsum += zi;
   }
   if(T<0)for(int i=0;i<ninit;++i){noft+=Hsz-i;}
   return noft;
   
}
     


