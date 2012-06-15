/* spectre.cpp
 *
 * Routines to calculate the IC Hamiltonian using the method in the Spectra/XTAL programs of Hannah Crosswhite.
 * References: https://chmwls.chm.anl.gov/downloads/spectra/about/matrix.html
 *             Carnall, Crosswhite, Crosswhite, Conway, J. Chem. Phys., v64, p3582, 1976
 *
 * Functions:
 *   iceig spectre_eig(&pars, elim)                                // Diagonalises a truncated Hamiltonian
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008-2011 Duc Le - duc.le@ucl.ac.uk
 */

#include "ic1ion.hpp"

// --------------------------------------------------------------------------------------------------------------- //
// Calculates a rotation matrix to transform the full IC Hamiltonian Hic into a basis in which the free-ion
//   Hamiltonian, Hc+Hso, is diagonal. 
//
// First the basis states are re-ordered with states having the same J-values following each other consecutively.
//   Then the Hc+Hso interactions are calculated in this basis, yielding a matrix of blocks along the diagonal
//   (Hc and Hso are diagonal in J,mJ). Each J-submatrix is diagonalised individually to save time, and the full
//   eigenvector matrix Vrot (also block-diagonal) is reconstructed from the individual J-submatrices.
//
// Vrot is returned by the function after rearranging such that it's rows are in the original (non-J-ordered)
//   basis states.
// --------------------------------------------------------------------------------------------------------------- //
void spectre_vrot(icpars &pars, double &elim, bool isreal=true, double *dVrot=0, complexdouble *zVrot=0)
{
   if(pars.e_units.find("meV")!=std::string::npos)    elim *= MEV2CM;
   else if(pars.e_units.find("K")!=std::string::npos) elim /= CM2K;

   // Calculates the Coulomb and Spin orbit Hamiltonian matrix in the |alpha,LSJ> basis
   sMat<double> Hcso = ic_Hcso(pars);

   // Finds indices to allow reordering the basis states so states with same J are grouped together
   unsigned int i, j, num_states=Hcso.nr(), ns=getdim(pars.n,pars.l), imax=num_states;
   int ii=0, jj, imj=0, J2;
   fconf conf(pars.n,0,pars.l);
   std::vector< std::vector<int> > JsortM, JsortMmj; std::vector<int> Jind, JsortV(ns,0);
   for(i=0; i<26; i++) { JsortM.push_back(Jind); JsortMmj.push_back(Jind); }
   for(i=0; i<conf.states.size(); i++) 
   { 
      J2=conf.states[i].J2; JsortM[J2].push_back(i); 
      for(j=0; j<(unsigned int)J2+1; j++) JsortMmj[J2].push_back(imj++); 
   }
   for(i=JsortM.size()-1; i<JsortM.size(); i--) 
      if(!JsortM[i].empty()) 
      {
         for(j=0; j<JsortM[i].size(); j++) Jind.push_back(i);
         for(j=0; j<JsortMmj[i].size(); j++) JsortV[ii++]=JsortMmj[i][j];
      }
   i=0; while(i<JsortM.size()) { if(JsortM[i].empty()) JsortM.erase(JsortM.begin()+i); else i++; }

   // Split Hc+Hso into J-submatrices
   double **submats = new double*[JsortM.size()];
   for(i=0; i<JsortM.size(); i++) {
      submats[i] = new double[JsortM[i].size()*JsortM[i].size()];
      memset(submats[i],0,JsortM[i].size()*JsortM[i].size()*sizeof(double)); }
   std::vector< std::vector<int> > nz = Hcso.find();
   std::vector<int> i1(num_states,0),i2(num_states,0),i3(num_states,0);
   for(i=0; i<JsortM.size(); i++) for(j=0; j<JsortM[i].size(); j++) { i1[JsortM[i][j]] = i; i2[JsortM[i][j]]=j; }
   for(unsigned int kk=0; kk<nz.size(); kk++) 
   {
      ii = nz[kk][0]; jj = nz[kk][1];
      if(i1[ii]==i1[jj]) submats[i1[ii]][i2[jj]*JsortM[i1[ii]].size()+i2[ii]] = Hcso(ii,jj);
   }

   // Diagonalises the submatrices, and reconstruct the full eigenvector matrix from the eigenvector 
   //   matrices of the submatrices, and expand to the |alpha,LSJmJ> basis
   double *Vunrot = new double[num_states*num_states], *Eunrot = new double[num_states]; 
   memset(Vunrot,0,num_states*num_states*sizeof(double));
   int idx = 0;
   iceig VEs;
   for(i=JsortM.size()-1; i<JsortM.size(); i--)
   {
      int isz = JsortM[i].size();
      VEs.calc(isz,submats[i]); delete[]submats[i];
      for(ii=0; ii<isz; ii++) { 
         for(int jj=0; jj<isz; jj++) Vunrot[(jj+idx)*num_states+ii+idx] = VEs.V(ii,jj); Eunrot[ii+idx] = VEs.E(ii); } 
      idx+=isz; 
   }
   delete[]submats;

   // Finds the indices of the sorted Energies of Hc+Hso, to determine the state cb at which to
   //    truncate the eigenvector matrix used to rotate Hic from the maximum energy limit elim
   std::vector<int> Esort(num_states,0); for(i=0; i<num_states; i++) Esort[i] = i;
   i=1;j=2;
   while(i<num_states)
   {
      if(Eunrot[Esort[i-1]]<=Eunrot[Esort[i]]) { i=j; j++; }
      else { ii=Esort[i-1]; Esort[i-1]=Esort[i]; Esort[i]=ii; i--; if(i==0) i=1; }
   }
   if(elim>=0) for(imax=0; imax<num_states; imax++) if((Eunrot[Esort[imax]]-Eunrot[Esort[0]])>elim) break; imax--;
   delete[]Eunrot;

   // Construct the rotation matrix from the full eigenvector matrix, but putting the columns in Energy-order
   unsigned int iy=0, jy=0, id, J2i, J2j;
   for(i=0; i<num_states; i++) 
   { 
      jy=0; J2i = Jind[i];
      for(j=0; j<num_states; j++) 
      { 
         J2j = Jind[Esort[j]];
         if(fabs(Vunrot[Esort[j]*num_states+i])>SMALL) {
            if(isreal) for(id=0; id<J2i+1; id++) dVrot[JsortV[iy+id]+(jy+id)*ns]   = Vunrot[Esort[j]*num_states+i];
            else       for(id=0; id<J2i+1; id++) zVrot[JsortV[iy+id]+(jy+id)*ns].r = Vunrot[Esort[j]*num_states+i]; }
         jy+=J2j+1;
         if(i==0) if(j==imax) elim=(double)jy;
      } 
      iy+=J2i+1;
   }
   delete[]Vunrot;
}

void spectre_vrot(icpars &pars, double *Vrot, double &elim) { spectre_vrot(pars,elim,true,Vrot); }
void spectre_vrot(icpars &pars, complexdouble *Vrot, double &elim) { spectre_vrot(pars,elim,false,0,Vrot); }

// --------------------------------------------------------------------------------------------------------------- //
// Calculates eigenvalues/vectors of the IC Hamiltonian by a method similar to that used in the XTALS84 programs
//   of Hannah Crosswhite. 
//
// A rotation matrix is first calculated using spectre_vrot(). This matrix is truncated at some energy [elim] 
//   corresponding to a level [cb] so that only cb columns of the matrix are used to rotate the full IC Hamiltonian 
//   by: V(:,1:cb)'*Hic*V(:,1:cb), yielding a cb*cb matrix which is diagonalised. The eigenvalues are the energies, 
//   and the eigenvectors in the original basis is given by (Vr*Vrot')' if Vr are eigenvectors in the rotated basis.
// --------------------------------------------------------------------------------------------------------------- //
iceig spectre_eig(sMat<double> Hic, double *Vrot, int cb)
{
   int i, j, ns = (int)Hic.nr();
   // Rotates and truncates Hic by using the eigenvectors of Hc+Hso: Hic_r = V(:,1:cb)'*Hic*V(:,1:cb)
   double *Jmat,*Hrot,*mt; mt = new double[ns*cb]; Hrot = new double[cb*cb]; Jmat = Hic.f_array();
   char notranspose='N',transpose='C',uplo='U',side='L'; double alpha=1, beta=0;
   F77NAME(dsymm)(&side,&uplo,&ns,&cb,&alpha,Jmat,&ns,Vrot,&ns,&beta,mt,&ns);
   F77NAME(dgemm)(&transpose,&notranspose,&cb,&cb,&ns,&alpha,Vrot,&ns,mt,&ns,&beta,Hrot,&cb); free(Jmat);

   // Diagonalises the truncated and rotated matrix
   iceig VEr; VEr.calc(cb,Hrot); delete[]Hrot;

   // Transform the eigenvectors of the rotated basis back into the original |alpha,LSJmJ> basis
   memset(mt,0.,ns*cb*sizeof(double)); double *Vr = new double[cb*ns*sizeof(double)]; 
   for(i=0; i<ns; i++) for(j=0; j<cb; j++) Vr[i*cb+j] = VEr.V(i,j);
   F77NAME(dgemm)(&notranspose,&transpose,&cb,&ns,&cb,&alpha,Vr,&cb,Vrot,&ns,&beta,mt,&cb);  // (Vr*Vrot')'
   delete[]Vr;
   iceig VE(ns,true); double *En = VE.E(); for(i=0; i<cb; i++) En[i] = VEr.E(i); for(i=cb; i<ns; i++) En[i] = -DBL_MAX;
   double *V = VE.V(0); for(i=0; i<ns; i++) for(j=0; j<cb; j++) V[j*ns+i] = mt[i*cb+j]; delete[]mt;

   return VE;
}
iceig spectre_eig(sMat<double> Hic, sMat<double> iHic, complexdouble *Vrot, int cb)
{
   int i, j, ns = (int)Hic.nr(); 
   iceig VE(ns,false);
   
   // Rotates and truncates Hic by using the eigenvectors of Hc+Hso: Hic_r = V(:,1:cb)'*Hic*V(:,1:cb)
   complexdouble *zJmat,*Hrot,*zmt; zmt = new complexdouble[ns*cb]; Hrot = new complexdouble[cb*cb];
   char notranspose='N',transpose='C',uplo='U',side='L'; 
   complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   zJmat = zmat2f(Hic,iHic);
   F77NAME(zhemm)(&side,&uplo,&ns,&cb,&zalpha,zJmat,&ns,Vrot,&ns,&zbeta,zmt,&ns);
   F77NAME(zgemm)(&transpose,&notranspose,&cb,&cb,&ns,&zalpha,Vrot,&ns,zmt,&ns,&zbeta,Hrot,&cb); free(zJmat);

   // Diagonalises the truncated and rotated matrix
   iceig VEr; VEr.calc(cb,Hrot); delete[]Hrot;

   // Transform the eigenvectors of the rotated basis back into the original |alpha,LSJmJ> basis
   memset(zmt,0.,ns*cb*sizeof(complexdouble));
   complexdouble *Vr = new complexdouble[cb*ns*sizeof(complexdouble)]; 
   for(i=0; i<ns; i++) for(j=0; j<cb; j++) Vr[i*cb+j] = VEr.zV(i,j);
   F77NAME(zgemm)(&notranspose,&transpose,&cb,&ns,&cb,&zalpha,Vr,&cb,Vrot,&ns,&zbeta,zmt,&cb);  // (Vr*Vrot')'
   delete[]Vr;
   double *En = VE.E(); complexdouble *zV = VE.zV(0); 
   for(i=0; i<cb; i++) En[i] = VEr.E(i); for(i=cb; i<ns; i++) En[i] = -DBL_MAX;
   for(i=0; i<ns; i++) for(j=0; j<cb; j++) { zV[j*ns+i].r = zmt[i*cb+j].r; zV[j*ns+i].i = -zmt[i*cb+j].i; }
   delete[]zmt; 

   return VE;
}
