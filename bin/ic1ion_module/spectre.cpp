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
// Calculates the eigenvalues/vectors of the IC Hamiltonian by a method similar to that used in the XTALS84 programs
//   of Hannah Crosswhite. First the basis states are re-ordered with states having the same J-values following
//   each other consecutively. Then the Hc+Hso interactions are calculated in this basis, yielding a matrix of
//   blocks along the diagonal (Hc and Hso are diagonal in J,mJ). Each J-submatrix is diagonalised individually
//   to save time, and the full eigenvector matrix Vrot (also block-diagonal) is reconstructed from the individual 
//   J-submatrices. This matrix is truncated at some energy [elim] or level [cb] (e.g. only cb columns of the 
//   matrix are used), and used to rotate the full IC Hamiltonian by: V(:,1:cb)'*Hic*V(:,1:cb), yielding a cb*cb
//   matrix which is diagonalised. The eigenvalues are the energies, and the eigenvector in the original 
//   basis is given by (Vr*Vrot')' if Vr are the eigenvectors in the rotated basis.
// --------------------------------------------------------------------------------------------------------------- //
iceig spectre_eig(icpars &pars, double elim)
{
   if(elim>=0) elim *= MEV2CM;

   // Calculates the Coulomb and Spin orbit Hamiltonian matrix in the |alpha,LSJ> basis
   sMat<double> Hcso = ic_Hcso(pars);

   // Finds indices to allow reordering the basis states so states with same J are grouped together
   unsigned int i, j, num_states=Hcso.nr(), ns=getdim(pars.n,pars.l), imax=num_states, cb=ns; 
   int ii=0, jj, imj=0, J2;
   fconf conf(pars.n,0,pars.l);
   std::vector< std::vector<int> > JsortM, JsortMmj; std::vector<int> JsortV, Jind, JsortVmj(ns,0);
   for(i=0; i<26; i++) JsortM.push_back(JsortV); JsortMmj = JsortM;
   for(i=0; i<conf.states.size(); i++) 
   { 
      J2=conf.states[i].J2; JsortM[J2].push_back(i); 
      for(j=0; j<(unsigned int)J2+1; j++) JsortMmj[J2].push_back(imj++); 
   }
   for(i=JsortM.size()-1; i<JsortM.size(); i--) 
      if(!JsortM[i].empty()) 
      {
         for(j=0; j<JsortM[i].size(); j++) { JsortV.push_back(JsortM[i][j]); Jind.push_back(i); }
         for(j=0; j<JsortMmj[i].size(); j++) JsortVmj[JsortMmj[i][j]]=ii++; 
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
      VEs.calc(isz,submats[i]);
      for(ii=0; ii<isz; ii++) { 
         for(int jj=0; jj<isz; jj++) Vunrot[(jj+idx)*num_states+ii+idx] = VEs.V(ii,jj); Eunrot[ii+idx] = VEs.E(ii); } 
      idx+=isz; 
   }
   for(i=0; i<JsortM.size(); i++) delete[]submats[i];
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

   // Construct the rotation matrix from the full eigenvector matrix, but putting the columns in Energy-order
   complexdouble *Vrot = new complexdouble[ns*ns];
   unsigned int iy=0, jy=0, id, J2i, J2j;
   for(i=0; i<num_states; i++) 
   { 
      jy=0; J2i = Jind[i];
      for(j=0; j<num_states; j++) 
      { 
         J2j = Jind[Esort[j]];
         if(fabs(Vunrot[Esort[j]*num_states+i])>SMALL) for(id=0; id<J2i+1; id++) Vrot[(iy+id)+(jy+id)*ns].r = Vunrot[Esort[j]*num_states+i];
         jy+=J2j+1;
         if(i==0) if(j==imax) cb=jy;
      } 
      iy+=J2i+1;
   }
   delete[]Vunrot;
   if(elim<0) cb = (unsigned int)elim;

   // Calculates the free-ion + crystal field Hamiltonian in the |alpha,LSJmJ> basis
   sMat<double> Hic,iHic; Hic = ic_hmltn(iHic,pars);
   sMat<double> zeroes; zeroes.zero(ns,ns);
   // Calculates the Zeeman term if magnetic field is not zero
   if(fabs(pars.Bx)>DBL_EPSILON || fabs(pars.By)>DBL_EPSILON || fabs(pars.Bz)>DBL_EPSILON)
   {
      std::vector<double> gjmbH(6,0.);
      if(fabs(pars.Bx)>DBL_EPSILON) { gjmbH[1]=-MUBc*pars.Bx; gjmbH[0]=GS*gjmbH[1]; }
      if(fabs(pars.By)>DBL_EPSILON) { gjmbH[3]=-MUBc*pars.By; gjmbH[2]=GS*gjmbH[3]; }
      if(fabs(pars.Bz)>DBL_EPSILON) { gjmbH[5]=-MUBc*pars.Bz; gjmbH[4]=GS*gjmbH[5]; }
      sMat<double> J,iJ; icmfmat mfmat(pars.n,pars.l,6,pars.save_matrices); mfmat.Jmat(J,iJ,gjmbH,pars.save_matrices); Hic+=J; iHic+=iJ;
   }
   Hic.reorder(JsortVmj); iHic.reorder(JsortVmj);

   // Rotates and truncates Hic by using the eigenvectors of Hc+Hso: Hic_r = V(:,1:cb)'*Hic*V(:,1:cb)
   complexdouble *zJmat,*Hrot,*zmt; zmt = new complexdouble[ns*cb]; Hrot = new complexdouble[cb*cb];
   char notranspose='N',transpose='C',uplo='U',side='L'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   int ins=(int)ns, icb=(int)cb;
   zJmat = zmat2f(Hic,iHic);
   F77NAME(zhemm)(&side,&uplo,&ins,&icb,&zalpha,zJmat,&ins,Vrot,&ins,&zbeta,zmt,&ins);
   F77NAME(zgemm)(&transpose,&notranspose,&icb,&icb,&ins,&zalpha,Vrot,&ins,zmt,&ins,&zbeta,Hrot,&icb); free(zJmat);

   // Diagonalises the truncated and rotated matrix
   iceig VEr; VEr.calc(icb,Hrot); delete[]Hrot;

   // Transform the eigenvectors of the rotated basis back into the original |alpha,LSJmJ> basis
   memset(zmt,0.,ns*cb*sizeof(complexdouble));
   complexdouble *Vr = new complexdouble[cb*ns*sizeof(complexdouble)]; 
   for(i=0; i<ns; i++) for(j=0; j<cb; j++) Vr[i*cb+j] = VEr.zV(i,j);
   F77NAME(zgemm)(&notranspose,&transpose,&icb,&ins,&icb,&zalpha,Vr,&icb,Vrot,&ins,&zbeta,zmt,&icb);  // (Vr*Vrot')'
   delete[]Vr;
   bool imagfl=false; for(i=0; i<ns; i++) for(j=0; j<cb; j++) if(zmt[i*cb+j].i>SMALL) { imagfl=true; break; }
   iceig VE(ns,!imagfl); double *En = VE.E(); for(i=0; i<cb; i++) En[i] = VEr.E(i); for(i=cb; i<ns; i++) En[i] = -DBL_MAX;
   if(imagfl) 
   {
      complexdouble *zV = VE.zV(0);
      for(i=0; i<ns; i++) for(j=0; j<cb; j++) { zV[j*ns+i].r = zmt[JsortVmj[i]*cb+j].r; zV[j*ns+i].i = -zmt[JsortVmj[i]*cb+j].i; }
   }
   else
   {
      double *V = VE.V(0); for(i=0; i<ns; i++) for(j=0; j<cb; j++) V[j*ns+i] = zmt[JsortVmj[i]*cb+j].r; 
   }

   return VE;
}

