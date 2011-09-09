/* so_cf.cpp
 *
 * Calculates the spin-orbit and crystal field interactions operator matrices, after Racah.
 *
 * Functions:
 *    sMat<double> racah_so(int n, double xi, orbital e_l=F);                   // Calculates the spin-orbit matrix
 *    sMat<double> racah_Umat(int n, int k, orbital e_l=F);                     // Calculates the reduced matrix U^k
 *    sMat<double> racah_ukq(int n, int k, int q, orbital e_l=F);               // Calculates the tensor operator U^k_q 
 *    sMat<double> fast_ukq(int n, int k, int q, orbital e_l=F);                // Calculates the tensor operator U^k_q 
 *        ("fast" algorithm... buggy! don't use!)
 *    sMat<double> racah_mumat(int n, int q, orbital e_l=F);                    // Calculates the magnetic moment operator
 *    void racah_mumat(int n,int q,sMat<double>&L,sMat<double>&S, orbital l=F); // Calculates the magnetic moment operator
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the spin orbit Hamilton matrix for a particular f^n configuration with strength xi
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_so(int n, double xi, orbital e_l)  // Defaults to f-electrons (see ic1ion.hpp - e_l=F by default)
{
   if(e_l!=S&&e_l!=P&&e_l!=D&&e_l!=F) {  std::cerr << "racah_so(): Only s-, p-, d- and f- configurations are implemented.\n"; exit(EXIT_FAILURE); }

   fconf conf(n,0,e_l);
   int num_states = (int)conf.states.size();
   std::vector<cfpls> cfpsi,cfpsj;
   int i,j,k,l,isz,jsz;
   double sumcfp;
   sMat<double> so(num_states,num_states);

   // Single electron configurations are programed directly in
   if(n==1)
   {
      int J2, S2, L2;
      S2 = conf.states[0].S2; L2 = 2*abs(conf.states[0].L);      // There is only 1 |LS> state in a 1-electron conf.
      double rmS=sqrt(3/2.), rmL=sqrt(e_l*(e_l+1)*(2*e_l+1));    // Reduce matrix elements, with s=1/2 substituted in.
      for(i=0; i<num_states; i++)                                // H_so is diagonal in J.
      {
         J2 = conf.states[i].J2;                                 // Equation 3-36 or 4-12 of Judd 1963
         so(i,i) = pow(-1.,(S2+L2+J2)/2.) * sixj(S2,S2,2,L2,L2,J2) * rmS * rmL;
      }
      return so*xi;
   }

   fconf confp(n-1,e_l);

   for(i=0; i<num_states; i++)
   {
      switch(e_l) {
         case P: cfpsi = racah_parents(n,conf.states[i].S2,conf.states[i].L); break;
         case D: cfpsi = racah_parents(n,conf.states[i].v,conf.states[i].S2,conf.states[i].L); break;
	 default:cfpsi = racah_parents(n,conf.states[i].v,conf.states[i].U,conf.states[i].S2,conf.states[i].L);  }
      for(j=i; j<num_states; j++)
      {
         if(abs(conf.states[i].S2-conf.states[j].S2)>2) continue;
         if(abs(abs(conf.states[i].L)-abs(conf.states[j].L))>1) continue;
         if(abs(conf.states[i].v-conf.states[j].v)==1) continue; // Selection rule due to group theory v-v'=0 or 2
         if(abs(conf.states[i].v-conf.states[j].v)>2) continue;  //    see Judd book, sections 8-3, and 8-4.
         if(conf.states[i].J2==conf.states[j].J2)
         {
            switch(e_l) {
               case P: cfpsj = racah_parents(n,conf.states[j].S2,conf.states[j].L); break;
               case D: cfpsj = racah_parents(n,conf.states[j].v,conf.states[j].S2,conf.states[j].L); break;
	       default:cfpsj = racah_parents(n,conf.states[j].v,conf.states[j].U,conf.states[j].S2,conf.states[j].L);  }
            sumcfp = 0.;
            isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(l=0; l<jsz; l++)
                  if(cfpsi[k].ind==cfpsj[l].ind)
//                   sumcfp += racahW(confp.states[cfpsi[k].ind].S2,conf.states[i].S2,1,2,1,conf.states[j].S2)
//                             * racahW(abs(confp.states[cfpsi[k].ind].L)*2,abs(conf.states[i].L)*2,2*e_l,2,2*e_l,abs(conf.states[j].L)*2)
//                             * cfpsi[k].cfp * cfpsj[l].cfp;
//          so(i,j) = -n*xi * racahW(conf.states[i].J2,abs(conf.states[i].L)*2,conf.states[j].S2,2,conf.states[i].S2,abs(conf.states[j].L)*2)
//                    * sqrt( (2.*abs(conf.states[i].L)+1.)*(2.*abs(conf.states[j].L)+1.)*(conf.states[i].S2+1.)*(conf.states[j].S2+1.) )
//                    * sqrt( (9./6)*e_l*(e_l+1)*(2*e_l+1) ) * sumcfp;
                     sumcfp += pow(-1.,abs(conf.states[i].L)+abs(confp.states[cfpsi[k].ind].L)+e_l + (conf.states[i].S2+confp.states[cfpsi[k].ind].S2+1.)/2.) 
                               * sixj(conf.states[i].S2,2,conf.states[j].S2,1,confp.states[cfpsi[k].ind].S2,1)
                               * sixj(abs(conf.states[i].L)*2,2,abs(conf.states[j].L)*2,2*e_l,abs(confp.states[cfpsi[k].ind].L)*2,2*e_l)
                               * cfpsi[k].cfp * cfpsj[l].cfp;
            so(i,j) = -n*xi * pow(-1.,abs(conf.states[i].L)+(conf.states[j].S2+conf.states[i].J2)/2.+1.) 
                      * sixj(conf.states[i].S2,2,conf.states[j].S2,abs(conf.states[j].L)*2,conf.states[i].J2,abs(conf.states[i].L)*2)
                      * sqrt( (2.*abs(conf.states[i].L)+1.)*(2.*abs(conf.states[j].L)+1.)*(conf.states[i].S2+1.)*(conf.states[j].S2+1.) )
                      * sqrt( (3./2)*e_l*(e_l+1.)*(2.*e_l+1.) ) * sumcfp;
            if(n>(2*e_l+1)) so(i,j) = -so(i,j);   // Phase factor difference for >half-filled shell. See Nielson/Koster or Racah III
            if(i!=j) so(j,i) = so(i,j);
         }
      }
   }
   return so;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the U^k reduced matrix which may be used to calculate the tensor operators U^k_q
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_Umat(int n, int k, orbital e_l)
{
   if(n==1) { sMat<double> U(1,1); U(0,0) = 1.; return U; }      // See Judd 1963, Eqn 5-13. with U^k=V^k/sqrt(2k+1)
// if(n==(4*e_l+1)) { sMat<double> U(1,1); U(0,0) = -1.; return U; }  // Error! Removed 21.11.10 After checking with Carnall paper.
   if(e_l!=P&&e_l!=D&&e_l!=F) { std::cerr << "racah_Umat(): Only p-, d- and f- configurations are implemented.\n"; exit(EXIT_FAILURE); }
   fconf conf(n,e_l);
   fconf confp(n-1,e_l);
   int num_states = (int)conf.states.size();
   std::vector<cfpls> cfpsi,cfpsj;
   int i,j,ii,jj,isz,jsz;
   double sumcfp,noncfpprod;
   sMat<double> U(num_states,num_states);
// if(k>e_l*2) return U;                      // Condition from reduced matrix element

   for(i=0; i<num_states; i++)
   {
      switch(e_l) {
         case P: cfpsi = racah_parents(n,conf.states[i].S2,conf.states[i].L); break;
         case D: cfpsi = racah_parents(n,conf.states[i].v,conf.states[i].S2,conf.states[i].L); break;
	 default:cfpsi = racah_parents(n,conf.states[i].v,conf.states[i].U,conf.states[i].S2,conf.states[i].L);  }
      for(j=i; j<num_states; j++)
      {
         if(abs(abs(conf.states[i].L)-abs(conf.states[j].L))>k) continue;
         if(conf.states[i].S2==conf.states[j].S2)
         {
            switch(e_l) {
               case P: cfpsj = racah_parents(n,conf.states[j].S2,conf.states[j].L); break;
               case D: cfpsj = racah_parents(n,conf.states[j].v,conf.states[j].S2,conf.states[j].L); break;
	       default:cfpsj = racah_parents(n,conf.states[j].v,conf.states[j].U,conf.states[j].S2,conf.states[j].L);  }
            sumcfp = 0.; 
//          noncfpprod = pow(-1.,-(double)e_l-abs(conf.states[j].L)) * sqrt( (2.*abs(conf.states[i].L)+1.)*(2.*abs(conf.states[j].L)+1.) );
            noncfpprod =                                               sqrt( (2.*abs(conf.states[i].L)+1.)*(2.*abs(conf.states[j].L)+1.) );
            isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(ii=0; ii<isz; ii++)
               for(jj=0; jj<jsz; jj++)
                  if(cfpsi[ii].ind==cfpsj[jj].ind)
//                   sumcfp += racahW(e_l*2,abs(conf.states[i].L)*2,e_l*2,abs(conf.states[j].L)*2,abs(confp.states[cfpsi[ii].ind].L)*2,k*2)
//                             * cfpsi[ii].cfp * cfpsj[jj].cfp * pow(-1.,(double)abs(confp.states[cfpsi[ii].ind].L)+k) * noncfpprod;
                     sumcfp += pow(-1.,(double)abs(confp.states[cfpsi[ii].ind].L)+abs(conf.states[i].L)+k+e_l)
                               * sixj(2*abs(conf.states[j].L),2*k,2*abs(conf.states[i].L),2*e_l,2*abs(confp.states[cfpsi[ii].ind].L),2*e_l)
                               * cfpsi[ii].cfp * cfpsj[jj].cfp * noncfpprod;
            if(fabs(sumcfp)!=0.) 
            {
             //if(n>(2*e_l+1)) sumcfp = -sumcfp; // Phase difference: [4l+2-n] = -(-1)^K [n] 
               U(i,j) = n * sumcfp;
               if(i!=j) U(j,i) = pow(-1.,abs(conf.states[i].L)-conf.states[i].S2/2.-abs(conf.states[j].L)+conf.states[j].S2/2.) * n * sumcfp;
            }
         }
      }
   }
 //char rmat[255]; strcpy(rmat,"results/ic1ion.umat"); rmzeros(U); mm_gout(U,rmat);
   return U;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix elements of the tensor operators U^k_q which describes (amongst other things) the CF.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_ukq(int n, int k, int q, orbital e_l)
{
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
 //int nn = n; if(nn>(2*e_l+1)) n = 4*e_l+2-n; 
   sMat<double> redmat = racah_Umat(n,k,e_l);
   int i,j,m,ns=0;
   int j2min,j2max;
   std::vector<int> L2,S2,J2,Jz2,irm;
   double rm;

   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L*2) - conf.states[i].S2);
      j2max =     abs(conf.states[i].L*2) + conf.states[i].S2;
      for(j=j2min; j<=j2max; j+=2)
      {
         for(m=-j; m<=j; m+=2)
         {
            L2.push_back(abs(conf.states[i].L*2));
            S2.push_back(    conf.states[i].S2 );
            irm.push_back(i);
            J2.push_back(j);
            Jz2.push_back(m);
	    ns++;
         }
      }
   }

   sMat<double> Ukq(ns,ns);

   for(i=0; i<ns; i++)
      for(j=0; j<ns; j++)
         if(S2[i]==S2[j])
         {
//          rm = pow(-1.,(S2[i]-L2[i]-J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * racahW(L2[i],J2[i],L2[j],J2[j],S2[i],2*k) * redmat(irm[i],irm[j]);
            rm = pow(-1.,(S2[i]+L2[j]+J2[i])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) *   sixj(J2[j],2*k,J2[i],L2[i],S2[i],L2[j]) * redmat(irm[i],irm[j]); // changed MR 26.1.10
//          rm = pow(-1.,(S2[i]-L2[i]-J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * racahW(L2[i],J2[i],L2[j],J2[j],S2[i],2*k) * redmat(irm[i],irm[j]);
          //if(nn>(2*e_l+1)) rm = -rm;
//          Ukq(i,j) = pow(-1.,(J2[i]+Jz2[i])/2.+k+q) * rm * wigner(J2[i],J2[j],0-Jz2[i],Jz2[j],2*k,-2*q) / sqrt(2.*k+1.);
//          Ukq(i,j) = pow(-1.,(J2[i]-Jz2[i])/2.) * threej(J2[j],2*k,J2[i],-Jz2[j],2*q,Jz2[i]) * rm; // changed MR 26.1.10
//          Ukq(i,j) = pow(-1.,(J2[i]+Jz2[i])/2.+k  ) * rm * wigner(J2[i],J2[j],0-Jz2[i],Jz2[j],2*k,-2*q) / sqrt(2.*k+1.);
            Ukq(i,j) = pow(-1.,(J2[i]-Jz2[i])/2.) * threej(J2[i],2*k,J2[j],-Jz2[i],2*q,Jz2[j]) * rm; // changed MR 26.1.10

          //if(i!=j) Ukq(j,i) = Ukq(i,j) * pow(-1.,(L2[i]-S2[i]-L2[j]+S2[j])/2.);
         }

   return Ukq;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix elements of the tensor operators U^k in the |LSJ> basis
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_uJ(int n, int k, orbital e_l)
{
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
// int nn = n; if(nn>(2*e_l+1)) n = 4*e_l+2-n; 
   sMat<double> redmat = racah_Umat(n,k,e_l);
   int i,j,ns=0;
   int j2min,j2max;
   std::vector<int> L2,S2,J2,irm;

   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L*2) - conf.states[i].S2);
      j2max =     abs(conf.states[i].L*2) + conf.states[i].S2;
      for(j=j2min; j<=j2max; j+=2)
      {
         L2.push_back(abs(conf.states[i].L*2));
         S2.push_back(    conf.states[i].S2 );
         irm.push_back(i);
         J2.push_back(j);
         ns++;
      }
   }

   sMat<double> Ukq(ns,ns);

   for(i=0; i<ns; i++)
      for(j=0; j<ns; j++)
         if(S2[i]==S2[j])
         {
//          Ukq(i,j) = pow(-1.,(S2[i]+L2[i]+J2[j])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * sixj(L2[i],J2[i],S2[i],J2[j],L2[j],2*k) * redmat(irm[i],irm[j]);
            Ukq(i,j) = pow(-1.,(S2[i]+L2[j]+J2[i])/2.+k) * sqrt((J2[i]+1.)*(J2[j]+1.)) * sixj(J2[j],2*k,J2[i],L2[i],S2[i],L2[j]) * redmat(irm[i],irm[j]);
         // if(nn>(2*e_l+1)) Ukq(i,j) = -Ukq(i,j);
         }

   return Ukq;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the matrix elements of U_k^q by a fast method -- currently buggy!
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> fast_ukq(int n, int k, int q, orbital e_l)
{
   // The function calculates first the reduce matrix elements <vULS||U^k||v'U'L'S> using the racah_Umat function. 
   // It then only calculates the subsequent J- and mJ-dependent matrix elements for the values of the reduced matrix
   // which is non-zero.  It does this by calculating a J-dependent only matrix <vULSJ||U^k||v'U'L'SJ'> for each 
   // non-zero reduced matrix element, and then an mJ-dependent matrix <vULSJmJ|Ukq|v'U'L'SJ'mJ> for each J-J' value. 
   // The function then slots each of these smaller "J-blocks" into the full Ukq matrix leaving the rest zeros.
   // In this way it avoids recalculating the J/mJ-dependent elements many times, resulting in faster computation.
   //
   // The equations to calculate the matrix elements may be found in Elliot, Judd, Runciman, Proc. R. Soc. Lon. A
   // vol. 240, no. 1223, pp. 509-523, equations 25(mJ-dependent part), 26(J-dependent part), 27(reduced m.e.).

   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
   sMat<double> redmat = racah_Umat(n,k,e_l);
   std::vector< std::vector<int> > nz = redmat.find();
   int sz = (int)nz.size();
   int i,j,ns,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int j2min,j2max,j2pmin,j2pmax,count,countp;
   int rmi,rmj,L2,L2p,S2,S2p,J2,J2p,Jz2,Jz2p,imJ_i,imJ_j;
   std::vector<int> indexJstart,indexJstop;
   sMat<double> mJmat_i(1,1);
   sMat<double> rmJ;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   double rm;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   ns = 0; minJ2 = 99; maxJ2 = 0;
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      if(j2min<minJ2) minJ2 = j2min; if(j2max>maxJ2) maxJ2 = j2max;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2)
         ns += j+1;
      indexJstop.push_back(ns);
   }

   sMat<double> Ukq(ns,ns);

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // Only calculates for the non-zero values of the reduced matrix elements <vULS||U^k||v'U'L'S>
   for(i=0; i<sz; i++)
   {
      rmi = nz[i][0]; rmj = nz[i][1];
      L2 = abs(conf.states[rmi].L)*2; L2p = abs(conf.states[rmj].L)*2; 
      S2 = conf.states[rmi].S2; S2p = conf.states[rmj].S2;
      if(S2!=S2p) continue;
      
      // Caculate the J-dependent matrix elements <J||U^k||J'> = <vULSJ||U^k||v'U'L'SJ'>/<vULS||U^k||v'U'L'S>
      j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2); j2pmax = L2p+S2;
      rmJ.clear(); count = 0;
      for(J2=j2min; J2<=j2max; J2+=2)
      {
         countp = 0; valJ_i = (J2-minJ2)/2;
         for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
         {
            rm = pow(-1.,(S2-L2-J2p)/2.+k) * sqrt((J2+1.)*(J2p+1.)) * racahW(L2,J2,L2p,J2p,S2,2*k);

            // Calculates the q- and Jz- dependent matrix elements <mJ|Ukq|mJ'> = <psi|Ukq|psi'>/<J||U^k||J'>
            valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {  
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p;
                     mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,(J2+Jz2)/2.+k+q) * wigner(J2,J2p,0-Jz2,Jz2p,2*k,-2*q) / sqrt(2.*k+1.);
                  }
               }
            }
            rmJ.pset(count+1,count+J2+1,countp+1,countp+J2p+1,mJmat[valJ_i][valJ_j]*rm);
            countp += J2p+1;
         }
         count += J2+1;
      }
      Ukq.pset(indexJstart[rmi],indexJstop[rmi],indexJstart[rmj],indexJstop[rmj],rmJ*redmat(rmi,rmj));
   }

   return Ukq;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the magnetic moment operator matrix mu_{x,y,z} = L + gS, in the |vSLJM> basis
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> racah_mumat(int n, int q, orbital e_l)
{
   const double g_s = 2.0023193043622;   // electronic g-factor, from NIST - http://physics.nist.gov/cuu/Constants/

   if(e_l!=S && e_l!=P && e_l!=D && e_l!=F) { std::cerr << "racah_mumat(): Only s-, p-, d- and f- configurations are implemented.\n"; }
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
   int i,j,ns,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int j2min,j2max,j2pmin,j2pmax,count,countp;
   int L2,L2p,S2,S2p,J2,J2p,Jz2,Jz2p,imJ_i,imJ_j;
   std::vector<int> indexJstart,indexJstop;
   sMat<double> mJmat_i(1,1);
   sMat<double> rmJ;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   double Lrm,Srm,rm;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   ns = 0; minJ2 = 99; maxJ2 = 0;
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      if(j2min<minJ2) minJ2 = j2min; if(j2max>maxJ2) maxJ2 = j2max;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2)
         ns += j+1;
      indexJstop.push_back(ns);
   }

   sMat<double> muiq(ns,ns);

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // Calculates the matrix for the operator L+gS for a particular q (q=0 is z, q=+/-1 is linear combination of x or y)
   for(i=0; i<num_states; i++)
    for(j=0; j<num_states; j++)
    {
      L2 = abs(conf.states[i].L)*2; L2p = abs(conf.states[j].L)*2; S2 = conf.states[i].S2; S2p = conf.states[j].S2;

      if(S2!=S2p || L2!=L2p) continue;
      if(e_l>1 && conf.states[i].v!=conf.states[j].v) continue;
      if(e_l>2 && conf.states[i].U!=conf.states[j].U) continue;
      
      // Caculate the J-dependent reduced matrix elements 
      j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2); j2pmax = L2p+S2;
      rmJ.clear(); count = 0;
      for(J2=j2min; J2<=j2max; J2+=2)
      {
         countp = 0; valJ_i = (J2-minJ2)/2;
         for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
         {
	    Lrm = pow(-1.,(S2p+L2p+J2)/2.)  * sqrt((L2p+1.)*(J2+1.)*(J2p+1.)*(L2p/2.)*(L2p/2.+1.)) * sixj(L2,J2,S2p,J2p,L2p,2);
	    Srm = pow(-1.,(S2p+L2p+J2p)/2.) * sqrt((S2p+1.)*(J2+1.)*(J2p+1.)*(S2p/2.)*(S2p/2.+1.)) * sixj(S2,J2,L2p,J2p,S2p,2);
//          Lrm /= pow(-1.,L2p/2.)*(L2p+1.); Srm /= pow(-1.,S2p/2.)*(S2p+1.);
//          Lrm /= pow(-1.,L2p/2.); Srm /= pow(-1.,S2p/2.);
//          Lrm /= (L2p+1.); Srm /= (S2p+1.);
            rm  = -( Lrm + g_s*Srm );

            // Calculates the q- and Jz- dependent matrix elements as given by the Wigner-Eckart theorem
            valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {  
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p;
                     mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,(J2-Jz2)/2.) * threej(J2,2,J2p,-Jz2,2*q,Jz2p);
                  }
               }
            }
            rmJ.pset(count+1,count+J2+1,countp+1,countp+J2p+1,mJmat[valJ_i][valJ_j]*rm);
            countp += J2p+1;
         }
         count += J2+1;
      }
      muiq.pset(indexJstart[i],indexJstop[i],indexJstart[j],indexJstop[j],rmJ);
    }

   return muiq;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the magnetic moment operator as separate matrices, L_{x,y,z} and S_{x,y,z}, in the |vSLJM> basis
// --------------------------------------------------------------------------------------------------------------- //
void racah_mumat(int n, int q, sMat<double> &L1q, sMat<double> &S1q, orbital e_l)
{
   if(e_l!=S && e_l!=P && e_l!=D && e_l!=F) { std::cerr << "racah_mumat(): Only s-, p-, d- and f- configurations are implemented.\n"; }
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
   int i,j,ns,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int j2min,j2max,j2pmin,j2pmax,count,countp;
   int L2,L2p,S2,S2p,J2,J2p,Jz2,Jz2p,imJ_i,imJ_j;
   std::vector<int> indexJstart,indexJstop;
   sMat<double> mJmat_i(1,1);
   sMat<double> rmL, rmS;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   double Lrm,Srm;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   ns = 0; minJ2 = 99; maxJ2 = 0;
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      if(j2min<minJ2) minJ2 = j2min; if(j2max>maxJ2) maxJ2 = j2max;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2)
         ns += j+1;
      indexJstop.push_back(ns);
   }

   L1q.zero(ns,ns); S1q.zero(ns,ns);

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // Calculates the matrix for the operator L+gS for a particular q (q=0 is z, q=+/-1 is linear combination of x or y)
   for(i=0; i<num_states; i++)
    for(j=0; j<num_states; j++)
    {
      L2 = abs(conf.states[i].L)*2; L2p = abs(conf.states[j].L)*2; S2 = conf.states[i].S2; S2p = conf.states[j].S2;

      if(S2!=S2p || L2!=L2p) continue;
      if(e_l>1 && conf.states[i].v!=conf.states[j].v) continue;
      if(e_l>2 && conf.states[i].U!=conf.states[j].U) continue;
      
      // Caculate the J-dependent reduced matrix elements 
      j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2); j2pmax = L2p+S2;
      rmL.clear(); rmS.clear(); count = 0;
      for(J2=j2min; J2<=j2max; J2+=2)
      {
         countp = 0; valJ_i = (J2-minJ2)/2;
         for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
         {
            Lrm = ( pow(-1.,(S2p+L2p+J2)/2.)  * sqrt((L2p+1.)*(J2+1.)*(J2p+1.)*(L2p/2.)*(L2p/2.+1.)) * sixj(J2p,2,J2,L2,S2,L2p) );
            Srm = ( pow(-1.,(S2p+L2p+J2p)/2.) * sqrt((S2p+1.)*(J2+1.)*(J2p+1.)*(S2p/2.)*(S2p/2.+1.)) * sixj(J2p,2,J2,S2,L2,S2p) );
//          Lrm /= pow(-1.,L2p/2.)*(L2p+1.); Srm /= pow(-1.,S2p/2.)*(S2p+1.);
//          Lrm /= pow(-1.,L2p/2.); Srm /= pow(-1.,S2p/2.);
//          Lrm /= (L2p+1.); Srm /= (S2p+1.);

            // Calculates the q- and Jz- dependent matrix elements as given by the Wigner-Eckart theorem
            valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {  
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p;
                     mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,(J2-Jz2)/2.) * threej(J2,2,J2p,-Jz2,2*q,Jz2p);
                  }
               }
            }
            rmL.pset(count+1,count+J2+1,countp+1,countp+J2p+1,mJmat[valJ_i][valJ_j]*Lrm);
            rmS.pset(count+1,count+J2+1,countp+1,countp+J2p+1,mJmat[valJ_i][valJ_j]*Srm);
            countp += J2p+1;
         }
         count += J2+1;
      }
      L1q.pset(indexJstart[i],indexJstop[i],indexJstart[j],indexJstop[j],rmL);
      S1q.pset(indexJstart[i],indexJstop[i],indexJstart[j],indexJstop[j],rmS);
    }
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the magnetic moment operator as separate matrices, L_{x,y,z} and S_{x,y,z}, in the |vSLJM> basis
// --------------------------------------------------------------------------------------------------------------- //
void chanlam_mumat(int n, int q, sMat<double> &mu, orbital e_l)
{
   if(e_l!=S && e_l!=P && e_l!=D && e_l!=F) { std::cerr << "racah_mumat(): Only s-, p-, d- and f- configurations are implemented.\n"; }
   fconf conf(n,e_l);
   int num_states = (int)conf.states.size();
   int i,j,ns,k,j2min,j2max,L2,L2p,S2,S2p;
   std::vector<int> index, J2, Jz2;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   ns = 0;
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      for(j=j2min; j<=j2max; j+=2)
         for(k=-j; k<=j; k+=2) { Jz2.push_back(k); J2.push_back(j); index.push_back(i); ns++; }
   }

   mu.zero(ns,ns);

   double M,J,denom,f,fp,g; double g_s = 2.0023193043622; // electronic g-factor
   // Calculates the matrix for the operator L+gS for a particular q (q=0 is z, q=+/-1 is linear combination of x or y)
   for(i=0; i<ns; i++)
      for(j=0; j<ns; j++)
      {
         L2 = abs(conf.states[index[i]].L)*2; L2p = abs(conf.states[index[j]].L)*2; S2 = conf.states[index[i]].S2; S2p = conf.states[index[j]].S2;

         if(S2!=S2p || L2!=L2p) continue;
         if(e_l>1 && conf.states[i].v!=conf.states[j].v) continue;
         if(e_l>2 && conf.states[i].U!=conf.states[j].U) continue;
      
         J = J2[i]/2.; M=Jz2[i]/2.; denom = (J*J*(2*J+1.)*(2*J-1.));
         if(fabs(denom)<DBL_EPSILON) f=0.; else f = sqrt( ((S2+L2+2*J)/2.+1)*((S2+L2-2*J)/2.+1)*((S2+2*J-L2)/2.)*((L2+2*J-S2)/2.) / denom );
         denom = ((J+1)*(J+1)*(2*(J+1)+1)*(2*(J+1)-1));
         if(fabs(denom)<DBL_EPSILON) fp=0.; else fp = sqrt( ((S2+L2)/2.+(J+1)+1)*((S2+L2)/2.-(J+1)+1)*((S2-L2)/2.+(J+1))*((L2-S2)/2.+(J+1)) / denom );
	 if (J2[i]!=0) g = 1 + (g_s-1) * (J*(J+1) - (L2/2.)*((L2/2.)+1) + (S2/2.)*((S2/2.)+1)) / (2*J*(J+1)); else g = 0;

         if(J2[i]==J2[j])
         {
           if (Jz2[i]==Jz2[j] && q==3)
             mu(i,j) = M*g;                                           // mu_z
           else if(Jz2[j]==(Jz2[i]+2))
           {
             if(q==1)      mu(i,j) = sqrt( (J+M+1)*(J-M) ) * g/2;     // mu_x 
             else if(q==2) mu(i,j) = sqrt(/*-*/(J+M+1)*(J-M) ) * g/2; // mu_y
           }
           else if(Jz2[j]==(Jz2[i]-2))
           {
             if(q==1)      mu(i,j) = sqrt( (J-M+1)*(J+M) ) * g/2;
             else if(q==2) mu(i,j) =-sqrt(/*-*/(J-M+1)*(J+M) ) * g/2;
           }
         }
         else if(J2[j]==(J2[i]-2))
         {
           if(Jz2[i]==Jz2[j] && q==3)
             mu(i,j) = (g_s-1) * sqrt(J*J-M*M) * f/2;
           else if(Jz2[j]==(Jz2[i]+2))
           {
             if(q==1)      mu(i,j) = (g_s-1) * sqrt( (J-M-1)*(J-M) ) * f/4;
             else if(q==2) mu(i,j) = (g_s-1) * sqrt(/*-*/(J-M-1)*(J-M) ) * f/4;
           }
           else if(Jz2[j]==(Jz2[i]-2))
           {
             if(q==1)      mu(i,j) =-(g_s-1) * sqrt( (J+M-1)*(J+M) ) * f/4;
             else if(q==2) mu(i,j) = (g_s-1) * sqrt(/*-*/(J+M-1)*(J+M) ) * f/4;
           }
         }
         else if(J2[j]==(J2[i]+2))
         {
           if(Jz2[i]==Jz2[j] && q==3)
             mu(i,j) = (g_s-1) * sqrt( (J+M+1)*(J-M+1) ) * fp/2;
           else if(Jz2[j]==(Jz2[i]+2))
           {
             if(q==1)      mu(i,j) =-(g_s-1) * sqrt( (J+M+1)*(J+M+2) ) * fp/4;
             else if(q==2) mu(i,j) =-(g_s-1) * sqrt(/*-*/(J+M+1)*(J+M+2) ) * fp/4;
           }
           else if(Jz2[j]==(Jz2[i]-2))
           {
             if(q==1)      mu(i,j) = (g_s-1) * sqrt( (J-M+1)*(J-M+2) ) * fp/4;
             else if(q==2) mu(i,j) =-(g_s-1) * sqrt(/*-*/(J-M+1)*(J-M+2) ) * fp/4;
           }
         }
      }
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ states.cpp cfp.cpp njsyms.cpp so_cf.cpp; ./a.out
// --------------------------------------------------------------------------------------------------------------- //
/*int main(int argc, char *argv[])
{
   int n,k,q;
   if(argc>1) n = atoi(argv[1]); else n = 2;
   if(argc>2) k = atoi(argv[2]); else k = 2;
   if(argc>3) q = atoi(argv[3]); else q = 0;

 //sMat<double> so = racah_so(n,1.);
 //std::cout << "SO matrix is:\n" << so;

   sMat<double> U2 = racah_Umat(n,6);
 //std::cout << "U^2 matrix is:\n" << U2;
 //std::cout << "x=" << U2.display_full();

 //sMat<double> U20 = fast_ukq(n,4,0);
   sMat<double> U20 = racah_ukq(n,k,q);
 //std::cout << "U^2_0 matrix is:\n" << U20;
 //std::cout << "x=" << U20.display_full();
   sMat<double> V20 = fast_ukq(n,k,q);
   sMat<double> comp = U20 - V20;
   std::cout << "sum(sum(abs( racah_ukq(" << n << "," << k << "," << q << ") - ";
   std::cout << "fast_ukq(" << n << "," << k << "," << q << ") ))) = ";
   std::cout << vsum(msum(mabs(comp))) << "\n"; // Checks that racah_ukq and fast_ukq agree!
 //std::cout << V20.display_full();
 //std::cout << "comp =\n" << comp;
 //std::cout << "\nfast_ukq =\n" << V20;

   return 0;
}*/
