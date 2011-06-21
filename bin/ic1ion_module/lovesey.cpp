/* lovesey.cpp
 *
 * Function:  bool lovesey_akk(aKK,K,Kp,n,l)                    //  Calculates the a(K,K') matrix
 *            bool lovesey_ckk(cKK,K,Kp,n,l)                    //  Calculates the c(K,K') matrix
 *            void lovesey_Qq(Qq,iQq,q,n,l,Jvec)                //  Calculates the scattering operator Q_q
 *            complexdouble spherical_harmonics(k,q,theta,phi)  //  Calculates Yq^k(theta,phi)
 *
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2009 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 * Ref: E. Balcar and S.W. Lovesey, Theory of Magnetic Neutron and Photon Scattering, OUP, 1989
 * The formulae for the matrices are reproduced below, but may also be obtained from equations 9.8-9.16 of
 * Lovesey and Rimmer, Rep. Prog. Phys. v32, p333-394, 1969
 *
 */

#include "ic1ion.hpp"
#include <fstream>

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the a(K,K') matrix - the routine returns a true if the matrix is real, and false if it is imaginary
//    The matrix is returned in the argument aKK
// --------------------------------------------------------------------------------------------------------------- //
bool lovesey_aKK(sMat<double> &aKK, int K, int Kp, int n, orbital l)
{
   // Loads a previously save matrix if it exists
   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
 //nstr[0] = (l==F?102:100); 
   switch(l) { case S: nstr[0]=115; break; case P: nstr[0]=112; break; case D: nstr[0]=100; break; default: nstr[0]=102; }
   if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 65;   // 65 is ASCII for "A", 100=="d" and 102=="f"
   nstr[1] = K+48; nstr[2] = Kp+48; nstr[3] = 0; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
   aKK = mm_gin(filename); if(!aKK.isempty()) return 1;

   if(l!=S&&l!=P&&l!=D&&l!=F) {  std::cerr << "lovesey_cKK(): Only s-, p-, d- and f- configurations are implemented.\n"; exit(EXIT_FAILURE); }
   int np = (n==1)?n:(n-1); 
   fconf confp(np,l);
   fconf conf(n,l);
   int num_states = (int)conf.states.size(), ns=0;
   int i,j,k,kk,v,vp,isz,jsz,iJ,jJ;
   int j2min,j2max,j2pmin,j2pmax;
   int L2,L2p,S2,S2p,J2,J2p,L2b;//,S2b;
   std::vector<int> indexJstart,indexJstop;
   std::vector<cfpls> cfpsi,cfpsj;
   double rmLS,sumcfp;

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2) ns++;
      indexJstop.push_back(ns);
   }

   aKK.zero(ns,ns);

   // Selection rules on K.
   if(K<0  || K>(2*abs(l))   || (K%2)!=0)      return 1;        // Eqn 4.2.2
   if(Kp<1 || Kp>=(2*abs(l)) || ((Kp+1)%2)!=0) return 1;        // Eqn 4.2.3
   if(K!=(Kp+1) && K!=(Kp-1))                  return 1;        // Eqn 4.2.4  (first 3j symbol, needs 1+K+K' even)

   // Calculates the matrix elements:
   // Eqn. 3.6.8 (d==delta)
   //
   //                                  K'+1           J'+S+L+L'+l  1/2    2                  1/2
   // A(K,K') = ( <j    > + <j    > ) i      d    (-1)            2    [l]  [L,L',K',K',K,J']
   //               K'-1      K'+1            SS'
   //
   //              \/   ( 1 K K' )  { 1  1  1  }  { J' K' J  } 
   //              /\   ( 0 0 0  )  { K' K  K' }  { L  S  L' }  A(K',K',l)
   //                                            _
   //                     ---     _   _          L                                _
   //              \/   n >   (t{|t) (t|}t') (-1)  { L' K' L }                Lb==L
   //              /\     ---                      { l  Lb l }
   //                      t

   // Where: (Eqn 3.4.9)
   //                                    1/2
   //                 l+1 [ (2l+3)(l+1) ]    ( l K' l+1 )  { 1  K'  K }
   // A(K,K',l) = (-1)    [ ----------- ]    ( 0 0  0   )  { l l+1  l }
   //                     [   (2l+1)    ]
   //
   
   // Calculate the coefficient A(K,K',l)
   double AKKl = pow(-1.,l+1.) * sqrt( (2.*l+3.)*(l+1.)/(2.*l+1.) ) * threej(2*l,2*Kp,2*(l+1),0,0,0) * sixj(2,2*Kp,2*Kp,2*l,2*(l+1),2*l);

   // Calculate the none state dependent part of the matrix elements
   double redmat = sqrt(2) * (2*l+1)*(2*l+1) * (2*Kp+1.) * sqrt(2*K+1.) * threej(2,2*K,2*Kp,0,0,0) * sixj(2,2,2,2*Kp,2*K,2*Kp) * AKKl;

   // Determines the factor i^(K'-1) - because K' is always odd, matrix is always real
   if((Kp-1)%4==0) redmat = -redmat;

   // Calculates the matrix elements at particular |J,M>, |J',M'>
   for(i=0; i<num_states; i++)
   {
      L2 = abs(conf.states[i].L)*2; S2 = conf.states[i].S2; v = conf.states[i].v;
      if(n!=1) 
      {
         switch(l) {
            case P: cfpsi = racah_parents(n,S2,conf.states[i].L); break;
            case D: cfpsi = racah_parents(n,v,S2,conf.states[i].L); break;
            default:cfpsi = racah_parents(n,v,conf.states[i].U,S2,conf.states[i].L);  }
      }

      for(j=0; j<num_states; j++)
      {
         L2p = abs(conf.states[j].L)*2; S2p = conf.states[j].S2; vp = conf.states[j].v;

         if(S2!=S2p) continue;                                  // delta_SS'
         if(abs(L2-L2p)>(2*Kp)) continue;                       // Triangular condition on 6j symbol in 3.6.9

	 if(n==1) 
            rmLS = -redmat;
         else
         {
            switch(l) {
               case P: cfpsj = racah_parents(n,S2p,conf.states[j].L); break;
               case D: cfpsj = racah_parents(n,vp,S2p,conf.states[j].L); break;
               default:cfpsj = racah_parents(n,vp,conf.states[j].U,S2p,conf.states[j].L);  }

            sumcfp = 0.; isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(kk=0; kk<jsz; kk++)
                  if(cfpsi[k].ind==cfpsj[kk].ind) {
                     L2b = 2*abs(confp.states[cfpsi[k].ind].L); //S2b = confp.states[cfpsi[k].ind].S2;
                     sumcfp += cfpsi[k].cfp*cfpsj[kk].cfp * pow(-1.,L2b/2) * sixj(L2p,2*Kp,L2,2*l,L2b,2*l); }
            rmLS = sqrt( (L2+1.)*(L2p+1.) ) * n * sumcfp * redmat;
         }

         // Caculate the J-dependent reduced matrix elements 
         j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2p); j2pmax = L2p+S2p;
         iJ=indexJstart[i]-1; jJ=indexJstart[j]-1;
         for(J2=j2min; J2<=j2max; J2+=2)
         {
            for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
            {
               if(abs(J2-J2p)<=(2*Kp))                          // Triangular condition on 6j symbol in 3.6.9
                  aKK(iJ,jJ) = pow(-1.,(J2p+S2+L2+L2p)/2.+l) * sqrt(J2p+1.) * sixj(J2p,2*Kp,J2,L2,S2,L2p) * rmLS;
               jJ++;
            }
            iJ++; jJ=indexJstart[j]-1;
         }
      }
   }

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out); FILEOUT.close();
   rmzeros(aKK); mm_gout(aKK,filename); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the c(K,K') matrix, returning true if matrix is real, false if imaginary
// --------------------------------------------------------------------------------------------------------------- //
bool lovesey_cKK(sMat<double> &cKK, int K, int Kp, int n, orbital l)
{
   // Loads a previously save matrix if it exists
   char nstr[6]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
 //nstr[0] = (l==F?102:100); 
   switch(l) { case S: nstr[0]=115; break; case P: nstr[0]=112; break; case D: nstr[0]=100; break; default: nstr[0]=102; }
   if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 67;   // 67 is ASCII for "C", 100=="d" and 102=="f"
   nstr[1] = K+48; nstr[2] = Kp+48; nstr[3] = 0; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
   cKK = mm_gin(filename); if(!cKK.isempty()) return 1;

   if(l!=S&&l!=P&&l!=D&&l!=F) {  std::cerr << "lovesey_cKK(): Only s-, p-, d- and f- configurations are implemented.\n"; exit(EXIT_FAILURE); }
   int np = (n==1)?n:(n-1);
   fconf confp(np,l);
   fconf conf(n,l);
   int num_states = (int)conf.states.size(), ns=0;
   int i,j,k,kk,j2min,j2max,j2pmin,j2pmax;
   int L2,L2p,S2,S2p,J2,J2p;
   std::vector<cfpls> cfpsi,cfpsj;
   std::vector<int> indexJstart,indexJstop,irm;
   sMat<double> mJmat_i(1,1);
   sMat<double> rmJ;
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   double rmLS,sumcfp;
   int v,vp,isz,jsz,L2b,S2b,iJ,jJ; 

   // Eqn. 3.6.11
   //
   //                  K      1/2                            1/2      1/2+S'+L'
   // C(K,K') = <j >  i  (1/2)     [l,l,S,S',L,L',J',K,K',K']     (-1)
   //             K                                                       _ _
   //                              { 1  K  K' }    ---     _   _          S+L
   //              \/   ( l K l )  { S' L' J' }  n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }
   //              /\   ( 0 0 0 )  { S  L  J  }    ---                        { s  Sb s  }  { l  Lb l  }
   //                                               t

   // Determines the L S J values for each matrix elements and the index of each J-J' block
   for(i=0; i<num_states; i++)
   {
      j2min = abs(abs(conf.states[i].L)*2-conf.states[i].S2); j2max = abs(conf.states[i].L)*2+conf.states[i].S2;
      indexJstart.push_back(ns+1);
      for(j=j2min; j<=j2max; j+=2) { ns++; irm.push_back(i); }
      indexJstop.push_back(ns);
   }

   cKK.zero(ns,ns);

   // The triangular conditions require that K<=2l and K=even (from 3j), and that (K-1)<=K'<=(K+1) (from 9j)
   if(K%2!=0 || abs(Kp-K)>1) return 1;

   // Calculates part of the prefactor which only depends on K,K',l
   double redmat = sqrt(1./2) * (2*l+1.) * sqrt(2*K+1.) * (2*Kp+1.) * threej(2*l,2*K,2*l,0,0,0); 
   if(fabs(redmat)<DBL_EPSILON) return 1;


   // Determines the factor i^K - because K is always even, matrix is always real
   if(K%4==2) redmat = -redmat;

   rmJ.zero(num_states,num_states);
   // Calculates the matrix elements at particular |J,M>, |J',M'>
   for(i=0; i<num_states; i++)
   {
      L2 = abs(conf.states[i].L)*2; S2 = conf.states[i].S2; v = conf.states[i].v;
      if(n!=1) 
      {
         switch(l) {
            case P: cfpsi = racah_parents(n,S2,conf.states[i].L); break;
            case D: cfpsi = racah_parents(n,v,S2,conf.states[i].L); break;
            default:cfpsi = racah_parents(n,v,conf.states[i].U,S2,conf.states[i].L);  }
      }

      for(j=0; j<num_states; j++)
      {
         L2p = abs(conf.states[j].L)*2; S2p = conf.states[j].S2; vp = conf.states[j].v;

         if(abs(S2-S2p)>2)     continue;                        // Triangular condition on 6j symbol in 3.6.11
         if(abs(L2-L2p)>(2*K)) continue;                        // Triangular condition on 6j symbol in 3.6.11

	 if(n==1)
	    rmLS = redmat;
	 else
	 {
            switch(l) {
               case P: cfpsj = racah_parents(n,S2p,conf.states[j].L); break;
               case D: cfpsj = racah_parents(n,vp,S2p,conf.states[j].L); break;
               default:cfpsj = racah_parents(n,vp,conf.states[j].U,S2p,conf.states[j].L);  }

            sumcfp = 0.; isz = (int)cfpsi.size(); jsz = (int)cfpsj.size();
            for(k=0; k<isz; k++)
               for(kk=0; kk<jsz; kk++)
                  if(cfpsi[k].ind==cfpsj[kk].ind) {
                     L2b = 2*abs(confp.states[cfpsi[k].ind].L); S2b = confp.states[cfpsi[k].ind].S2;
                     sumcfp += cfpsi[k].cfp*cfpsj[kk].cfp * pow(-1.,(1+S2p+L2p+L2b+S2b)/2.) 
                                  * sixj(S2,2,S2p,1,S2b,1)*sixj(L2,2*K,L2p,2*l,L2b,2*l); }

            rmLS = sqrt( (S2+1.)*(S2p+1.)*(L2+1.)*(L2p+1.) ) * n * sumcfp * redmat;
         }

         // Caculate the J-dependent reduced matrix elements 
         j2min = abs(L2-S2); j2max = L2+S2; j2pmin = abs(L2p-S2p); j2pmax = L2p+S2p;
         iJ=indexJstart[i]-1; jJ=indexJstart[j]-1;
         for(J2=j2min; J2<=j2max; J2+=2)
         {
            for(J2p=j2pmin; J2p<=j2pmax; J2p+=2)
            {
               if(abs(J2-J2p)<=(2*Kp))                          // Triangular condition on 9j symbol in 3.6.11
                  cKK(iJ,jJ) = sqrt(J2p+1.) * ninej(2,2*K,2*Kp,S2p,L2p,J2p,S2,L2,J2) * rmLS;
               jJ++;
            }
            iJ++; jJ=indexJstart[j]-1;
         }
      }
   }

   std::fstream FILEOUT; FILEOUT.open(filename, std::fstream::out); FILEOUT.close();
   rmzeros(cKK); mm_gout(cKK,filename); return 1;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the transition operator Q_q
// --------------------------------------------------------------------------------------------------------------- //
void lovesey_Qq(std::vector< sMat<double> > &Qq, int q, int n, orbital l, std::vector<double> &Jvec)
{
   int i,j,ns=getdim(n,l);                                      // Number of states = ^{4l+2}C_{n}

   std::string errormsg("lovesey_Qq(): Unable to calculate the A(K,K') or B(K,K') matrix\n");
   double theta = Jvec[0], phi = Jvec[1], J[]={Jvec[2],0,Jvec[3],0,Jvec[4],0,Jvec[5]};
   int K,Kp,Q,Qp;
   sMat<double> Akk,Bkk,ckk,Qmat,QLmat,QSmat;
   Qmat.zero(ns,ns); Qq.clear(); for(i=0; i<6; i++) Qq.push_back(Qmat);

   bool iA,ic;
   complexdouble Ykq; double Tj,Tj2;

   int k,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int imJ_i,imJ_j,J2,J2p,Jz2,Jz2p;
   std::vector< std::vector<int> > nz;

   sMat<double> mJmat_i(1,1);
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   std::vector<int> iJst,iJen;
   fconf conf(n,0,l); minJ2=99; maxJ2=0; k=0;
   for(i=0; i<(int)conf.states.size(); i++) {
      j = conf.states[i].J2; if(j<minJ2) minJ2 = j; if(j>maxJ2) maxJ2 = j; iJst.push_back(k+1); k+=j+1; iJen.push_back(k); }

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // The scattering operator is given by: (eqn. 3.7.9 and 3.7.11 for the orbital and spin parts respectively)
   //
   //              ____ ---  {  K'-1 ^    2K'+1                                                  K'                      }
   // <i|Q |j> = \/4*pi >    { Y    (k) [ ------ ] [ A(K'-1,K') + B(K'-1,K') ] ( K'-1 K' 1  ) + Y  B(K',K') ( K' K' 1  ) }
   //     q             ---  {  Q    -     K'+1                                ( Q    Q' -q )    Q          ( Q  Q' -q ) }
   //                  K'QQ'
   //                     -J'+q+M  ( K' J'  J )
   //                 (-1)      *  ( Q' M' -M )                  with |i>==|avUSLJM> and |j>==|a'v'U'S'L'J'M'>
   //
   // The matrices A(K,K') and B(K,K') are:                                                              1/2 
   //                                                                K+2 [                        ( K+1 )            ]
   //  A(K,K') = ( <j    > + <j    > ) a(K,K')           B(K,K+1) = ---- [ <j > c(K,K+1) + <j   > ( --- ) c(K+2,K+1) ]
   //                K'-1      K'+1                                 2K+3 [   K               K+2  ( K+2 )            ]
   //
   //                                                                                                    1/2
   //                                                                K-1 [                        (  K  )            ]
   //  B(K,K) = <j > c(K,K)                              B(K,K-1) = ---- [ <j > c(K,K-1) + <j   > ( --- ) c(K-2,K-1) ]
   //             K                                                 2K-1 [   K               K-2  ( K-1 )            ]
   //
   //
   //  The 3j symbols with the matrices a(K,K') and c(K,K') restrict K to even integers up to 2l, and K' to odd integers
   //  less than 2l. I.e. for f-electrons, K=0,2,4,6 and K'=1,3,5. The summation over K' in the first case however is over
   //  all integers up to 2l, so that in the case of even K', the first term (with the A+B) is zero and in the case of odd
   //  K', the second term (with just the B term) is zero.
   //
   //  The matrices a(K,K') and c(K,K') are solely functions of |i> and |j> and do not depend on the bessel functions
   //  <j_k>, or spherical harmonics Y^K_Q, and are calculated in the functions lovesey_aKK() and loveset_cKK() above.

   for(Kp=0; Kp<=(2*l+1); Kp++)                                 // K may take even values between [0,2l]; K' odd values [1,(2l+1)]
   {
      if(Kp%2==1)                                               // Calculates the matrices for the first term with A+B
      {
         iA = lovesey_aKK(Akk,Kp-1,Kp,n,l); if(iA) Akk*= (J[Kp-1]+J[Kp+1])*((2*Kp+1)/(Kp+1.)); else { std::cerr << errormsg; return; }
         ic = lovesey_cKK(ckk,Kp-1,Kp,n,l); if(ic) Bkk = ckk*J[Kp-1];                          else { std::cerr << errormsg; return; }
         ic = lovesey_cKK(ckk,Kp+1,Kp,n,l); if(ic) Bkk+= ckk*(J[Kp+1]*sqrt(Kp/(Kp+1.)));       else { std::cerr << errormsg; return; }
         Bkk *= ((Kp+1.)/(2*Kp+1.)) * ( (2*Kp+1)/(Kp+1.) ); K = Kp-1; ckk = Bkk+Akk;
      }
      else                                                      // Calculates the matrices for the second term with just B
      {
         ic = lovesey_cKK(ckk,Kp,Kp,n,l); if(ic) ckk *= J[Kp]; else { std::cerr << errormsg; return; } K = Kp;
      }

      for(Q=-K; Q<=K; Q++)
      {
         Qp = -(Q-q); if(Qp<-Kp || Qp>Kp) continue;             // 3j symbol requires Q+Q'-q = 0
 
         Ykq = spherical_harmonics(K,Q,theta,phi); 
         Tj = pow(-1.,K-Kp+q) * sqrt(3) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q);
         if((fabs(Ykq.r)<DBL_EPSILON && fabs(Ykq.i)<DBL_EPSILON) || fabs(Tj)<DBL_EPSILON) continue;

         Qmat.zero(ns,ns); QLmat.zero(ns,ns); QSmat.zero(ns,ns);

         nz = ckk.find();
         for(i=0; i<(int)nz.size(); i++)                        // The matrices above already contain the dependence on |vULSJ>
         {                                                      //    we now use the W-E theorem to add the Jz, Q, dependence
            J2 = conf.states[nz[i][0]].J2;
            J2p = conf.states[nz[i][1]].J2;
            valJ_i = (J2-minJ2)/2; valJ_j = (J2p-minJ2)/2;
            if(mJmat[valJ_i][valJ_j].isempty())
            {
               mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
               for(imJ_i=0; imJ_i<=J2; imJ_i++)
               {
                  Jz2 = imJ_i*2-J2;
                  for(imJ_j=0; imJ_j<=J2p; imJ_j++)
                  {
                     Jz2p = imJ_j*2-J2p; Tj2 = threej(2*Kp,J2p,J2,2*Qp,Jz2p,-Jz2);
                     if(fabs(Tj2)>DBL_EPSILON) mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,Kp+(-J2p+Jz2)/2.) * sqrt(J2+1.) * Tj2;
                  }
               }
            }
            Qmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*ckk(nz[i][0],nz[i][1]));
            if(Kp%2==1)
            {
               QLmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Akk(nz[i][0],nz[i][1]));
               QSmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*Bkk(nz[i][0],nz[i][1]));
            }
         }
         for(i=0; i<=valJ; i++) for(j=0; j<=valJ; j++) mJmat[i][j].clear();

         Qq[0] += Qmat*(Ykq.r*Tj); Qq[1] += Qmat*(Ykq.i*Tj);
         if(Kp%2==1) { Qq[2] += QSmat*(Ykq.r*Tj); Qq[3] += QSmat*(Ykq.i*Tj);  Qq[4] += QLmat*(Ykq.r*Tj); Qq[5] += QLmat*(Ykq.i*Tj); }
         else        { Qq[2] += Qmat*(Ykq.r*Tj); Qq[3] += Qmat*(Ykq.i*Tj); }
      }
   }

   for(i=0; i<6; i++) { rmzeros(Qq[i]); Qq[i] *= sqrt(4*PI); }
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the coefficient of the spin density operator M^S_q(r), after Balcar, J.Phys.C. v8, p1581, 1975, eqn 10.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> balcar_MSq(int q, int K, int Q, int n, orbital l)  
{
   std::string errormsg("balcar_MSq(): Unable to calculate the c(K,K') matrix\n");
   int i,j,ns=getdim(n,l);                                      // Number of states = ^{4l+2}C_{n}
   sMat<double> ckk,Mmat,Qmat; Mmat.zero(ns,ns);

   // Eqn 10 of Balcar 1975:
   //
   //         S                     -2uB     ---  K     2    ---   [     1/2+M+q-J'+L'+S'     1/2
   // <vSLJM|M (r)|v'S'L'J'M'> =  ---------  >   Y (r) U (r) >     [ (-1)                (3/2)
   //         q                   sqrt(4PI)  ---  Q          ---   [
   //                                        K,Q             K',Q'
   //
   //              \/                                1/2  ( l K l )  { 1  K  K' }
   //              /\    [l,l,S,S',L,L',J,J',K,K',K']     ( 0 0 0 )  { S' L' J' }
   //                                                                { S  L  J  }
   //                                            _ _
   //                     ---     _   _          S+L                                                         ]
   //              \/   n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }    (  J K' J' ) ( K K'  1 )  ]
   //              /\     ---                        { s  Sb s  }  { l  Lb l  }    ( -M Q' M' ) ( Q Q' -q )  ]
   //                      t
   //
   // In this function we calculate only the sum over the square brackets, leaving the K,Q sum to chrgplt.
   // The first 3j symbol imposes the selection rule that K=even,K<=2l. The 9j symbol means (K-1)<=K'<=(K+1), whilst
   // the final 3j symbol means that Q+Q'=q, that is Q'=q-Q, and that -K'<=Q'<=K', so there is only ever one Q' term
   // in the sum.

   if(K%2!=0 || K>(2*l) || Q<-K || Q>K || q<-1 || q>1) return Mmat;
   int Kp,Qp=q-Q;

   int k,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int imJ_i,imJ_j,J2,J2p,Jz2,Jz2p;
   std::vector<int> iJst,iJen;
   fconf conf(n,0,l); minJ2=99; maxJ2=0; k=0;
   for(i=0; i<(int)conf.states.size(); i++) {
      j = conf.states[i].J2; if(j<minJ2) minJ2 = j; if(j>maxJ2) maxJ2 = j; iJst.push_back(k+1); k+=j+1; iJen.push_back(k); }

   bool ic; 
   double Tj,Tj2,phase = 1.; /*if(q%2==1) phase = -1.;*/ if(K%4==2) phase = -phase;
   std::vector< std::vector<int> > nz;

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   sMat<double> mJmat_i(1,1);
   std::vector< std::vector< sMat<double> > > mJmat;
   std::vector< sMat<double> > mJmat_row;
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // Simplifying, using the C(K,K') coeficients previously calculated in lovesey_cKK()  (Balcar+Lovesey Book, eqn 3.8.5)
   //
   //          S                    -2uB    ---  K     2    ---   [             M+q-J'      1/2  -K                           ]
   //  <vSLJM|M (r)|v'S'L'J'M'> = --------- >   Y (r) U (r) >     [ C(K,K') (-1)      (3[J])    i   (  J K' J' ) ( K K'  1 )  ]
   //          q                  sqrt(4PI) ---  Q          ---   [                                 ( -M Q' M' ) ( Q Q' -q )  ]
   //                                       K,Q             K',Q'
   //

   // Eqn. 3.6.11
   //
   //                  K      1/2                            1/2      1/2+S'+L'
   // C(K,K') = <j >  i  (1/2)     [l,l,S,S',L,L',J',K,K',K']     (-1)
   //             K                                                       _ _
   //                              { 1  K  K' }    ---     _   _          S+L
   //              \/   ( l K l )  { S' L' J' }  n >   (t{|t) (t|}t') (-1)    { S  1  S' }  { L  K  L' }
   //              /\   ( 0 0 0 )  { S  L  J  }    ---                        { s  Sb s  }  { l  Lb l  }
   //                                               t

   for(Kp=(K-1); Kp<=(K+1); Kp++)
   {
      if(Qp<-Kp || Qp>Kp) continue; Qmat.zero(ns,ns);
      Tj = -phase * sqrt(3/PI) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q);
      ic = lovesey_cKK(ckk,K,Kp,n,l); if(ic) ckk *= Tj; else { std::cerr << errormsg; return Qmat; }
      
      nz = ckk.find();
      for(i=0; i<(int)nz.size(); i++)                        // The matrices above already contain the dependence on |vULSJ>
      {                                                      //    we now use the W-E theorem to add the Jz, Q, dependence
         J2 = conf.states[nz[i][0]].J2;
         J2p = conf.states[nz[i][1]].J2;
         valJ_i = (J2-minJ2)/2; valJ_j = (J2p-minJ2)/2;
         if(mJmat[valJ_i][valJ_j].isempty())
         {
            mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
            for(imJ_i=0; imJ_i<=J2; imJ_i++)
            {
               Jz2 = imJ_i*2-J2;
               for(imJ_j=0; imJ_j<=J2p; imJ_j++)
               {
                  Jz2p = imJ_j*2-J2p; Tj2 = threej(J2,2*Kp,J2p,-Jz2,2*Qp,Jz2p);
                  if(fabs(Tj2)>DBL_EPSILON) mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,(-J2p+Jz2)/2.+q) * sqrt(J2+1.) * Tj2;
               }
            }
         }
         Qmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*ckk(nz[i][0],nz[i][1]));
      }
      for(i=0; i<=valJ; i++) for(j=0; j<=valJ; j++) mJmat[i][j].clear();

      Mmat += Qmat;
   }
   return Mmat;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the coefficient of the orbital magnetic density operator M^L_q(r), after Balcar 1975, eqn 12.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> balcar_MLq(int q, int K, int Q, int n, orbital l)  
{
   std::string errormsg("balcar_MLq(): Unable to calculate the c(K,K') matrix\n");
   int i,j,ns=getdim(n,l);                                      // Number of states = ^{4l+2}C_{n}
   sMat<double> akk,Mmat,Qmat; Mmat.zero(ns,ns);

   // Eqn 12 of Balcar 1975:
   //                                                       infty
   //         L                     -2uB     ---  K     1  /     2     ---   [     q+M+L+L'+S
   // <vSLJM|M (r)|v'S'L'J'M'> =  ---------  >   Y (r) --- | dX U (X)  >     [ (-1)            d
   //         q                   sqrt(4PI)  ---  Q     r  /           ---   [                  SS'
   //                                        K,Q            r          K',Q'
   //
   //              \/                             1/2         1/2  ( l K l )  { l  K' l  }  { K' L' L  }
   //              /\    [l,l,l,L,L',J,J',K,K',K']    [l(l+1)]     ( 0 0 0 )  { K  l  1  }  { S  J  J' }
   //                                                                                      
   //                                            _
   //                     ---     _   _          L                                           ]
   //              \/   n >   (t{|t) (t|}t') (-1)  { L  K  L' }    (  J K' J' ) ( K K'  1 )  ]
   //              /\     ---                      { l  Lb l  }    ( -M Q' M' ) ( Q Q' -q )  ]
   //                      t
   //
   // In this function we calculate only the sum over the square brackets and include the factor -2/sqrt(4pi), leaving the K,Q sum to chrgplt.
   // The first 3j symbol imposes the selection rule that K=even,0<=K<=2l. The 9j symbol means (K-1)<=K'<=(K+1), whilst
   // the final 3j symbol means that Q+Q'=q, that is Q'=q-Q, and that -K'<=Q'<=K', so there is only ever one Q' term
   // in the sum.

   if(K%2!=0 || K>(2*l) || Q<-K || Q>K || q<-1 || q>1) return Mmat;
   int Kp,Qp=q-Q;

   int k,minJ2,maxJ2,valJ,valJ_i,valJ_j;
   int imJ_i,imJ_j,J2,J2p,Jz2,Jz2p;
   std::vector<int> iJst,iJen;
   fconf conf(n,0,l); minJ2=99; maxJ2=0; k=0;
   for(i=0; i<(int)conf.states.size(); i++) {
      j = conf.states[i].J2; if(j<minJ2) minJ2 = j; if(j>maxJ2) maxJ2 = j; iJst.push_back(k+1); k+=j+1; iJen.push_back(k); }

   bool iA;
   double Tj,Tj2,phase = 1.;/* if(q%2==1) phase = -1.;*/ if(K%4==2) phase = -phase;
   std::vector< std::vector<int> > nz;

   // Initialises a cell array of matrices of q- and Jz- dependent matrices so that we don't have to calculate
   //    each (2J+1)x(2J'+1) matrix more than once, for each J and J' values.
   sMat<double> mJmat_i(1,1);
   std::vector< sMat<double> > mJmat_row;
   std::vector< std::vector< sMat<double> > > mJmat;
   valJ = (maxJ2-minJ2)/2;
   for(i=0; i<=valJ; i++)
      mJmat_row.push_back(mJmat_i);
   for(i=0; i<=valJ; i++)
      mJmat.push_back(mJmat_row);

   // Noting that the e(K,K') coefficient is, from Balcar and Lovesey, eqn 3.6.4:
   //
   //                    i^K                            1/2          1/2      S+L+L'+J'  ( l K l ) { 1 K' K } { L  K' L' }
   //  e(K,K') = d     --------  [l,l,l,K,K',K',J',L,L']     [l(l+1)]     (-1)           ( 0 0 0 ) { l l  l } { J' S  J  }
   //             SS'  2 sqrt(3)
   //                                            _
   //                     ---     _   _          L
   //              \/   n >   (t{|t) (t|}t') (-1)  { L  K  L' }
   //              /\     ---                      { l  Lb l  }
   //                      t
   //
   // We can simplify the matrix element above, eqn 3.8.9. Note the 2/r rather than 1/r in the previous equation.
   //
   //          L                    -2uB    ---  K      infty
   //  <vSLJM|M (r)|v'S'L'J'M'> = --------- >   Y (r) 2  /     2    ---   [          -K       1/2    M-J'+q                           ]
   //          q                  sqrt(4PI) ---  Q   --- | dX U(X)  >     [ e(K,K') i  (3[J])    (-1)       (  J K' J' ) ( K  K'  1 ) ]
   //                                       K,Q       r  /          ---   [                                 ( -M Q' M' ) ( Q  Q' -q ) ]
   //                                                    r          K',Q'
   //
   // Finally, we use the relation 4.2.14 of Balcar and Lovesey,  e(K,K') = +/- 0.5(2K'+1) a(K,K')  where + is for K'=K+1, - for K'=K-1

   for(Kp=(K-1); Kp<=(K+1); Kp+=2)       // The 3j symbol in A(K,K') means that K==Kp gives zero... (NB. Does this apply to E(K,K') too?)
   {
      if(Qp<-Kp || Qp>Kp) continue; Qmat.zero(ns,ns);
      if(Kp==(K-1)) Tj = -phase * sqrt(3/PI) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q) * -(2*Kp+1);
      else          Tj = -phase * sqrt(3/PI) * threej(2*K,2*Kp,2,2*Q,2*Qp,-2*q) *  (2*Kp+1);
      iA = lovesey_aKK(akk,K,Kp,n,l); if(iA) akk *= Tj; else { std::cerr << errormsg; return Qmat; }
      
      nz = akk.find();
      for(i=0; i<(int)nz.size(); i++)                        // The matrices above already contain the dependence on |vULSJ>
      {                                                      //    we now use the W-E theorem to add the Jz, Q, dependence
         J2 = conf.states[nz[i][0]].J2;
         J2p = conf.states[nz[i][1]].J2;
         valJ_i = (J2-minJ2)/2; valJ_j = (J2p-minJ2)/2;
         if(mJmat[valJ_i][valJ_j].isempty())
         {
            mJmat[valJ_i][valJ_j].zero(J2+1,J2p+1);
            for(imJ_i=0; imJ_i<=J2; imJ_i++)
            {
               Jz2 = imJ_i*2-J2;
               for(imJ_j=0; imJ_j<=J2p; imJ_j++)
               {
                  Jz2p = imJ_j*2-J2p; Tj2 = threej(J2,2*Kp,J2p,-Jz2,2*Qp,Jz2p);
                  if(fabs(Tj2)>DBL_EPSILON) mJmat[valJ_i][valJ_j](imJ_i,imJ_j) = pow(-1.,(-J2p+Jz2)/2.+q) * sqrt(J2+1.) * Tj2;
               }
            }
         }
         Qmat.pset(iJst[nz[i][0]],iJen[nz[i][0]],iJst[nz[i][1]],iJen[nz[i][1]],mJmat[valJ_i][valJ_j]*akk(nz[i][0],nz[i][1]));
      }
      for(i=0; i<=valJ; i++) for(j=0; j<=valJ; j++) mJmat[i][j].clear();

      Mmat += Qmat;
   }
   return Mmat;
}

complexdouble * balcar_Mq(int xyz, int K, int Q, int n, orbital l)
{
#define NSTR(K,Q) nstr[3] = K+48; nstr[4] = Q+48; nstr[5] = 0
#define MSTR(K,Q) nstr[3] = K+48; nstr[4] = 109;  nstr[5] = Q+48; nstr[6] = 0
#define NPOS std::string::npos
   int Hsz = getdim(n,l);
   sMat<double> retval_r(Hsz,Hsz),retval_i(Hsz,Hsz),qpp,qmp,qpm,qmm;
   char nstr[7]; char filename[255]; char basename[255]; strcpy(basename,"results/mms/");
   nstr[0] = (l==F?102:100); if(n<10) { nstr[1] = n+48; nstr[2] = 0; } else { nstr[1] = 49; nstr[2] = n+38; nstr[3] = 0; }
   strcat(basename,nstr); strcat(basename,"_"); nstr[0] = 77;   // ASCII codes: 77="M", 83=="S", 100=="d", 102=="f", 109=="m", 112="p"
   if(xyz==1||xyz==2)
   {
      nstr[1]=83;
      if(Q!=0)
      {
         NSTR(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MSq(1,K,abs(Q),n,l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTR(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MSq(1,K,-abs(Q),n,l); rmzeros(qmp); mm_gout(qmp,filename); }
         NSTR(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MSq(-1,K,abs(Q),n,l); rmzeros(qpm); mm_gout(qpm,filename); }
         MSTR(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmm = mm_gin(filename); if(qmm.isempty()) { qmm = balcar_MSq(-1,K,-abs(Q),n,l); rmzeros(qmm); mm_gout(qmm,filename); }
	 // sum to coeff of Zlm (neglecting a 1/sqrt(2) factor) 
         if(Q<0) { if(Q%2==0) { qmp -= qpp; qmm -= qpm; } else { qmp += qpp; qmm += qpm; } }
         else    { if(Q%2==0) { qmp += qpp; qmm += qpm; } else { qmp -= qpp; qmm -= qpm; } }
         // add spherical components and multiply by addition factor 1/sqrt(2)which was neglected in the line above
         if(xyz==1) { if(Q<0) retval_i = (qmm-qmp)/(-2.);    else retval_r = (qmm-qmp)/2.; }
         if(xyz==2) { if(Q<0) retval_r = (qmm+qmp)/2.; else retval_i = (qmm+qmp)/2.; }// changed MR 25.5.2010 // Q<0 signs changed MR 28.5.2010
      }
      else
      {
         NSTR(K,0); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MSq(1,K,0,n,l); rmzeros(qpp); mm_gout(qpp,filename); }
         NSTR(K,0); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MSq(-1,K,0,n,l); rmzeros(qpm); mm_gout(qpm,filename); }      
         if(xyz==1) {  retval_r = (qpp-qpm)/(-sqrt(2.)); }
         if(xyz==2) {  retval_i = (qpp+qpm)/sqrt(2.); }// changed MR 25.5.2010
      }
   }
   else if(xyz==3)
   {
      nstr[1]=83;
      if(Q!=0)
      {
         NSTR(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MSq(0,K,abs(Q),n,l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTR(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MSq(0,K,-abs(Q),n,l); rmzeros(qmp); mm_gout(qmp,filename); }
         if(Q<0) {  if(Q%2==0) qmp -= qpp; else qmp += qpp;  retval_i = qmp/(-sqrt(2.));}
         else    {  if(Q%2==0) qmp += qpp; else qmp -= qpp;  retval_r = qmp/sqrt(2.);}// changed by MR 28.5.2010
      }
      else
      {
         NSTR(K,0); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         retval_r = mm_gin(filename); if(retval_r.isempty()) { retval_r = balcar_MSq(0,K,0,n,l); rmzeros(retval_r); mm_gout(retval_r,filename); }
      }
   }
   else if(xyz==-1||xyz==-2)
   {
      nstr[1]=76;
      if(Q!=0)
      {
         NSTR(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MLq(1,K,abs(Q),n,l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTR(K,abs(Q)); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MLq(1,K,-abs(Q),n,l); rmzeros(qmp); mm_gout(qmp,filename); }
         NSTR(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MLq(-1,K,abs(Q),n,l); rmzeros(qpm); mm_gout(qpm,filename); }
         MSTR(K,abs(Q)); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmm = mm_gin(filename); if(qmm.isempty()) { qmm = balcar_MLq(-1,K,-abs(Q),n,l); rmzeros(qmm); mm_gout(qmm,filename); }
	 // sum to coeff of Zlm (neglecting a 1/sqrt(2) factor) 
         if(Q<0) { if(Q%2==0) { qmp -= qpp; qmm -= qpm; } else { qmp += qpp; qmm += qpm; } }
         else    { if(Q%2==0) { qmp += qpp; qmm += qpm; } else { qmp -= qpp; qmm -= qpm; } }
         // add spherical components and multiply by addition factor 1/sqrt(2)which was neglected in the line above
         if(xyz==-1) { if(Q<0) retval_i = (qmm-qmp)/(-2.);    else retval_r = (qmm-qmp)/2.; }
         if(xyz==-2) { if(Q<0) retval_r = (qmm+qmp)/2.; else retval_i = (qmm+qmp)/2.; }// changed MR 25.5.2010 // Q<0 signs changed MR 28.5.2010
      }
      else
      {
         NSTR(K,0); nstr[2]=112; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MLq(1,K,0,n,l); rmzeros(qpp); mm_gout(qpp,filename); }
         NSTR(K,0); nstr[2]=109; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpm = mm_gin(filename); if(qpm.isempty()) { qpm = balcar_MLq(-1,K,0,n,l); rmzeros(qpm); mm_gout(qpm,filename); }
         if(xyz==-1) {  retval_r = (qpp-qpm)/(-sqrt(2.)); }
         if(xyz==-2) {  retval_i = (qpp+qpm)/sqrt(2.); }// changed MR 25.5.2010
      }
   }
   else if(xyz==-3)
   {
      nstr[1]=76;
      if(Q!=0)
      {
         NSTR(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qpp = mm_gin(filename); if(qpp.isempty()) { qpp = balcar_MLq(0,K,abs(Q),n,l); rmzeros(qpp); mm_gout(qpp,filename); }
         MSTR(K,abs(Q)); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         qmp = mm_gin(filename); if(qmp.isempty()) { qmp = balcar_MLq(0,K,-abs(Q),n,l); rmzeros(qmp); mm_gout(qmp,filename); }
         if(Q<0) {  if(Q%2==0) qmp -= qpp; else qmp += qpp; retval_i = qmp/(-sqrt(2.));}
         else    {  if(Q%2==0) qmp += qpp; else qmp -= qpp; retval_r = qmp/sqrt(2.);}// changed by MR 28.5.2010
      }
      else
      {
         NSTR(K,0); nstr[2]=48; strcpy(filename,basename); strcat(filename,nstr); strcat(filename,".mm");
         retval_r = mm_gin(filename); if(retval_r.isempty()) { retval_r = balcar_MLq(0,K,0,n,l); rmzeros(retval_r); mm_gout(retval_r,filename); }
      }
   }
   return zmat2f(retval_r,retval_i);
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile:
//    g++ -g states.cpp cfp.cpp njsyms.cpp so_cf.cpp maths.cpp mmio.cpp lovesey.cpp -llapack -lblas; ./a.out
/* --------------------------------------------------------------------------------------------------------------- //
int main(int argc, char *argv[])
{
   int n,K,Kp;
   orbital l=F;
   if(argc>1) n = atoi(argv[1]); else n = 2;
   if(argc>2) K = atoi(argv[2]); else K = 2;
   if(argc>3) Kp = atoi(argv[3]); else Kp = 1;

   clock_t start,end; start = clock();

   sMat<double> aKK,cKK;

   bool ia = lovesey_aKK(aKK,K,Kp,n,l);
   end = clock(); std::cerr << "Time to calculate a(K,K') = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "a(" << K << "," << Kp << ") is "; if(ia) std::cerr << "real\n"; else std::cerr << "imaginary\n";
   std::cout << aKK << "\naKK=x; clear x;\n";

   start = clock();
   bool ic = lovesey_cKK(cKK,K,Kp,n,l);
   end = clock(); std::cerr << "Time to calculate c(K,K') = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   std::cerr << "c(" << K << "," << Kp << ") is "; if(ia) std::cerr << "real\n"; else std::cerr << "imaginary\n";
   std::cout << cKK << "\ncKK=x; clear x;\n";

   start = clock();
   std::vector<double> Jvec(6,1.);
   lovesey_Qq(aKK,cKK,0,n,l,Jvec);
   end = clock(); std::cerr << "Time to calculate Q = " << (double)(end-start)/CLOCKS_PER_SEC << "s.\n";
   std::cout << aKK << "Qq=x; clear x;\n" << cKK << "iQq=x; clear x;\n";

   return 0;
}*/
