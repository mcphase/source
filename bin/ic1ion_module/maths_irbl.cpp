/* maths_irbl.cpp
 *
 * This file is a C++ implementation of the Implicitly Restarted Block Lanczos method to solve the sparse Hermitian
 * eigenvalue problem. The original work was by James Baglama, Daniela Calvetti and Lothar Reichel, and was implemented
 * as a MATLAB program: irbleigs.m
 *
 * Reference: IRBL: An Implicitly Restarted Block Lanczos Method for large-scale Hermitian eigenproblems,
 *            J. Baglama, D. Calvetti, and L. Reichel, SIAM J. Sci. Comput. v24, No. 5, (2003), pp. 1650-1677 
 *            http://www.math.uri.edu/~jbaglama/
 *
 * Class:					 		//
 * 								//
 * Function:							//
 *
 * This file is part of the ic1ionmodule of the McPhase package, calculating the single-ion properties of a rare
 * earth or actinide ion in intermediate coupling.
 *
 * (c) 2008 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 *
 */

#include "maths.hpp"

// --------------------------------------------------------------------------------------------------------------- //
// Computes the  QR factorization for a singular off diagonal block of the block tridiagonal matrix T
//    !!Untested!! - H_IC matrices don't seem to have singular blocks - at least that I could find easily...
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> qrsingblk(sMat<double> & F,     // 
                       sMat<double> & V,     // 
		       sMat<double> & eigvec,// 
		       sMat<int> & I,        // 
		       int & blsz,           // Block size
		       double & sqrteps,     // sqrt(eps)
		       sMat<double> & R)     // The R matrix of Q*R
{
   sMat<double> Fp = F;                      //
   sMat<double> Mt(blsz,blsz);               // Temporary Matrix
   std::vector<int> Vt;                      // Temporary Vector
   sMat<double> dotVF,doteF,v1,v2;
   double dotFF;
   sMat<int> It;
   int n = F.nr();
   int i,j,k,kk;
   bool allflag;

   // Initialize off-diagonal block to zero
   R = Mt;
   
   for(k=0; k<blsz; k++)
   {
      allflag = true; It = I-k; Vt = all(It); for(kk=0; kk<(int)Vt.size(); kk++) if(Vt[kk]==0) allflag = false;
      if(!allflag)  // if ~all(I-k)
      {
         // Set F(:,k) to a random vector and orthogonalizes F(:,k) to all previous Lanczos vectors and any converged eigenvectors.
         F.pset(-1,-1,k,k,  randn(n,1));
         Mt = V; Mt.transpose();
         dotVF = Mt*F.setp(-1,-1,k,k); F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k) - V*dotVF);
         if(k>1) 
         {
            Mt = F.setp(-1,-1,1,k-1); Mt.transpose();
            dotVF = Mt*F.setp(-1,-1,k,k); F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k) - F.setp(-1,-1,1,k-1)*dotVF);
         }
         // Iterative refinement.
         Mt = V; Mt.transpose(); dotVF = Mt*F.setp(-1,-1,k,k); F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k) - V*dotVF);
         if(k>1) 
         {
            Mt = F.setp(-1,-1,1,k-1); Mt.transpose();
            dotVF = Mt*F.setp(-1,-1,k,k); F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k) - F.setp(-1,-1,1,k-1)*dotVF);
         }
         // Orthogonalize F(:,k) against all converged eigenvectors.
         if(!eigvec.isempty())
         {
	    Mt = F.setp(-1,-1,k,k); Mt.transpose(); Mt *= eigvec; Mt.transpose();
            F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k) - eigvec*Mt); Mt = F.setp(-1,-1,k,k); Mt.transpose(); doteF = Mt*eigvec; doteF.transpose();
            if(norm(doteF) > sqrteps) F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k) - eigvec*doteF);
         }
         Mt = F.setp(-1,-1,k,k); F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k)/norm(Mt));
         // Set diagonal element to zero.
         R.del(k,k);
      }
      else
      {
         Mt = F.setp(-1,-1,k,k); R(k,k) = norm(Mt);
         F.pset(-1,-1,k,k,   F.setp(-1,-1,k,k)/R(k,k));
      }
      // Modified Gram-Schmidt orthogonalization.
      for(j=k+1; j<blsz; j++)
      {
         //Mt = F.setp(-1,-1,k,k); Mt.transpose();
         //R(k,j) = Mt*F.setp(-1,-1,j,j);
	 dotFF = 0.; v1 = F.setp(-1,-1,k,k); v2 = F.setp(-1,-1,j,j);
	 for(i=0; i<F.nr(); i++)
            dotFF += v1(i,0)*v2(i,0);
	 R(k,j) = dotFF; Mt = F.setp(-1,-1,k,k);
         F.pset(-1,-1,j,j,   F.setp(-1,-1,j,j) - Mt*R(k,j));
      } 
      // Iterative refinement
      for(j=k+1; j<blsz; j++)
      {
         //Mt = F.setp(-1,-1,k,k); Mt.transpose();
         //dotFF = Mt*F.setp(-1,-1,j,j);
	 dotFF = 0.; v1 = F.setp(-1,-1,k,k); v2 = F.setp(-1,-1,j,j);
	 for(i=0; i<F.nr(); i++)
            dotFF += v1(i,0)*v2(i,0);
         R(k,j) = R(k,j) + dotFF; Mt = F.setp(-1,-1,k,k);
         F.pset(-1,-1,j,j,   F.setp(-1,-1,j,j) - Mt*dotFF);
      } 
   } 

   Mt = F; F = Fp;
   return Mt;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the Block Lanczos decomposition: A*V = V*T + F*E^T
// --------------------------------------------------------------------------------------------------------------- //
bool blanz(sMat<double> & A,                 // Matrix to be decomposed
           int K,                            // Number of eigenvalues desired
	   sMat<double> & V,                 // Matrix of starting eigenvectors
	   int & blsz,                       // Block size
	   sMat<double> & eigvec,            // Matrix of converged eigenvectors
	   int & mprod,                      // Number of matrix products
	   int n,                            // Order of A
	   int nbls,                         // Number blocks in decomposition
	   std::vector<int> & singblk,       // ??
	   double sqrteps,                   // sqrt(eps)
	   double tol,                       // tolerance
	   sMat<double> & F,                 // Matrix of decomposition
	   sMat<double> & T)                 // Matrix of decomposition
{
   F.zero(n,blsz); T.clear(); int J = 1;     // Main loop count
   singblk.clear(); singblk.push_back(0);
   bool retval = true;
   int Jblsz,Jm1blszp1,i;

   sMat<double> Mt,Mt2,doteV,D,dotFV,dotFV2; // Temporary matrices
   sMat<double> doteF;
   double stol,dtmp;
   std::vector<double> Vt;
   bool Iflag = false;
   sMat<int> I;
   int sizevec,Isz = 0;

   // Checks size of V and eigvec
   if(V.nr()!=n || blsz < 1)
   {
      std::cerr << "Incorrect size of starting vector V in blanz.";
      retval = false;
   }
   if(V.nc()!=blsz) blsz = V.nc();
   if(!eigvec.isempty())
      if(eigvec.nr()!=n)
      {
         std::cerr << "Incorrect size of EIGVEC in blanz.";
         retval = false;
      }

   if(retval)
   {
      // First orthogonalisation of starting vectors.
      V = orth(V);
    //std::cout << "V_orth = " << V.display_full(); 

      // Orthogonalize V against all converged eigenvectors.
      if(!eigvec.isempty())
      {
         if(eigvec.nc() < V.nc())
         {
            Mt = eigvec; Mt.transpose();
            V = V - eigvec*(Mt*V); doteV = Mt*V;
         }
         else
         {
            Mt = V; Mt.transpose(); Mt *= eigvec; Mt.transpose();
            V = V - eigvec*Mt; 
            Mt = V; Mt.transpose(); doteV = Mt*eigvec; doteV.transpose();
         }
         if(norm(doteV) > sqrteps) V = V - eigvec*doteV;
      }
   
      // First check of linear dependence of starting vector(s). If starting vectors(s) are linearly dependent
      //    then add normalised random vectors and reorthogonalised.
      if(rank(V,sqrteps*n*blsz) < blsz)
      {
         V = V + randn(n,blsz); V = orth(V);
         if(!eigvec.isempty())
         {
            if(eigvec.nc() < V.nc())
            {
               Mt = eigvec; Mt.transpose();
	       V = V - eigvec*(Mt*V); doteV = Mt*V;
	    }
	    else
            {
	       Mt = V; Mt.transpose(); Mt *= eigvec; Mt.transpose();
               V = V - eigvec*Mt; 
               Mt = V; Mt.transpose(); doteV = Mt*eigvec; doteV.transpose();
            }
            if(norm(doteV) > sqrteps) V = V - eigvec*doteV;
         }
      }

      // If needed second orthogonalization of starting vectors.
      if(!eigvec.isempty()) V = orth(V);

      // Second check of linear dependence of starting vector(s). If starting vector(s) are linearly dependent then add
      //    normalized random vectors and reorthogonalize.
      if(V.nc()<blsz)
      {
         V.pset(1,A.nr(),A.nc()+1,blsz,randn(n,blsz));   // V = [V,randn(n,blsz-size(V,2))];
         if(!eigvec.isempty())
         {
            if(eigvec.nc() < V.nc())
            {
	       Mt = eigvec; Mt.transpose();
               V = V - eigvec*(Mt*V); doteV = Mt*V;
	    }
            else
	    {
	       Mt = V; Mt.transpose(); Mt *= eigvec; Mt.transpose();
               V = V - eigvec*Mt; 
	       Mt = V; Mt.transpose(); doteV = Mt*eigvec; doteV.transpose();
            }
	    if(norm(doteV) > sqrteps) V = V - eigvec*doteV;
         }
         // If needed third orthogonalization of starting vectors.
         V = orth(V);
         // Third check of linear dependence of starting vector(s). If starting vector(s) are linearly dependent fatal error return.
         if(V.nc()<blsz)
         {
            std::cerr << "Dependent starting vectors in block Lanczos\n";
	    retval = false;
         }
      }
   } // if(retval)
   if(retval)
   {
      if(blsz*nbls >= n) nbls = (int)floor((double)n/(double)blsz);
      while(J<=nbls)
      {
         // Values used for indices.
         Jblsz = J*blsz; Jm1blszp1 = blsz*(J-1)+1;
	 // Standard eigenvalue problem.
         F = A*V.setp(-1,-1,Jm1blszp1,Jblsz);  // F = A*V(:,Jm1blszp1:Jblsz);
	 // Count the number of matrix vector products.
	 mprod += blsz;

	 // Orthogonalize F against the previous Lanczos vectors.
	 if(J>1)
	 {
	    Mt = T.setp(Jm1blszp1,Jblsz,blsz*(J-2)+1,Jm1blszp1-1); Mt.transpose();
	    F = F - V.setp(-1,-1,blsz*(J-2)+1,Jm1blszp1-1) * Mt; 
	 }

	 // Compute the diagonal block of T.
	 Mt = F; Mt.transpose();
	 D = Mt * V.setp(-1,-1,Jm1blszp1,Jblsz);

	 // One step of the block classical Gram-Schmidt process.
	 F = F - V.setp(-1,-1,Jm1blszp1,Jblsz)*D;

	 // Full reorthogonalization step.
	 Mt = F; Mt.transpose(); dotFV = Mt*V; dotFV.transpose(); // dotFV = (F'*V)';
	 F = F - V*dotFV;
	 if(norm(dotFV)>sqrteps)
	 {
	    Mt = F; Mt.transpose(); dotFV2 = Mt*V; dotFV2.transpose();
	 }
	 for(i=1; i<=J; i++)
	    D = D + dotFV.setp(blsz*(i-1)+1,blsz*i,-1,-1);

	 // Orthogonalize F against all converged eigenvectors.
	 if(!eigvec.isempty())
	 {
	    if(eigvec.nc() < F.nc())
	    {
	       Mt = eigvec; Mt.transpose();
               F = F - eigvec*(Mt*F); doteF = Mt*F;
	    }
	    else
	    {
	       Mt = F; Mt.transpose(); Mt *= eigvec; Mt.transpose();
               F = F - eigvec*Mt; 
               Mt = F; Mt.transpose(); doteF = Mt*eigvec; doteF.transpose();
            }
         }

	 // To ensure a symmetric matrix T is computed.
	 Mt = tril(D); Mt.transpose();
         T.pset(Jm1blszp1,Jblsz,Jm1blszp1,Jblsz,/* = */ tril(D,-1) + Mt);

         // Compute QR factorization and off diagonal block of T.
         if(J<nbls)
         {
          //V.pset(-1,-1,Jblsz+1,blsz*(J+1),   qr(F,Mt,0)); T(Jblsz+1,blsz*(J+1),Jm1blszp1,Jblsz,   Mt);
	    T.pset(Jblsz+1,blsz*(J+1),Jm1blszp1,Jblsz,   qr(F,Mt,0));     // [V(:,Jblsz+1:blsz*(J+1)),
            V.pset(-1,-1,Jblsz+1,blsz*(J+1),             Mt);             //  T(Jblsz+1:blsz*(J+1),Jm1blszp1:Jblsz)] = qr(F,0);
            
            // Check for linearly dependent vectors among the blocks.
            Vt = diag(T.setp(Jblsz+1,blsz*(J+1),Jm1blszp1,Jblsz)); Vt = vabs(Vt); dtmp = DBL_EPSILON*vmax(Vt);
            stol = (sqrteps<tol) ? sqrteps : tol; stol = (stol>dtmp) ? stol : dtmp;
            for(i=0; i<(int)Vt.size(); i++)
               if(Vt[i] <= stol)
               {
                  I(i,i) = 1; Iflag = true; Isz++;
               }

            // Linearly dependent vectors detected in the current block.
            if(Iflag)
            {
	       // Exit. Convergence or not enough vectors to continue to build up the space.
               if(!eigvec.isempty()) sizevec = eigvec.nc(); else sizevec = 0; 
               if( (Isz==blsz && T.nc()>=K) || (sizevec+T.nc())>=n )
	       {
                  // Resize T and V and exit.
                  T.resize(1,T.nc(),1,T.nc()); V.resize(-1,-1,1,T.nc());
                  return retval;
               }
 
               // Full Reorthogonalization step to ensure orthogonal vectors.
               V.resize(-1,-1,1,Jblsz);
	       Mt = F; Mt.transpose(); dotFV = Mt*V; dotFV.transpose(); F = F - V*dotFV;
	       if(norm(dotFV)>sqrteps) 
               {
                  Mt = F; Mt.transpose(); dotFV2 = Mt*V; dotFV2.transpose();
                  dotFV = dotFV+dotFV2; F = F - V*dotFV2;
               }
               for(i=0; i<J; i++)
                  T.pset(Jm1blszp1,Jblsz,Jm1blszp1,Jblsz, T.setp(Jm1blszp1,Jblsz,Jm1blszp1,Jblsz) + dotFV.setp(blsz*(i-1)+1,blsz*i,-1,-1) );

               // Orthogonalize F against all converged eigenvectors.
	       if(!eigvec.isempty())
               {
                  if(eigvec.nc()<F.nc())
                  {
                     Mt = eigvec; Mt.transpose(); F = F - eigvec*(Mt*F); doteF = Mt*F;
                  }
                  else
                  {
                     Mt = eigvec; Mt.transpose(); Mt *= eigvec; Mt.transpose();
                     F = F - eigvec*Mt; Mt = F; Mt.transpose(); doteF = Mt*eigvec; doteF.transpose();
                  }
                  if(norm(doteF)>sqrteps) F = F - eigvec*doteF;
               }

               // To ensure a symmetric matrix T is computed.
               Mt = T.setp(Jm1blszp1,Jblsz,Jm1blszp1,Jblsz); Mt = tril(Mt); Mt.transpose();
               T.pset(Jm1blszp1,Jblsz,Jm1blszp1,Jblsz, tril(T.setp(Jm1blszp1,Jblsz,Jm1blszp1,Jblsz),-1) + Mt);

               // Re-compute QR with random vectors.
               Mt2 = V.setp(-1,-1,1,Jblsz);
               V.pset(-1,-1,Jblsz+1,blsz*(J+1), qrsingblk(F,Mt2,eigvec,I,blsz,sqrteps, Mt));
               T.pset(Jblsz+1,blsz*(J+1),Jm1blszp1,Jblsz, Mt);

               // Set the singular block indicator to true along with which vector(s) are linearly dependent.
               singblk[0] = 1; for(i=0; i<I.nr(); i++) singblk.push_back(Jblsz+I(i,i));
            }

            // Set off diagonal blocks to be equal.
            Mt = T.setp(Jblsz+1,blsz*(J+1),Jm1blszp1,Jblsz); Mt.transpose();
            T.pset(Jm1blszp1,Jblsz,Jblsz+1,blsz*(J+1),Mt);
         }
         J++;
      }

   } // if(retval)

   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Checks the convergence of Ritz values and Ritz vectors
// --------------------------------------------------------------------------------------------------------------- //
void convtests(bool computvec,int &conv, bool &deflate, std::vector<double> &eigval, sMat<double> &eigvec,
               std::vector<double> &eigresdisp, int iter, int K, std::vector<double> &pritz,
               std::vector<double> residuals, eigVE<double> &ritz, bool &ritzconv, std::vector<int> &singblk,
               double sqrteps, double tol, int Tsz, sMat<double> V, double &Rmax)
{
   unsigned int i;
   std::vector<int> I,Vt1,Vt2,Vt3;                    // I are index of Jrv indices. Vt are temporary vectors
   sMat<double> Mt,Mt2;                               // Temporary matrices
   double RVTol;                                      // Ritz value tolerance

   // Initialization of local variables.
   std::vector<int> dif;                              // Place holder for the difference of sets Jr(Jre) and Jrv(Jrev).
   std::vector<int> Jr; Jr.reserve(pritz.size());     // Place holder for which Ritz values converged.
   std::vector<int> Jre; Jre.reserve(pritz.size());   // Place holder for which desired Ritz values onverged.
   std::vector<int> Jrv; Jrv.reserve(pritz.size());   // Place holder for which Ritz vectors converged.
   std::vector<int> Jrev; Jrev.reserve(pritz.size()); // Place holder for which desired Ritz vectors converged.
   std::vector<int> ST; ST.reserve(pritz.size());     // Place holder for which Ritz values stagnated.

   // Compute maximum Ritz value to estimate ||A||_2.
   if(Rmax==-1.)                                                  // if isempty(Rmax)
   {  Rmax = fabs(ritz.E[Tsz-1]); }                               //    Rmax = abs(ritz(Tsz));
   else                                                           // else
   {  if(fabs(ritz.E[Tsz-1])>Rmax) Rmax = fabs(ritz.E[Tsz-1]); }  //    Rmax = max(Rmax,abs(ritz(Tsz)));
                                                                  // end
   if(pow(DBL_EPSILON,2./3.)>Rmax) Rmax = pow(DBL_EPSILON,2./3.); // Rmax = max(eps^(2/3),Rmax);

   // Compute tolerance to determine when a Ritz Vector has converged.
   RVTol = (sqrteps>tol) ? tol : sqrteps;                         // RVTol = min(sqrteps,tol);

   // Check for stagnation of Ritz values. 
   //    eps*100 is used for a tolerance to determine when the desired Ritz values are stagnating.
   if(iter > 1)
   {
      if(pritz.size()==ritz.E.size())                             // if length(pritz) == length(ritz)
         for(i=0; i<pritz.size(); i++)
            if(fabs(pritz[i]-ritz.E[i])<(DBL_EPSILON*100))        //    ST = find(abs(pritz-ritz) < eps*100);
	       ST.push_back(i);
   }                                                              // end

   // Check for convergence of Ritz vectors.
   for(i=0; i<residuals.size(); i++)
      if(residuals[i]<(RVTol*Rmax)) Jrv.push_back(i);             // Jrv  = find(residuals < RVTol*Rmax);
   for(i=0; i<Jrv.size(); i++)
      if(Jrv[i]<(K-conv)) Jrev.push_back(i);                      // Jrev = find(Jrv <= K-conv);
   
   // Check for convergence of Ritz values.
   for(i=0; i<residuals.size(); i++)
      if(residuals[i]<(tol*Rmax)) Jr.push_back(i);
   Jr = setunion(Jr,ST);                                          // Jr = union(find(residuals < tol*Rmax),ST);
   for(i=0; i<Jr.size(); i++)
      if(Jr[i]<(K-conv)) Jre.push_back(i);                        // Jre =  find(Jr  <= K-conv);

 //std::cerr << "Jr=" << dispvect(Jr) << "\n";
 //std::cerr << "Jre=" << dispvect(Jre) << "\n";
 //std::cerr << "Jrev=" << dispvect(Jrev) << "\n";
 //std::cerr << "Jrv=" << dispvect(Jrv) << "\n";
     
   // Remove common values in Jre and Jrev. Common values indicate a desired Ritz pair has converged and will be deflated.
   if(!Jre.empty())                                               // if ~isempty(Jr(Jre))
   {
      for(i=0; i<Jre.size(); i++)  Vt1.push_back(Jr[Jre[i]]); 
      for(i=0; i<Jrev.size(); i++) Vt2.push_back(Jrv[Jrev[i]]); 
    //std::cerr << "Vt1=" << dispvect(Vt1) << "\n";
    //std::cerr << "Vt2=" << dispvect(Vt2) << "\n";
      dif = setdiff(Vt1,Vt2);                                     //    dif = setdiff(Jr(Jre),Jrv(Jrev))
    //std::cerr << "setdiff(Vt1,Vt2)=" << dispvect(dif) << "\n";
   }                                                              // end

   // Determine the number of converged desired Ritz vectors.
   conv += (int)Jrev.size();                                      // conv = conv + length(Jrev);

 //std::cerr << "conv = " << conv << "\n";
 //std::cerr << "Rmax = " << Rmax << "\n";
   // Determine if the requested number of desired Ritz values have converged.
   if((conv+(int)dif.size()) >= K)                                // if conv+length(dif) >= K
   {
      Vt3 = setunion(Vt1,Vt2);
      for(i=0; i<Vt3.size(); i++)
         eigval.push_back(ritz.E[Vt3[i]]);                        //    eigval = [eigval;ritz(union(Jr(Jre),Jrv(Jrev)))];

      // If eigenvectors are requested then compute eigenvectors.
      if(computvec)                                               //    if strcmp(computvec,'T')
      {
         Mt.zero();
         for(i=0; i<Vt3.size(); i++)
            Mt.pset(-1,-1,i+1,i+1, ritz.V.setp(-1,-1,Vt3[i]+1,Vt3[i]+1));
         if(eigvec.isempty())                                     //       if isempty(eigvec)
            eigvec = V*Mt;                                        //          eigvec = V*ritzvec(:,union(Jr(Jre),Jrv(Jrev)));
         else                                                     //       else
            eigvec.pset(-1,-1,eigvec.nc()+1,-1, V*Mt);            //          eigvec = [eigvec,V*ritzvec(:,union(Jr(Jre),Jrv(Jrev)))];
      }                                                           //    end

      // Set convergence to true and exit.
      ritzconv = true; return;                                    //    ritzconv = 'T'; return;
   }                                                              // end

   // Store previous values of Ritz values to check for stagnation.
   for(i=0; i<ritz.E.size(); i++) I.push_back(i);                 // I = [1:length(ritz)];
   I = setxor(I,Jrv); pritz.clear();                              // I = setxor(I,Jrv);
   for(i=0; i<I.size(); i++) pritz.push_back(ritz.E[I[i]]);       // pritz = ritz(I);

   // Compute converged Ritz vectors so they can be deflated.    
   if(!Jrv.empty())                                               // if length(Jrv) > 0
   {
      for(i=0; i<Jrv.size(); i++)
      {  eigval.push_back(ritz.E[Jrv[i]]);                        //    eigval = [eigval;ritz(Jrv)];
         eigresdisp.push_back(residuals[Jrv[i]]); }               //    eigresdisp = [eigresdisp,residuals(Jrv)];
      Mt.zero(ritz.V.nr(),(int)Jrv.size());
      for(i=0; i<Jrv.size(); i++)
         Mt.pset(-1,-1,i+1,i+1, ritz.V.setp(-1,-1,Jrv[i]+1,Jrv[i]+1));
      if(eigvec.isempty())                                        //    if isempty(eigvec)
         eigvec = V*Mt;                                           //       eigvec = V*ritzvec(:,Jrv);
      else                                                        //    else
         eigvec.pset(-1,-1,eigvec.nc()+1,-1, V*Mt);               //       eigvec = [eigvec,V*ritzvec(:,Jrv)];
                                                                  //    end
      singblk[0] = 0;                                             //    singblk(1) = 0;

      // Check for convergence of Ritz vectors that do not occur in the desired part of the spectrum. The end 
      //    points of the dampening interval(s) will be reset in the main iteraion loop.

      if(abs((int)Jrev.size()-(int)Jrv.size())!=0)                //    if abs(length(Jrev)-length(Jrv)) ~= 0, deflate = 1; end
         deflate = true;
   }                                                              // end
}

// --------------------------------------------------------------------------------------------------------------- //
// This function computes the fast Leja points on the interval [-2,2].
// References: "Fast Leja Points", J. Baglama, D. Calvetti, and L. Reichel, ETNA (1998) Volume 7.
// --------------------------------------------------------------------------------------------------------------- //
void fastleja(std::vector<double> &fLejapts, int numshfts, std::vector<double> &rcandpts, sMat<int> &rindex, 
              std::vector<double> &rprd)
{
   int i,k,maxk=0;
   double maxval,d;

   // Number of previously generated fast Leja points.
   int nfLejapts = (int)fLejapts.size();               // nfLejapts = length(fLejapts);

   // Checking for errors in input argument(s).
   if(numshfts<1) { std::cerr << "Number of shifts is 0 or (< 0) in fast Leja routine."; return; }

   // Initialize the array fLejapts.
   if(nfLejapts<2)
   {
      if(nfLejapts<1) fLejapts.push_back(2.); else
      fLejapts[0] = 2.;                                // fLejapts(1) = 2;
      if((nfLejapts==0) & (numshfts==1)) return;       // if ((nfLejapts == 0) & (numshfts == 1)), return; end
      fLejapts.push_back(-2.);                         // fLejapts(2) = -2;
      if(rcandpts.empty()) rcandpts.push_back(0); else //
      rcandpts[0] = 0;                                 // rcandpts(1) = (fLejapts(1)+fLejapts(2))/2;
      if(rprd.empty()) rprd.push_back(-4.); else       //
      rprd[0] = -4.; rindex.zero(1,2);                 // rprd(1) = prod(rcandpts(1)-fLejapts);
      rindex(0,0) = 1; rindex(0,1) = 0;                // rindex(1,1) = 2; rindex(1,2) = 1;
      if((nfLejapts==0) & (numshfts==2)) return;       // if ((nfLejapts == 0) & (numshfts == 2)), return; end
      if((nfLejapts==1) & (numshfts==1)) return;       // if ((nfLejapts == 1) & (numshfts == 1)), return; end
      nfLejapts = 2; numshfts -= 2;                    // nfLejapts = 2; numshfts = numshfts-2;
   }

   if((int)fLejapts.size()<(numshfts+nfLejapts))
      fLejapts.insert(fLejapts.end(),(numshfts+nfLejapts)-fLejapts.size(),0.);
   if((int)rcandpts.size()<(numshfts+nfLejapts-1))
      rcandpts.insert(rcandpts.end(),(numshfts+nfLejapts-1)-rcandpts.size(),0.);
   if((int)rprd.size()<(numshfts+nfLejapts-1))
      rprd.insert(rprd.end(),(numshfts+nfLejapts-1)-rprd.size(),0.);
   // Main loop to calculate the fast Leja points.
   for(k=nfLejapts; k<(numshfts+nfLejapts); k++)       // for k = 1+nfLejapts:numshfts+nfLejapts
   {  maxval = vmax(vabs(rprd),maxk);                  //   [maxval,maxk]  = max(abs(rprd));
      fLejapts[k] = rcandpts[maxk];                    //   fLejapts(k)    = rcandpts(maxk);
      rindex(k-1,0) = k; rindex(k-1,1)=rindex(maxk,1); //   rindex(k-1,1)  = k; rindex(k-1,2) = rindex(maxk,2); 
      rindex(maxk,1) = k; rindex.reshape(k,2);         //   rindex(maxk,2) = k;
      rcandpts[maxk] = (fLejapts[rindex(maxk,0)]+      //   rcandpts(maxk) = (fLejapts(rindex(maxk,1))+
                        fLejapts[rindex(maxk,1)])/2.;  //                     fLejapts(rindex(maxk,2)))/2;
      rcandpts[k-1]  = (fLejapts[rindex(k-1,0)]+       //   rcandpts(k-1)  = (fLejapts(rindex(k-1,1))+
                        fLejapts[rindex(k-1,1)])/2.;   //                     fLejapts(rindex(k-1,2)))/2;
      d=1.; for(i=0;i<=(k-1);i++) d*=(rcandpts[maxk]-  //   rprd(maxk)     = prod(rcandpts(maxk)-fLejapts(1:k-1));
                        fLejapts[i]); rprd[maxk] = d;  //
      d=1.; for(i=0;i<=(k-1);i++) d*=(rcandpts[k-1]-   //   rprd(k-1)      = prod(rcandpts(k-1)-fLejapts(1:k-1));
                        fLejapts[i]); rprd[k-1] = d;   //
      for(i=0; i<(int)rprd.size(); i++)                //
         rprd[i] *= rcandpts[i]-fLejapts[k];           //   rprd           = rprd.*(rcandpts-fLejapts(k));
   }                                                   // end

 //sMat<int> rindtmp = rindex+1; 
 //std::cerr << "fLejapts = " << dispvect(fLejapts) << "\nrcandpts = " << dispvect(rcandpts);
 //std::cerr << "\nrprd = " << dispvect(rprd) << "\nrindex = " << rindtmp.display_full();
}

// --------------------------------------------------------------------------------------------------------------- //
// This function determines the endpoints of the dampening intervals, zeros (weighted Leja points or mapped Leja 
// points or mapped Leja points) and applies the zeros as shifts.
// --------------------------------------------------------------------------------------------------------------- //
void applyshifts(int blsz, sMat<double> &F, std::vector<double> &fLejapts, int iter, int K, std::vector<double> lcandpts,
                 double &leftendpt, double &leftintendpt, std::vector<double> &lprd, double &leftLejapt, double &lindex,
                 int &maxdpol, int &nbls, double &norlpol, std::vector<double> &rcandpts, double &rightendpt, double
                 &rightintendpt, eigVE<double> &ritz, std::vector<double> &rprd, double &rightLejapt, sMat<int> &rindex,
                 std::vector<int> &singblk, int sizint, int Tsz, sMat<double> &T, sMat<double> &V, int &flcount)
{
   int i,J;
   std::vector<double> pshifts; pshifts.reserve(nbls);
   sMat<double> Q,R,Mt;
   std::vector<int> Tsize,k;

   // Determine intervals for dampening.
   if(rightendpt==-1.)
   {
      rightintendpt = ritz.E[Tsz-sizint-1]; rightendpt = ritz.E[Tsz-1];
   }
   else
   {
      if(singblk[0]==0)
      {
         if(rightintendpt>ritz.E[Tsz-sizint-1]) 
            rightintendpt = ritz.E[Tsz-sizint-1];      // rightintendpt = min(ritz(Tsz-sizint),rightintendpt);
         if(rightendpt<ritz.E[Tsz-1]) 
            rightendpt = ritz.E[Tsz-1];                // rightendpt = max(ritz(Tsz),rightendpt);
      }
   }

   // Choose a value to scale the Leja polynomial to avoid underflow/overflow.
   if(fLejapts.empty()) norlpol = ritz.E[0];

   // Determine shifts. Compute fast Leja points on [-2,2] and map to the dampening interval.
   if( (flcount+nbls) > maxdpol ) 
   {
      flcount = 0; maxdpol = (int)fLejapts.size();     // flcount = 0; maxdpol = length(fLejapts);
   }
   if(flcount==(int)fLejapts.size())
      fastleja(fLejapts,nbls,rcandpts,rindex,rprd);    // [fLejapts,rcandpts,rindex,rprd] 
                                                       //     = fastleja(fLejapts,nbls,rcandpts,rindex,rprd);
   for(i=flcount; i<(nbls+flcount); i++)
      pshifts.push_back(fLejapts[i]);                  // pshifts = fLejapts(flcount+1:nbls+flcount);
   flcount += (int)pshifts.size();                     // flcount = flcount + length(pshifts);

   for(i=0; i<(int)pshifts.size(); i++)                // pshifts = (rightendpt-rightintendpt)/4*(pshifts-2) + rightendpt;
      pshifts[i] = (rightendpt-rightintendpt)/4.*(pshifts[i]-2.) + rightendpt;

   // Set-up Q to accumulate the Q product.
   Tsize = T.size(); sMat<double> Q_p(Tsize[0],Tsize[1]);
   for(i=0; i<vmax(Tsize); i++) Q_p(i,i) = 1.;         // Q_p = eye(size(T));

   // Apply nbls-1 shifts.
   for(J=nbls; J>=2; J--)
   {
      Mt = T; for(i=0; i<Tsz; i++) Mt(i,i) = Mt(i,i) - pshifts[nbls-J];
      R = qr(Mt,Q);                                    // [Q,R] = qr(T-pshifts(nbls-J+1)*eye(Tsz));
      Mt = Q; Mt.transpose(); T = Mt*(T*Q);            // T = Q'*(T*Q);
      T = tril(triu(T,-blsz),blsz);                    // T=tril(triu(T,-blsz),blsz);
      Q_p *= Q;                                        // Q_p = Q_p*Q;
   }
   V *= Q_p;                                           // V = V*Q_p;
   F = V.setp(-1,-1,blsz+1,2*blsz)*T.setp(blsz+1,2*blsz,1,blsz)  // F = V(:,blsz+1:2*blsz)*T(blsz+1:2*blsz,1:blsz)...
       + F*Q_p.setp(Tsz-blsz+1,Tsz,1,blsz);                      //       + F*Q_p(Tsz-blsz+1:Tsz,1:blsz);

   // Apply the last shift.  
   Mt = T.setp(1,blsz,1,blsz); for(i=0; i<blsz; i++); Mt(i,i) = Mt(i,i) - pshifts[nbls];
   V = F + V.setp(-1,-1,1,blsz)*Mt;                    // V = F + V(:,1:blsz)*(T(1:blsz,1:blsz)-pshifts(nbls)*eye(blsz));

   // Compute a random vector if a singular block occurs and no eigenvectors converged and modify starting vector(s). 
   // Details of singular blocks in the IRBL method are given in the paper "Dealing With Linear Dependence during the 
   // Iterations of the Restarted Block Lanczos Methods", J. Baglama, Num. Algs., 25, (2000) pp. 23-36. 
   // Modify starting vector (case when sigma='SE' and sigma='LE') if a singular block occurs and no eigenvectors converged.
   if(singblk[0]==1 && blsz>1)
   {
      // Compute which vector(s) of V need to be modified.
      for(i=1; i<(int)singblk.size(); i++) 
         k.push_back((singblk[i]-1)%blsz);             // k = rem(singblk(2:length(singblk))-1,blsz)+1; k = union(k,k);
      // Sorts and removes duplicates 
      k = setunion(k,k);
      Mt = randn(V.nr(),(int)k.size());
      for(i=0; i<(int)k.size(); i++)                   // V(:,k) = randn(size(V,1),length(k));
         V.pset(-1,-1,k[i],k[i], Mt.setp(-1,-1,i+1,i+1));
   }
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates a few of the lowest eigen[values,vectors] of a matrix using the block Lanczos method
// --------------------------------------------------------------------------------------------------------------- //
eigVE<double> irbleigs(sMat<double> & M, int K)
{
   // This function is a translation from MATLAB to C++ of the irbleigs() function by James Baglama, Daniela Calvetti
   //    and Lothar Reichel [ IRBL: An Implicitly Restarted Block Lanczos Method for large-scale Hermitian eigenproblems,
   //    J. Baglama, D. Calvetti, and L. Reichel, SIAM J. Sci. Comput. v24, No. 5, (2003), pp. 1650-1677 ]

   eigVE<double> retval;
   if(!M.issymm()) return retval;                               // Checks that matrix is symmetric, otherwise returns empty
   double eps = DBL_EPSILON;

   int blsz = 3;                                                // Block size of the Lanczos tridiagonal matrix
   // bool cholM = false;                                       // Cholesky factorisation not available
   int maxit = 2500;                                            // Number of iterations
   // sigma = 'SE'                                              // Finds smallest eigenvalues
   // endpt = 'MON'                                             // Interior end pt. location. Set so dampening interval inc.
   int maxdpol = 200;                                           // Maximum degree of dampening
   int sizint = 1;                                              // Size of dampening interval (1=consecutive)
   int nbls = (int)ceil((K+sizint+blsz)/(double)blsz+0.1);      // Number of blocks in Lanczos matrix
   double tol = 1e-6;                                           // Tolerance for convergence: ||Ax-Ex||^2 <= tol*||A||^2
   // zertyp = 'ML'                                             // Use mapped Leja points on [-2,2]

   int n = M.nc();
   sMat<double> V = randn(n,blsz);                              // Starting vectors
 //retval.V.zero(n,K);                                          // Vector of converged eigenvectors (eigvec)
   sMat<double> F(n,blsz);                                      // Block Lanczos decomposition: M*V = V*T + F*E.transpose()
   sMat<double> T(n,n);                                         //   T is tridiagonal, V is eigenvalues, E eigenvectors

   // Initialization and description of local variables.
   int conv = 0;                        // Number of desired eigenvalues that have converged.
   bool deflate = false;                // Value used to determine if an undesired eigenvector has converged.
   bool computvec = true;               // Determines whether to compute vectors or not
   std::vector<double> eigresdisp;      // Holds the residual values of converged eigenvalues.
 //std::vector<double> eigval;          // Array used to hold the converged desired eigenvalues. (retval.E)
   int iter = 1;                        // Main loop iteration count.
   std::vector<double> fLejapts;        // Stores fast Leja points.
   std::vector<double> lcandpts;        // Stores "left" candidate points.
   double lindex;                       // Index for lcandpts.
   std::vector<double> lprd;            // Stores the Leja polynomial product for "left" candidate points.
   double leftendpt = -1.;              // Stores the left most endpoint of the dampening interval.
   double leftintendpt = -1.;           // Stores the interior left endpoint of the dampening interval.
   double leftLejapt;                   // Place holder in fast leja routines.
   int mprod = 0;                       // The number of matrix vector products.
   double norlpol;                      // Normalizes the Leja polynomial to avoid underflow/overflow.
   int numbls = nbls;                   // Initial size of the Lanczos blocks.
   std::vector<double> pritz;           // Holds previous iteration of Ritz values. Used to determine stagnation.
   bool ritzconv = false;               // Boolean to determine if all desired eigenvalues have converged.
   std::vector<double> rcandpts;        // Stores "right" candidate points.
   sMat<int> rindex;                    // Index for rcandpts.
   std::vector<double> rprd;            // Stores the Leja polynomial product for "right" candidate points.
   double rightendpt = -1.;             // Stores the right most endpoint of the dampening interval.
   double rightintendpt = -1.;          // Stores the interior right endpoint of the dampening interval.
   double rightLejapt;                  // Place holder in fast Leja routines.
   std::vector<int> singblk;            // Integer values used to indicate singular block(s) in blanz.
   double sqrteps = sqrt(eps);          // Square root of machine tolerance used in convergence testing.
   int flcount = 0;                     // Count used to determine when the maximum number of shifts is reached.
   int Tsz;                             // Size of tridiagonal matrix
   eigVE<double> ritz;                  // Ritz values
   std::vector<int> J;                  // Index for sorted ritz values
   sMat<double> Mt;                     // Temporary matrix
   std::vector<double> residuals;       // Residuals of Ritz values
   double Rmax = -1.;                   // Holds the absolute maximum Ritz value

   int i,tmp;

   while(iter <= maxit)
   {
    //std::cerr << "iter = " << iter << "\tconv = " << conv << "\n";
      // Check and deflate the number of Lanczos blocks if possible. Do not deflate if searching for interior 
      // eigenvalues or if the users has inputed eigenvectors. This causes nbls to increase.
      if(nbls>2)
      {
         nbls = numbls - (int)floor(retval.V.nc()/blsz);
         if(sizint >= nbls*blsz-(fabs(K-retval.V.nc())))
         {
            tmp = (int)floor((sizint+fabs(K-retval.V.nc()))/blsz);
            nbls = (tmp>numbls) ? tmp : numbls;
         }
      }

      // Compute the block Lanczos decomposition.
      blanz(M,K,V,blsz,retval.V,mprod,n,nbls,singblk,sqrteps,tol,F,T);

      // Determine number of blocks and size of the block tridiagonal matrix T.
      Tsz = T.nr(); nbls = Tsz/blsz;
      if(floor(nbls) != nbls)
      {
         // Reset the starting matrix V and re-compute the block Lanczos decomposition if needed.
         Mt = randn(n,blsz); blanz(M,K,Mt,blsz,retval.V,mprod,n,nbls,singblk,sqrteps,tol,F,T);
         Tsz = T.nr(); nbls = Tsz/blsz;
         if(floor(nbls) != nbls)
         {
            std::cerr << "blanz in irbleigs.m returns an incorrect matrix T.";
            return retval;
         }
      }

      // Compute eigenvalues and eigenvectors of the block tridiagonal matrix T, with sorted eigenvalues
      ritz = eigT(T); 

      /* Sort the eigenvalues from smallest to largest.
      ritz.E = vsort(ritz.E,J);                        //   [ritz,J] = sort(real(ritz)); ritzvec = ritzvec(:,J);
      Mt.zero(ritz.V.nr(),ritz.V.nc());
      for(i=0; i<(int)J.size(); i++) Mt.pset(-1,-1,i+1,i+1,ritz.V.setp(-1,-1,J[i],J[i])); */

      // Reached maximum number of iterations, exit main loop.
      if(iter>=maxit) break;

      // Compute the residuals for all ritz values. 
      Mt = ( F*ritz.V.setp(Tsz-(blsz-1),Tsz,-1,-1) ); Mt ^= Mt; // ^ is element-wise multiplication
      residuals = msum(Mt); for(i=0; i<(int)residuals.size(); i++) residuals[i] = sqrt(residuals[i]);

      // To test convtests()
/*    std::cout << "eigval = " << dispvect(retval.E) << ";\neigresdisp = " << dispvect(eigresdisp) << ";\n";
      if(retval.V.isempty()) std::cout << "eigvec = [];\n"; else std::cout << "eigvec = " << retval.V.display_full();
      std::cout << "pritz = " << dispvect(pritz) << ";\nresiduals = " << dispvect(residuals) << ";\n";
      std::cout << "ritz = " << dispvect(ritz.E) << ";\nsingblk = " << dispvect(singblk) << ";\n";
      if(ritz.V.isempty()) std::cout << "ritzvec = [];\n"; else std::cout << "ritzvec = " << ritz.V.display_full();
      if(V.isempty()) std::cout << "V = [];\n"; else std::cout << "V = " << V.display_full();
      std::cout << "[conv,deflate,eigval,eigvec,eigresdisp,pritz,ritzconv,singblk] = convtests('T'"; // computvec
      std::cout << ","<<conv<<","<<deflate<<",0,eigval(:),eigvec,eigresdisp,"<<iter<<","<<K<<",[],pritz(:),residuals(:)";
      std::cout << ",ritz(:),"<<ritzconv<<",ritzvec,'SE',singblk(:),"<<sqrteps<<","<<tol<<","<<Tsz<<",V);\n";

*/    // Convergence tests.
      convtests(computvec,conv,deflate,retval.E,retval.V,eigresdisp,iter,K,pritz,residuals,ritz,ritzconv,
                singblk,sqrteps,tol,Tsz,V,Rmax);

/*    if(deflate) std::cout << "clear global;\n";
      std::cout << "myeigval = " << dispvect(retval.E) << ";\nmyeigvec = " << retval.V.display_full();
      std::cout << "myeigresdisp = " << dispvect(eigresdisp) << ";\nmypritz = " << dispvect(pritz) << ";\n";
      std::cout << "mysingblk = " << dispvect(singblk) << ";\n";
      std::cout << "disp('deflate\teigval\teigvec\teigresdisp\tpritz\tsingblk')\n[deflate-("<<deflate<<") ";
      std::cout << "sum(abs(eigval(:)-myeigval(:))) sum(sum(abs(eigvec-myeigvec))) sum(abs(eigresdisp(:)-myeigresdisp(:))) ";
      std::cout << "sum(abs(pritz(:)-mypritz(:))) sum(abs(singblk(:)-mysingblk(:)))]\nclear;\n";

*/    // If all desired Ritz values converged then exit main loop.
      if(ritzconv) break;

/*    // To test applyshifts() and fastleja() 
      sMat<int> rindtm = rindex + 1;
      std::cout << "F = " << F.display_full() << "fLejapts = " << dispvect(fLejapts) << ";\nlcandpts = " << dispvect(lcandpts) << ";\n";
      std::cout << "lprd = " << dispvect(lprd) << ";\nrcandpts = " << dispvect(rcandpts) << ";\nritz = " << dispvect(ritz.E) << ";\n";
      std::cout << "ritzvec = " << ritz.V.display_full() << "rprd = " << dispvect(rprd) << ";\nrindex = " << rindtm.display_full(); 
      std::cout << "singblk = " << dispvect(singblk) << ";\nT = " << T.display_full() << "V = " << V.display_full(); 
      std::cout << "[fLejapts,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,nval,rcandpts,rightendpt,"
      std::cout << rightintendpt,rprd,rightLejapt,rindex,V,flcount] = ...\n";
      std::cout << "applyshifts(3,'MON',F,fLejapts,"<<iter<<","<<K<<",lcandpts,"; 
      if(leftendpt==-1) std::cout << "[],"; else std::cout << leftendpt << ",";
      if(leftintendpt==-1) std::cout << "[]"; else std::cout << leftintendpt; std::cout << ",lprd,"<<leftLejapt<<",[],"<<maxdpol;
      std::cout << ","<<nbls<<","<<norlpol<<",[],rcandpts,"; if(rightendpt==-1) std::cout << "[],"; else std::cout <<rightendpt << ",";
      if(rightintendpt==-1) std::cout << "[]"; else std::cout <<rightintendpt; 
      std::cout << ",ritz,ritzvec,rprd,"<<rightLejapt<<",rindex,'ML','SE',";
      std::cout << "singblk,"<<sizint<<","<<sqrteps<<","<<Tsz<<",T,V,"<<flcount<<");\n";
*/    // Determine dampening intervals, Leja points, and apply Leja zeros as shifts.
      applyshifts(blsz,F,fLejapts,iter,K,lcandpts,leftendpt,leftintendpt,lprd,leftLejapt,lindex,maxdpol,nbls,norlpol,rcandpts,
                  rightendpt,rightintendpt,ritz,rprd,rightLejapt,rindex,singblk,sizint,Tsz,T,V,flcount);
/*    std::cout << "myfLejapts = " << dispvect(fLejapts) << ";\nmyrcandpts = " << dispvect(rcandpts) << ";\n";
      rindtm = rindex + 1; std::cout << "myrprd = " << dispvect(rprd) << ";\nmyrindex = " << rindtm.display_full(); 
      std::cout << "disp('fLejapts\tnorlpol\trcandpts\trightenpt\trightintendpt\trprd\trightLejapt\trindex\tV\tflcount')\n";
      std::cout << "[sum(sum(abs(fLejapts-myfLejapts))) norlpol-("<<norlpol<<") sum(sum(abs(rcandpts-myrcandpts))) rightendpt-("<<rightendpt;
      std::cout << ") rightintendpt-("<<rightintendpt<<") sum(sum(abs(rprd-myrprd))) rightLejapt-("<<rightLejapt;
      std::cout << ") sum(sum(abs(rindex-myrindex))) flcount-("<<flcount<<")]\n";
*/  //std::cout << "T = " << T.display_full() << "V = " << V.display_full(); 

      // Check to see if an undesired Ritz vector has been chosen to be deflated. If so,
      // reset the endpoints of the dampening interval(s) on the next interation.
      if(deflate)
      {
         deflate = false; Rmax = -1.; leftintendpt = -1; rightintendpt = -1; leftendpt = -1; rightendpt = -1;
      }

      // Update the main iteration loop count.
      iter++;
   }

   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: 
//    g++ states.cpp cfp.cpp njsyms.cpp so_cf.cpp coulomb.cpp ic_hmltn.cpp maths_irbl.cpp && ./a.out
// --------------------------------------------------------------------------------------------------------------- //
/*sMat<double> ic_hmltn(int n, std::vector<double> F, double xi, std::vector< std::vector<double> > B);
int main(int argc, char *argv[])
{
   int n;
   if(argc>1) { n = atoi(argv[1]); } else { n = 2; }

   std::vector<double> Fs(4,1.);  // F = [1 1 1 1];
   std::vector<double> B2(5,0.);  
   std::vector<double> B4(9,0.);
   std::vector<double> B6(13,0.);
   std::vector< std::vector<double> > B;

   B4[4] = 1.; B4[8] = sqrt(5./14.); B6[6] = 1.; B6[10] = -sqrt(7./2.);
   B.push_back(B2); B.push_back(B4); B.push_back(B6);
   sMat<double> Z;
   sMat<double> A = ic_hmltn(Z,n,Fs,1.,B);

*/ // K=6; opt.sigma='SE'; opt.k=(Nlev); opt.blsz=6; opt.nbls=ceil((Nlev+7)/6+0.1); opt.maxit=2500;
 /*int K = 6; int blsz = 6; int nbls = (int)ceil((K+7.)/6.+0.1); double tol = 1e-6; //int sizint = 1;

   // Testing blanz() - run ./ic_hmltn > out.m; and then run >> out in matlab
   int nc = A.nc();
   sMat<double> V = randn(nc,blsz);
   sMat<double> eigvec,F,T;
   int mprod = 0;
   std::vector<int> singblk;
   std::cout << "A = " << A.display_full(); std::cout << "V = " << V.display_full();
   blanz(A,K,V,blsz,eigvec,mprod,nc,nbls,singblk,sqrt(DBL_EPSILON),tol,F,T);
   std::cout << "myF = " << F.display_full(); std::cout << "myT = " << T.display_full();
   std::cout << "[F,T,V,blsz,mprod,singblk] = blanz(A," << K << ",[],[],V," << blsz << ",[],[]," << mprod << "," << nc << "," << nbls << ",[],"<< sqrt(DBL_EPSILON) << "," << tol << ",V_orth)\n";
   std::cout << "[sum(sum(abs(F-myF)))/sum(sum(abs(F))) sum(sum(abs(T-myT)))/sum(sum(abs(T)))]";
 
   // Testing irbleigs() - run ./ic_hmltn > out.m; and then run >> out in matlab
*//*int K = 6;
   std::cout << "A=ic_hmltn('f" << n << "',[1 1 1 1],1,{zeros(1,5) [0 0 0 0 1 0 0 0 sqrt(5/14)] [0 0 0 0 0 0 1 0 0 0 -sqrt(7/2) 0 0]});\n";
   std::cout << "opt.sigma='SE'; opt.k=" << K << "; opt.blsz=3; opt.nbls=ceil(("<<K<<"+1+3)/3+0.1); opt.maxit=2500;\n";
   eigVE<double> dg = irbleigs(A,K);
   std::cout << "myE = " << dispvect(dg.E) << ";\nmyV = " << dg.V.display_full();
   std::cout << "[V,E,prginf] = irbleigs(A,opt); [sum(sum(abs(myV-V))) sum(sum(abs(myE'-diag(E))))]";

   // Test times!
 //eigVE<double> dg = eig(A);

   return 0;
}*/
