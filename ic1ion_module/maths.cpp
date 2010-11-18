/* maths.cpp
 *
 * This file implements mathematical and matrix classes with similar syntax to Matlab, and may be used to
 * translate programs from Matlab to C++
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
#include<ctime>
#include<cstdlib>                      // For rand()

// --------------------------------------------------------------------------------------------------------------- //
// Uses bit manipulation to determine sign of a floating point number
// --------------------------------------------------------------------------------------------------------------- //
float sign(float val)
{
/* // Clears all bits except sign bit
   *((int*)&val) &= 0x80000000;
   // Sets mantissa to 1. and exponent to 0;
   *((int*)&val) |= 0x3f800000;
   return val; */
   return val<0?-1.:1.;
}
double sign(double val)
{
/* // Clears all bits except sign bit
   *((long int*)&val) &= 0x8000000000000000;
   // Sets mantissa to 1. and exponent to 0;
   *((long int*)&val) |= 0x3ff0000000000000;
   return val; */        // Without optimisation, takes about 11.6s for 1e9 iterations. Fails with optimisation!
   return val<0?-1.:1.;  // Without optimisation, takes about 9.2s for 1e9 iterations! 5.5s with optimisation...
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the eigenvalues and eigenvectors of a symmetric matrix using the LAPACK routine dsyev
// --------------------------------------------------------------------------------------------------------------- //
eigVE<double> eig(const sMat<double> & M)
{
   eigVE<double> retval;
   if(!M.issymm()) return retval;     // Checks that matrix is symmetric first otherwise returns empty

   // We must pass all parameters to Fortran routines by reference. Hence we must first set up variables for each
   /* parameters. The following are parameters to work out a good block size for the algorithm using ILAENV
   int ispec = 1;          // Calculates optimum block size
   char *name = "SSYTRD";  // Calculates using the tri-diagonal reduction of a real symmetric matrix
   int unused = -1;        // For unused parameters N3, N4. We use opts==jobz, N1==lda, N2==n (see below) */
   
   // Paramaters for DSYEV itself
   int lda = M.nc();
   int n = M.nr();
   char jobz = 'V';        // Job type - 'N' is eigenvalues only; 'V' calculates eigenvectors as well
   char uplo = 'U';        // States whether to use the upper ('U') or lower triangle ('L')
   int info = 0;           // Success flag on return. 0 is success. <0 is wrong input. >0 means not converged
   double *work;           // A C-array for the workspace
   int lwork;              // The length of the workspace

   // Declaration for the C-array holding the matrix on input and the eigenvectors on output
   double *m = M.f_array();
   double *eigval = (double*)malloc(n*sizeof(double));

   /* Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   lwork = ( F77NAME(ilaenv)(&ispec, name, &jobz, &lda, &n, &unused, &unused) +2 ) * n; */
   lwork = 4*n;
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dsyev)(&jobz, &uplo, &n, m, &lda, eigval, work, &lwork, &info);

   if(info==0)
   {
      retval.E = f2vec(eigval,n);
      retval.V = f2mat(m,lda,n);
   }
   else
   {
      std::cerr << "eig(): DSYEV return error code " << info << "\n";
   }

   free(eigval); free(m); free(work);
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the single value decomposition of a real matrix
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> svd(const sMat<double> & M)
{
   std::vector<double> retval;

   // The singular values q[] of the matrix A is given by:		       T
   // 			 T     T				A = U diag(q) V
   // 	 with:		U U = V V = I

   // Paramaters for DGESVD
   char jobu = 'S';        // Job type - 'A' returns all columns of U or V^T [jobu/jobvt]
   char jobvt = 'N';       //          - 'S' first min(m,n) columns of U or V^T ("economy" mode)
                           //          - 'O' fisrt min(m,n) columns of U or V^T overwritten on A
                           //          - 'N' don't compute U or V^T
   int m = M.nr();         // Number of rows of input matrix A
   int n = M.nc();         // Number of columns of input matrix A
   double *a = M.f_array();// Input matrix A (C-array version of M)
   int lda = m;            // Leading dimension of A >=max(1,m)
   double *s;              // On output the singular values of A, dimension min(m,n)
   double *u;              // On output the left singular vectors, U
   int ldu = m;            // Leading dimension of U
 //double *vt;             // On output the right singular vectors V^T
   int ldvt = n;           // Leading dimension of V^T
   double *work;           // A C-array for the workspace
   int lwork;              // The length of the workspace, LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
   int info = 0;           // Success flag on return. 0 is success. <0 is wrong input. >0 means not converged

   // Allocates memory to output arrays
   s  = (double*)malloc( (n<m?n:m) * sizeof(double));    // dimension min(m,n)
   u  = (double*)malloc( m*m * sizeof(double));          // dimension m*m

   // Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   lwork = 3*(m<n?m:n)+(m>n?m:n); lwork = (lwork>(5*(m>n?m:n))) ? lwork : (5*(m>n?m:n));
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, u, &ldvt, work, &lwork, &info);

   if(info==0)
      retval = f2vec(s,n<m?n:m);
   else
      std::cerr << "svd(): DGESVD return error code " << info << "\n";

   // Frees memory of output arrays
   free(s); free(u); free(work);
   return retval;
}
std::vector<double> svd(const sMat<double> & M, sMat<double> & V)
{
   std::vector<double> retval;

   // Calculates the right singular vectors as well this time

   // Paramaters for DGESVD
   char jobu = 'N',jobvt = 'A';
   int m = M.nr(), n = M.nc(), lda = m, ldu = m, ldvt = n, info = 0;
   double *a = M.f_array();
   double *s,*vt,*work;
   // Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   int lwork = 3*(m<n?m:n)+(m>n?m:n); lwork = (lwork>(5*(m>n?m:n))) ? lwork : (5*(m>n?m:n));

   // Allocates memory to arrays
   s  = (double*)malloc( (n<m?n:m) * sizeof(double));    // dimension min(m,n)
   vt = (double*)malloc( n*n * sizeof(double));          // dimension n*n
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, vt, &ldu, vt, &ldvt, work, &lwork, &info);

   if(info==0)
   {  retval = f2vec(s,n<m?n:m); V = f2mat(vt,n,n); }
   else
      std::cerr << "svd(): DGESVD return error code " << info << "\n";

   // Frees memory of output arrays
   free(s); free(vt); free(work);
   return retval;
}
std::vector<double> svd(const sMat<double> & M, sMat<double> & U, sMat<double> & V)
{
   std::vector<double> retval;

   // Calculates both the left and right singular vectors this time

   // Paramaters for DGESVD
   char jobu = 'A',jobvt = 'A';
   int m = M.nr(), n = M.nc(), lda = m, ldu = m, ldvt = n, info = 0;
   double *a = M.f_array();
   double *s,*u,*vt,*work;
   // Works out the length of the workspace using the block size calculated by ILAENV: (2*blz)*n
   int lwork = 3*(m<n?m:n)+(m>n?m:n); lwork = (lwork>(5*(m>n?m:n))) ? lwork : (5*(m>n?m:n));

   // Allocates memory to arrays
   s  = (double*)malloc( (n<m?n:m) * sizeof(double));    // dimension min(m,n)
   u  = (double*)malloc( m*m * sizeof(double));          // dimension m*m
   vt = (double*)malloc( n*n * sizeof(double));          // dimension n*n
   work = (double*) malloc(lwork*sizeof(double));

   F77NAME(dgesvd)(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info);

   if(info==0)
   {  retval = f2vec(s,n<m?n:m); U = f2mat(u,m,m); V = f2mat(vt,n,n); }
   else
      std::cerr << "svd(): DGESVD return error code " << info << "\n";

   // Frees memory of output arrays
   free(s); free(u); free(vt); free(work);
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the QR decomposition of a rectangular matrix
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> qr(const sMat<double> & M, sMat<double> & Q)
{
   sMat<double> retval;    // The matrix R to be returned.

   // Parameters for DGEQRF
   int m = M.nr();         // Number of rows of input matrix A
   int n = M.nc();         // Number of columns of A
   double *a = M.f_array();// Input matrix A as an array
   int lda = m;            // Leading dimension of A
   double *tau;            // On output, the scalar factors of elementary reflectors
   double *work;           // Workspace
   int lwork = -1;         // Length of workspace. If lwork=-1, calculates optimal workspace
   int info = 0;           // Success flag on return. 0 is success. <0 is wrong input. >0 means not converged
   // Extra parameters for DORGQR
   int k = m<n?m:n;        // The number of elementary reflectors

   // Allocates memory to array initially
   tau = (double*)malloc( (m<n?m:n) * sizeof(double) );
   work = (double*)malloc( 4*n * sizeof(double) );      // Initial value only. We use the workspace query to find optimum
   // Calculates optimum workspace size
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
   if(info==0) lwork = (int)work[0]; else { std::cerr << "qr(): DGEQRF could not determine optimum workspace size\n"; return retval; }
   work = (double*)realloc( work, lwork*sizeof(double) );
   if(work==NULL) { std::cerr << "qr(): error reallocating workspace memory\n"; return retval; }

   // Calculates the QR factorisation
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

   if(info==0)
   {
      retval = triu(f2mat(a,m,n));
      // Calculates the orthogonal matrix Q from the output of DGEQRF
      F77NAME(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
      if(info==0) Q = f2mat(a,m,n); else std::cerr << "qr(): DGEQRF return error code " << info << "\n";
   }
   else
      std::cerr << "qr(): DGEQRF return error code " << info << "\n";

   free(tau); free(work);
   return retval;
}
sMat<double> qr(const sMat<double> & M, sMat<double> & Q, int zero)
{
   sMat<double> retval;    // The matrix R to be returned.

   // Parameters for DGEQRF and DORGQR
   int m = M.nr(), n = M.nc(), lda = m, lwork = -1, info = 0, k = m<n?m:n;
   double *a, *tau, *work;
   if(m>n) { m = n; lda = n; retval = M; retval.resize(1,n,1,n); a = retval.f_array(); } else a = M.f_array();
   // Allocates memory to array initially
   tau = (double*)malloc( (m<n?m:n) * sizeof(double) );
   work = (double*)malloc( 4*n * sizeof(double) );      // Initial value only. We use the workspace query to find optimum
   // Calculates optimum workspace size
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
   if(info==0) lwork = (int)work[0]; else { std::cerr << "qr(): DGEQRF could not determine optimum workspace size\n"; return retval; }
   work = (double*)realloc( work, lwork*sizeof(double) );
   if(work==NULL) { std::cerr << "qr(): error reallocating workspace memory\n"; return retval; }

   // Calculates the QR factorisation
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

   if(info==0)
   {
      retval = triu(f2mat(a,m,n));
      // Calculates the orthogonal matrix Q from the output of DGEQRF
      F77NAME(dorgqr)(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
      if(info==0) Q = f2mat(a,m,n); else std::cerr << "qr(): DGEQRF return error code " << info << "\n";
   }
   else
      std::cerr << "qr(): DGEQRF return error code " << info << "\n";

   return retval;
}
sMat<double> qr(sMat<double> & M, int zero)
{
   sMat<double> retval;    // The matrix R to be returned.

   // Parameters for DGEQRF
   int m = M.nr(), n = M.nc(), lda = m, lwork = -1, info = 0;
   double *a, *tau, *work;
   if(m>n) { m = n; lda = n; retval = M; retval.resize(1,n,1,n); a = retval.f_array(); } else a = M.f_array();
   // Allocates memory to array initially
   tau = (double*)malloc( (m<n?m:n) * sizeof(double) );
   work = (double*)malloc( 4*n * sizeof(double) );      // Initial value only. We use the workspace query to find optimum
   // Calculates optimum workspace size
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);
   if(info==0) lwork = (int)work[0]; else { std::cerr << "qr(): DGEQRF could not determine optimum workspace size\n"; return retval; }
   work = (double*)realloc( work, lwork*sizeof(double) );
   if(work==NULL) { std::cerr << "qr(): error reallocating workspace memory\n"; return retval; }

   // Calculates the QR factorisation
   F77NAME(dgeqrf)(&m, &n, a, &lda, tau, work, &lwork, &info);

   if(info==0)
      retval = triu(f2mat(a,m,n));
   else
      std::cerr << "qr(): DGEQRF return error code " << info << "\n";

   free(tau); free(work);
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the absolute values of each non zero element in a matrix
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> mabs(const sMat<double> & M)
{
   std::vector< std::vector<int> > nz;
   int i,sz;
   sMat<double> r = M;

   nz = M.find(); sz = (int)nz.size();
   for (i=0; i<sz; i++)
      r(nz[i][0],nz[i][1]) = fabs(M(nz[i][0],nz[i][1]));

   return r;
}
std::vector<double> vabs(const std::vector<double> & V)
{
   std::vector<double> r(V.size(),0.);
   int i;

   for (i=0; i<(int)V.size(); i++)
      r[i] = fabs(V[i]);

   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Calculates the norm of a matrix
// --------------------------------------------------------------------------------------------------------------- //
double norm(const sMat<double> & M)
{
   return vmax<double>(svd(M));
}
double norm(const sMat<double> & M, normp p)
{
   sMat<double> Mt = M; Mt.transpose();
   switch(p)                                                       // To be compatible with the matlab function
   {
      case 1: return vmax<double>(msum<double>(mabs(M))); break;   // The 1-norm
      case 2: return norm(M); break;                               // Largest singular value
      case 3: return vmax<double>(msum<double>(mabs(Mt)));         // The infinity norm (normp inf==3, see maths.hpp)
              break;
      case 4: return sqrt(vsum<double>(diag<double>(Mt*M)));       // The frobenius norm (normp fro==4, see maths.hpp)
              break;
      default: return 0.;
   }
   return 0.;
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns a (mxn) matrix of normally distributed random numbers
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> randn(int m, int n)
{
   sMat<double> M(m,n);
   int i,j,k=0;
   srand(time(NULL));
   double x1=0,x2=0,y=0;
   double maxrand = RAND_MAX/2.;

   // This function uses the C++ standard function rand() to generate uniformly distributed random numbers
   //    between [0,1], and then the polar Box-Muller transformation to generate pairs of normally distributed
   //    random numbers (corresponding to sin(angle) and cos(angle) respectively). Hence we need some trickery to
   //    use both the pair, and not waste half the computing time! [ mean=0, s.d.=1 ]

   for(i=0; i<m; i++)
      for(j=0; j<n; j++)
      {
         if(k%2 == 0)		// First of pair, generate new pair.
	 {
	    do
	    {
               x1 = rand()/maxrand - 1.;
               x2 = rand()/maxrand - 1.;
               y  = x1*x1 + x2*x2;
            } while (y>=1.);
	    y = sqrt( (-2. * log(y)) / y);
	    M(i,j) = x1*y;
	 }
	 else			// Second of pair.
	    M(i,j) = x2*y;
	 k++;
      }

   return M;
}

// --------------------------------------------------------------------------------------------------------------- //
// Returns an orthonormal basis for a matrix. Translation of the MATLAB routine using svd.
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> orth(const sMat<double> & M)
{
   sMat<double> U,Q,V;
   int i,r=0;
   sMat<double> Mt = M; Mt.transpose();
   int m = M.nr(); int n = M.nc();                       // [m,n] = size(A);
   std::vector<double> s = svd(Mt,U); U.resize(1,m,1,n); // [U,S,V] = svd(A,0);
   double tol;

   if(m>1) {}                                            // if m > 1, s = diag(S);
      else if(m==1) { s.erase(s.begin()+1,s.end()); }    //    elseif m == 1, s = S(1);
      else { s.clear(); s.push_back(0.); }               //    else s = 0;
                                                         // end
   tol = (double)vmax<int>(M.size()) * vmax<double>(s)   // tol = max(m,n) * max(s) * eps(class(A)); 
           * DBL_EPSILON;
   for(i=0; i<(int)s.size(); i++)
      if(s[i]>tol) r++;                                  // r = sum(s > tol);

   Q = U.setp(-1,-1,1,r);                                // Q = U(:,1:r);

   return Q;
}

// --------------------------------------------------------------------------------------------------------------- //
// Estimates the number of linearly independent rows or columns of a matrix. Translation of the MATLAB function.
// --------------------------------------------------------------------------------------------------------------- //
int rank(const sMat<double> & M, double tol)
{
   std::vector<double> s = svd(M);
   int i,r=0;
   for(i=0; i<(int)s.size(); i++)
      if(s[i]>tol) r++;

   return r;
}

// --------------------------------------------------------------------------------------------------------------- //
// Removes entries in matrix less than machine epsilon
// --------------------------------------------------------------------------------------------------------------- //
void rmzeros(sMat<double> & M)
{
   std::vector< std::vector<int> > nz;
   int i,sz;
   nz = M.find(); sz = (int)nz.size();
   for (i=0; i<sz; i++)
      if(fabs(M(nz[i][0],nz[i][1]))<DBL_EPSILON*1000) M.del(nz[i][0],nz[i][1]);
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts a 1D C-array into a std::vector
// --------------------------------------------------------------------------------------------------------------- //
std::vector<double> f2vec(double *v, int n)
{
   int i;
   std::vector<double> retval(n,0.);
   for (i=0; i<n; i++) retval[i] = v[i];
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts a 2D C-array into an sMat
// --------------------------------------------------------------------------------------------------------------- //
sMat<double> f2mat(double *M, int m, int n)
{
   int i,j;
   sMat<double> retval(m,n);
   for(j=0; j<n; j++)
      for(i=0; i<m; i++)
         if(fabs(M[m*j+i])>DBL_EPSILON*1000) retval(i,j) = M[m*j+i];
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// Converts two double precision sMat to a complex C-array
// --------------------------------------------------------------------------------------------------------------- //
complexdouble* zmat2f(sMat<double> &r, sMat<double> &i)
{
   sMat<double> tmp = r+i;
   std::vector< std::vector<int> > u = tmp.find/*upper*/();
   int j, rows=tmp.nr();
   // Allocates an _r*_c array and initiallises all elements to zero.
   complexdouble x, *retval = (complexdouble*)calloc(rows*tmp.nc(),sizeof(complexdouble));
   for (j=0; j<(int)u.size(); j++)
   {
      x.r = r(u[j][0],u[j][1]); x.i = i(u[j][0],u[j][1]);
      retval[rows*u[j][1]+u[j][0]] = x;
   }
   return retval;
}

// --------------------------------------------------------------------------------------------------------------- //
// For testing the rest of the code! - Uncomment and compile: g++ maths.cpp && ./a.out
// --------------------------------------------------------------------------------------------------------------- //
/*int main(int argc, char *argv[])
{
   sMat<double> l,m,n,a,b,c;//(3,3);
   std::vector< std::vector<int> > nz;
   std::vector<int> size;
   std::vector<double> v,v2;
   int i,sz;
   eigVE<double> dg;
  
   std::cout << "Empty matrix is " << m; // << "\n";

   // Tests mset("") member function
   m.mset("1 0 0 ; 0 1 0 ; 0 0 1");
   nz = m.find();
   std::cout << "m.mset(\"1 0 0 ; 0 1 0 ; 0 0 1\"); m = " << m;
   std::cout << "Setting (2,3) to 2:\n"; m(1,2) = 2; std::cout << m;

   // Test find() member function.
   std::cout << "Testing find(): nz = m.find(); gives:";
   sz = (int)nz.size(); //std::cout << sz << "\n";
   for (i=0; i<sz; i++)
   {
      std::cout << "(" << nz[i][0] << "," << nz[i][1] << ")\t" << "\n";
   }
   m.mset("1/3 sqrt(2) 2 ; 2 1/2 2 ; sqrt(2/4) 2 1");
   nz = m.find(); sz = (int)nz.size();
   for (i=0; i<sz; i++)
   {
      std::cout << "(" << nz[i][0] << "," << nz[i][1] << ")\t" << m(nz[i][0],nz[i][1]) << "\n";
   }

   // Tests the iostream output friend function.
   n.mset("1*3 sqrt(2) 2/6 ; 2 1/2 2 ; sqrt(2/4) 2+1 1-6");
   std::cout << "Testing iostream: cout << n gives:\n" << n;
   //n >> l;

   // Tests size() and issquare()
   size = m.size(); std::cout << "Matrix size is " << size[0] << "x" << size[1] << "\n";
   std::cout << "Matrix issquare is " << (m.issquare() ? "true" : "false" ) << "\n";
  
   // Tests overloaded operators
   m+=n;
   std::cout << m; 
   m-=m;
   std::cout << m; 
   m.mset("1 0 2 ; 0-1 3 1"); n.mset("3 1 ; 2 1 ; 1 0");
   std::cout << m; std::cout << " x \n"; std::cout << n;
   std::cout << "null value is " << m(0,1) << "\n";
   m*=n;
   std::cout << m;
   //m.mset("1 0 2 ; 0-1 3 1"); n.mset("3 ; 2 ; 1");
   //std::cout << n;
   //m*=n; std::cout << m;
   m.mset("1 0 2 ; 0-1 3 1"); std::cout << "Full display gives:\n" << m.display_full();
   m(0,1) = 11; std::cout << "m(0,1)=11 gives " << m(0,1) << " gives:\n" << m.display_full();
   v.push_back(3); v.push_back(2); v.push_back(1);
   v = m*v; std::cout << "V size is " << v.size() << ". V = [" << v[0] << "," << v[1] << "]\n";
   m*=5; std::cout << "m=\n" << m.display_full(); 
   std::cout << "n=\n" << n.display_full();
   n.transpose(); a = m % n; std::cout << "a = m % n' =\n" << a.display_full();
   a = m ^ n; std::cout << "a = m ^ n' =\n" << a.display_full();
   
   // Tests numdenom
   //numdenom nd,nd3;
   //std::cout << "before: " << nd << " after: ";
   //nd = 1; std::cout << "nd = 1 is " << nd << "\n";
   //std::cout << "(double)nd gives "<< (double)nd << "\n";
   //numdenom nd2(1,2,2,1); std::cout << "nd = sqrt(2)/2 is " << nd2 << "\n";
   //std::cout << "double(nd) gives "<< double(nd2) << "\n";
   //nd = 2; std::cout << "nd = 2 is " << nd << "\n";
   //nd3 = sqrt(nd);
   //std::cout << "nd = sqrt(2) gives "<< nd3 << " = " << double(nd3) << "\n";

   std::cout << "#rows in m is: " << m.nr() << "\n";

   // Tests pset
   a.mset("6 0-2 5 7 ; 8 0-7 4 2 ; 1 0-3 1 5 ; 2 5 3 1"); b.mset("6 0-3 ; 2 9");
   std::cout << "a =\n" << a.display_full();
   std::cout << "b =\n" << b.display_full();
   a.pset(2,3,2,3,b); std::cout << "a.pset(2,3,2,3,b) gives:\n" << a.display_full();

   // Tests diagonalisation routines
   c = mabs(a); std::cout << "mabs(a) =\n" << c.display_full();
   v = mmax(a); std::cout << "mmax(a) = " << dispvect(v) << "\t vmax(mmax(a)) = " << vmax(v) << "\n";
   v = msum(a); std::cout << "msum(a) = " << dispvect(v) << "\t vsum(msum(a)) = " << vsum(v) << "\n";
   c = a; c.transpose(); std::cout << "a.transpose() =\n" << c.display_full();
   std::cout << "a.issymm() = " << a.issymm() << "\n";
   b = triu(a); c = tril(c); c = b+c; std::cout << "c = triu(a)+ tril(a.transpose) =\n" << c.display_full();
   std::cout << "c.issymm() = " << c.issymm() << "\n";
   // Matrix from the Handbook of Automatic Computation, pp.223
   a.mset("10 1 2 3 4 ; 1 9 0-1 2 0-3 ; 2 0-1 7 3 0-5 ; 3 2 3 12 0-1 ; 4 0-3 0-5 0-1 15");
   b = triu(a,1); std::cout << "triu(a,1) =\n" << b.display_full();
   b.mset("5 1 0-2 0 0-2 5 ; 1 6 0-3 2 0 6 ; 0-2 0-3 8 0-5 0-6 0 ; 0 2 0-5 5 1 0-2 ; 0-2 0 0-6 1 6 0-3 ; 5 6 0 0-2 0-3 8");
   std::cout << "a =\n" << a.display_full();
   dg = eig(a); std::cout << "dg=eig(a); dg.E = " << dispvect(dg.E) << "\n" << "dg.V =\n" << dg.V.display_full();
   std::cout << "b =\n" << b.display_full();
   dg = eig(b); std::cout << "dg=eig(b); dg.E = " << dispvect(dg.E) << "\n" << "dg.V =\n" << dg.V.display_full();
   v = svd(a); std::cout << "svd(a) = " << dispvect(v) << "\n";
   v = svd(b); std::cout << "svd(b) = " << dispvect(v) << "\n";
   a.mset("1 1 9 ; 5 8 5"); std::cout << "a is now\n" << a.display_full();
   v2 = svd(a); std::cout << "svd(a) = " << dispvect(v) << "\n";
   a.mset("5 2 ; 4 4 ; 5 6"); std::cout << "a is now\n" << a.display_full();
   v = svd(a,m,n); std::cout << "svd(a,u,v) = " << dispvect(v) << "\nU =\n" << m.display_full() << "V =\n" << n.display_full();
   a.transpose(); std::cout << "a is now\n" << a.display_full();
   v = svd(a,m); std::cout << "svd(a,u) = " << dispvect(v2) << "\nU =\n" << m.display_full();
   std::cout << "norm(a) = " << norm(a) << "\n";
   std::cout << "rank(a) = " << rank(a,(double)vmax<int>(a.size())*vmax<double>(v)*DBL_EPSILON) << "\n";
   m = orth(a); std::cout << "orth(a) =\n" << m.display_full();
   a.transpose(); std::cout << "a is now\n" << a.display_full();
   m = orth(a); std::cout << "orth(a) =\n" << m.display_full();
   b.mset("7 ; 8 ; 9"); a.pset(1,a.nr(),a.nc()+1,a.nc()+1,b); std::cout << "a is now\n" << a.display_full();
   a = randn(2,3); std::cout << "randn(2,3) =\n" << a.display_full(); 
   v.resize(5,0.); v[0]=-1; v[1]=0; v[2]=2; v[3]=4; v[4]=6; v2.resize(4,0.); v2[0]=-1; v2[1]=0; v2[2]=1; v2[3]=3;
   std::cout << "setunion(" << dispvect(v) << "," << dispvect(v2) << ") = " << dispvect(setunion(v,v2)) << "\n";
   std::cout << "setdiff(" << dispvect(v) << "," << dispvect(v2) << ") = " << dispvect(setdiff(v,v2)) << "\n";
   std::cout << "setxor(" << dispvect(v) << "," << dispvect(v2) << ") = " << dispvect(setxor(v,v2)) << "\n";
   a((int)(rand()%2),(int)(rand()%3)) = 0; std::cout << "a is now\n" << a.display_full();
   v = all(a); std::cout << "all(a) = " << dispvect(v) << "\n";
   v = all(a,2); std::cout << "all(a,2) = " << dispvect(v) << "\n";
   b = qr(a,c); std::cout << "qr(a,c) =\n" << b.display_full() << "c =\n" << c.display_full();
   //a.mset("5 1 0-2 0 0-2 5 ; 1 6 0-3 2 0 6 ; 0-2 0-3 8 0-5 0-6 1 ; 0 2 0-5 5 1 0-2 ; 0-2 0 0-6 1 6 0-3 ; 5 6 0 0-2 0-3 8");
   a.mset("1 2 3 ; 1 5 6 ; 1 8 9 ; 1 11 12");
   std::cout << "a is now\n" << a.display_full();
   b = qr(a,c); std::cout << "qr(a,c) =\n" << b.display_full() << "c =\n" << c.display_full();
   a.mset("1 2 0 0 ; 2 -3 5 0 ; 0 5 4 -1 ; 0 0 -1 2");
   std::cout << "a is now\n" << a.display_full();
   dg = eigT(a); std::cout << "dg=eigT(a); dg.E = " << dispvect(dg.E) << "\n dg.V =\n" <<dg.V.display_full();
}*/
