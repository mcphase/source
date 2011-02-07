#include<myev.h>

#define VERYSMALL 1e-04
#ifndef MAXNOFCHARINLINE
#define MAXNOFCHARINLINE 1000
#endif

void checkHerm(ComplexMatrix & M)
{double d,max=0;
 int i1,j1,n=0;
 for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=M.Chi();++j1){
   d=abs(M(i1,j1)-conj(M(j1,i1)));
   if (d>VERYSMALL){n++;if(d>max)max=d;}
   }}
  if(n>0)
   {fprintf(stderr,"myEigenSystemHermitean: ERROR- %ix%i matrix not hermitian\n %i offdiagonal elements deviate by more than %g\n - largest deviation of offdiagonal elements: %g\n\n Press enter to ignore and continue, press p to print matrix\n?\n",M.Rhi(),M.Rhi(),n,VERYSMALL,max);
    //printout matrix
    if(getchar()=='p'){getchar();myPrintComplexMatrix(stderr,M);fprintf(stderr,"press enter to continue\n");getchar();}
   }

}

// subs to be able to check and directly diagonalize hermitean
// matrizes, inverse a nearly singular matrix
void myPrintComplexMatrix(FILE * file,ComplexMatrix & M)
{int i1,j1;
 fprintf (file,"#Real Part\n");
   
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
    for (j1=M.Clo();j1<=M.Chi();++j1) fprintf (file,"%6.3g ",real(M(i1,j1)));
    fprintf (file,"\n");
    }
    fprintf (file,"#Imaginary Part\n");
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
   for (j1=M.Clo();j1<=M.Chi();++j1)fprintf (file,"%6.3g ",imag(M(i1,j1)));
    fprintf (file,"\n");
    }
}    

int myReadComplexMatrix (FILE * file, ComplexMatrix & M)
{int i1,j1;char instr[MAXNOFCHARINLINE];
 float *numbers;
  numbers = new float[M.Rhi()-M.Rlo()+2];
  numbers[0]=M.Rhi()-M.Rlo()+2;
 
     //read comment line 
     if(fgets(instr,MAXNOFCHARINLINE,file)==false) {fprintf (stderr, "ERROR reading complex matrix - comment line before real part\n");return false;}

    // read real part
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
    j1=inputline(file,numbers);if(j1!=M.Chi()-M.Clo()+1) {fprintf (stderr, "ERROR reading complex matrix - number of columns read (%i) does not match matrix dimension (%i)\n",j1,M.Chi()-M.Clo()+1);return false;}
    for (j1=M.Clo();j1<=M.Chi();++j1)M(i1,j1)=complex <double> (numbers[j1-M.Clo()+1],0);
    }

     //read comment line 
     if(fgets(instr,MAXNOFCHARINLINE,file)==false) {fprintf (stderr, "ERROR reading complex matrix - comment line before imaginary part\n");return false;}
    // read imaginary part
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
    j1=inputline(file,numbers);if(j1!=M.Chi()-M.Clo()+1) {fprintf (stderr, "ERROR reading complex matrix - number of columns read (%i) does not match matrix dimension (%i)\n",j1,M.Chi()-M.Clo()+1);return false;}
    for (j1=M.Clo();j1<=M.Chi();++j1)M(i1,j1)+=complex <double> (0,numbers[j1-M.Clo()+1]);
    }
delete []numbers;
return true;
}

void myPrintMatrix(FILE * file,Matrix & M)
{int i1,j1;
   for (i1=M.Rlo();i1<=M.Rhi();++i1){
    for (j1=M.Clo();j1<=M.Chi();++j1) fprintf (file,"%6.3g ",M(i1,j1));
    fprintf (file,"\n");
    }
}    

void myPrintVector(FILE * file,Vector & M)
{int j1;
 fprintf (file,"#Components:\n");
   
    for (j1=M.Lo();j1<=M.Hi();++j1) fprintf (file,"%6.3g ",M(j1));
    fprintf (file,"\n");    
}    

void myPrintComplexVector(FILE * file,ComplexVector & M)
{int j1;
 fprintf (file,"#Components:\n");
   
    for (j1=M.Lo();j1<=M.Hi();++j1) fprintf (file,"%6.3g %+6.3g i ",real(M(j1)),imag(M(j1)));
    fprintf (file,"\n");    
}    


void myEigenValuesHermitean (ComplexMatrix & M,Vector & lambda,int & sort,int & maxiter)
{ // this sub diagonalizes M and puts eigenvalues to lambda
  Matrix mat1(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
  int i1,j1;
  checkHerm(M);

  // put matrix to format needed for library diagonalize function
   for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=M.Chi();++j1){
    mat1(j1,i1)=imag(M(i1,j1)); 
    mat1(i1,j1)=real(M(i1,j1));
   }}
   EigenValuesHermitean (mat1,lambda,sort,maxiter);
   return;

}

void myEigenSystemHermitean (ComplexMatrix & M,Vector & lambda,ComplexMatrix & l,int & sort,int & maxiter)
{ // this sub diagonalizes M and puts eigenvalues to lambda, the eigenvectors to l
  Matrix mat1(M.Rlo(),M.Rhi(),M.Clo(),M.Chi());
  Matrix lr(l.Rlo(),l.Rhi(),l.Clo(),l.Chi());
  Matrix li(l.Rlo(),l.Rhi(),l.Clo(),l.Chi());
  complex<double> ii(0,1);  
  int i1,j1;
  //check if M it is hermitean
   checkHerm(M);

  // put matrix to format needed for library diagonalize function
   for(i1=M.Rlo();i1<=M.Rhi();++i1){for(j1=M.Clo();j1<=M.Chi();++j1){
    mat1(j1,i1)=imag(M(i1,j1)); 
    mat1(i1,j1)=real(M(i1,j1));
   }}
 // setup matrix and diagonalize
 //  Driver routine to compute the  eigenvalues and normalized eigenvectors 
 //  of a complex Hermitian matrix z.The real parts of the elements must be
 //  stored in the lower triangle of z,the imaginary parts (of the elements
 //  corresponding to the lower triangle) in the positions
 //  of the upper triangle of z[lo..hi,lo..hi].The eigenvalues are returned
 //  in d[lo..hi] in ascending numerical  order if the sort flag is set  to
 //  True, otherwise  not ordered for sort = False. The real  and imaginary
 //  parts of the eigenvectors are  returned in  the columns of  zr and zi. 
 //  The storage requirement is 3*n*n + 4*n complex numbers. 
 //  All matrices and vectors have to be allocated and removed by the user.
 //  They are checked for conformance !
 // void  EigenSystemHermitean (Matrix& z, Vector& d, Matrix& zr, Matrix& zi, 
 // 			   int sort, int maxiter)
   EigenSystemHermitean (mat1,lambda,lr,li,sort,maxiter);
  // l=li;l*=ii;l+=lr;
  l=ComplexMatrix(li,lr);
   return;
}


void myEigenSystemHermiteanGeneral (ComplexMatrix& a, ComplexMatrix& b, Vector & e, ComplexMatrix & T, int & sort, int & maxiter)
{ // this sub diagonalizes M and puts eigenvalues to lambda, the eigenvectors to l
  Matrix mata(a.Rlo(),a.Rhi(),a.Clo(),a.Chi());
  Matrix matb(b.Rlo(),b.Rhi(),b.Clo(),b.Chi());
  Matrix zr(T.Rlo(),T.Rhi(),T.Clo(),T.Chi());
  Matrix zi(T.Rlo(),T.Rhi(),T.Clo(),T.Chi());
  ComplexVector x(T.Rlo(),T.Rhi());
  complex<double> ii(0,1);  
  int i1,j1;

  //check if a,b it is hermitean
  checkHerm(a);
  checkHerm(b);


  // put matrix to format needed for library diagonalize function
   for(i1=a.Rlo();i1<=a.Rhi();++i1){for(j1=a.Clo();j1<=a.Chi();++j1){
    mata(j1,i1)=imag(a(i1,j1)); 
    mata(i1,j1)=real(a(i1,j1));
   }}

   for(i1=b.Rlo();i1<=b.Rhi();++i1){for(j1=b.Clo();j1<=b.Chi();++j1){
    matb(j1,i1)=imag(b(i1,j1)); 
    matb(i1,j1)=real(b(i1,j1));
   }}

//
//  Driver routine for the generalized hermitan eigenvalue problem:
//
//                     A * z = e * B * z
//
//  The real part of the  complex  hermitean matrices a[lo..hi,lo..hi] and
//  b[lo..hi,lo..hi]  must be stored in  the lower triangle, the imaginary
//  parts must be stored in the strict upper triangle. The eigenvalues are
//  returned in e[lo..hi] in ascending numerical order if the sort flag is 
//  set to  True,  otherwise not  ordered  for sort = False. The real  and 
//  imaginary parts of the eigenvectors are returned in  the columns of zr
//  and zi. C.f. comments on Chreduce(), Chreback() and the 
//  EigenSystemHermitean()  routine for the complex Hermitean problem 
//  (file cheigen.c)
//  All matrices and vectors have to be allocated and removed by the user.
//  They are checked for conformance !
//  
//  References:
//
//  B.T.Smith et al: Matrix Eigensystem Routines
//  EISPACK Guide,Springer,Heidelberg,New York 1976.
EigenSystemHermiteanGeneral (mata, matb, e,zr, zi,sort, maxiter);
  T=ComplexMatrix(zi,zr);
// normalize eigenvectors ( this is not automatically done);
  for(i1=T.Clo();i1<=T.Chi();++i1)
   {x=T.Column(i1);
    x=x/Norm(x);
    for(j1=T.Rlo();j1<=T.Rhi();++j1) {T(j1,i1)=x(j1);}
   }

   return;

}
