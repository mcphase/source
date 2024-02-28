#define EPS DBL_EPSILON

void jjjpar::cluster_ini_Imat() // to be called on initializing the cluster module
{// vector of matrix pointers which index the I1,I2,...In matrices of the cluster
 zsMat<double> * Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+3*(*clusterpars).cs.nofatoms+1];
 // initialize these matrices
 dim=1; Vector Hxc(1,(*clusterpars).cs.nofcomponents);Vector Hext(1,3);
 dnn= new int [(*clusterpars).cs.nofatoms+1];
 // determine dimension of H matrix
 for (int n=1;n<=(*clusterpars).cs.nofatoms;++n)
 {dnn[n]=(*(*clusterpars).jjj[n]).opmat(1,Hxc,Hext).Rhi();
  dim*=dnn[n];
 }
 printf("# cluster module - number of matrices: %i  dimension=%i  \n",(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+3*(*clusterpars).cs.nofatoms,dim);

 // Before initializing (allocating) matrices for other operators, check to see if we are using own parser, and if
 //   so get the sequence of operations and perform them here directly for each operator without allocating all the
 //   matrices  -  MDL 131024
 useperl=true; truncate=0; feast=0; arpack=0; oldeig=false; FILE *fin=fopen(sipffilename,"rb"); char instr[MAXNOFCHARINLINE];
 while(feof(fin)==0) if(fgets(instr,MAXNOFCHARINLINE,fin)!=NULL) {
    if(strncmp(instr,"#!noperl",8)==0) useperl=false; 
    else if(strncmp(instr,"#!eigsolver",11)==0) {      // Parses eigensolver string - can be either: truncate, arpack or feast
       char chomp[MAXNOFCHARINLINE], ts=0, fs=0, as=0; // New string without whitespaces
       char *ii=instr, *ic=chomp; while(*ii!=0) { if(*ii!=11||*ii!=' ') { *ic=tolower(*ii); ic++; } ii++; } *ic=0;
       ic=strchr(chomp,'=');
       do 
       {
          ic++;
          if(strncmp(ic,"truncate",8)==0) { 
             ii=strchr(ic,'='); if(ii==NULL) truncate=DBL_MAX; else { ii++; if(*ii=='e') truncate=DBL_EPSILON; else if(*ii=='f'||*ii=='a') ts=*ii; else truncate=atof(ii); } }
          else if(strncmp(ic,"feast",5)==0) {
             ii=strchr(ic,'='); if(ii==NULL) feast=DBL_MAX;    else { ii++; if(*ii=='e') feast=DBL_EPSILON;    else if(*ii=='t'||*ii=='a') fs=*ii; else feast=atof(ii); } }
          else if(strncmp(ic,"arpack",6)==0) {
             ii=strchr(ic,'='); if(ii==NULL) arpack=DBL_MAX;   else { ii++; if(*ii=='e') arpack=DBL_EPSILON;   else if(*ii=='t'||*ii=='f') as=*ii; else arpack=atof(ii); } }
          else if(strncmp(ic,"oldeig",6)==0||strncmp(ic,"matpack",7)==0) { oldeig = true; }
       }
       while( (ic=strchr(ic,','))!=NULL );
       if(ts!=0) { if(ts=='f') truncate=feast;  else truncate=arpack; }
       if(fs!=0) { if(fs=='t') feast=truncate;  else feast=arpack; }
       if(as!=0) { if(as=='t') arpack=truncate; else arpack=feast; }
       if(feast!=0 && arpack!=0) { fprintf(stderr,"WARNING: cluster_module cannot use both FEAST and ARPACK together. FEAST will be used.\n"); arpack=0; }
       else if(arpack!=0) { 
          fprintf(stderr,"\n -----------\n WARNING: ARPACK eigensolver selected. Note that this option is buggy. Proceed at you own risk!\n ----------- \n"); fflush(stderr); }
       if(arpack==DBL_MAX) arpack = 0.1; if(truncate==DBL_MAX) truncate = 0.1;
       printf("truncate=%g\tfeast=%g\tarpack=%g\n",truncate,feast,arpack); /*truncate=0;*/ feast=0; /*arpack=0;*/
    }
  //else if(strncmp(instr,"#!truncate",10)==0) { char *valsto = strchr(instr,'=')+1; truncate = atof(valsto); }
 }

 // initialize matrix and vector to cache Hamiltonian matrix from runs where Hext is constant
 clusterH = new zsMat<double>(dim,dim); oldHext = new Vector(1,3); justinit=true;
 // Allocates workspace for iterative eigensolvers
 workspace = new iterwork(dim,dim,dim);

 Iaa[0]=new zsMat<double>(dim,dim); (*Iaa[0])=1.; //Iaa[0]=new ComplexMatrix(1,dim,1,dim);(*Iaa[0])=1; 
 for(int a=1;a<=(*clusterpars).cs.nofatoms;++a)
   for(int i=1;i<=(*clusterpars).cs.nofcomponents;++i) {
     Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i] = new zsMat<double>(dim,dim);//ComplexMatrix(1,dim,1,dim); 
     cluster_Iaa(Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i],a,i); }

  // her we should make possible the perlparsing the cluster sipf file
  // in order to flexibly change the Ia operator sequence ...
  double *numbers[1]; numbers[0]=&gJ;
  char numnam[]="gJ\0";
  char * numbernames[2];numbernames[0]=numnam;numbernames[1]=NULL;

  char * strings[1];char * stringnames[1];
  stringnames[0]=NULL;

 char * operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+3*(*clusterpars).cs.nofatoms+2];
 // char opnam[]="one\0";   
  operatornames[0]=new char[5];strcpy(operatornames[0],"one\0");
 for(int a=1;a<=(*clusterpars).cs.nofatoms;++a)
 for(int i=1;i<=(*clusterpars).cs.nofcomponents;++i)
 {operatornames[(a-1)*(*clusterpars).cs.nofcomponents+i]=new char[10];
  snprintf(operatornames[(a-1)*(*clusterpars).cs.nofcomponents+i],MAXNOFCHARINLINE,"I%i_%i",a,i);
 }     

 int **cluster_seq, *cluster_seqconstInd, *cluster_opassignflags, cluster_nlines=0; double **cluster_seqconst; char **cluster_statements;
 int cluster_Ia_ind0, cluster_M_ind0;
 if(!useperl) { 
    cluster_nlines=0;
    // Allocates arrays to store the sequences from the non-perl parsing function
    cluster_seq = new int*[999]; cluster_seqconstInd = new int[999]; cluster_opassignflags = new int[999]; 
    cluster_seqconst = new double*[999]; cluster_statements = new char*[999];
    // Define all the operator names
    int Ia0=(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+1; cluster_M_ind0=Ia0+nofcomponents;
    for(int i=1;i<=nofcomponents;++i) { 
       operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+i]=new char[5];
       snprintf(operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+i],MAXNOFCHARINLINE,"I%i",i); }
    for(int i=1;i<=3;++i) { 
       operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+i]=new char[5];
       snprintf(operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+i],MAXNOFCHARINLINE,"M%i",i); }
    int index=0; cluster_Ia_ind0=(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+1;
    for(int a=1;a<=(*clusterpars).cs.nofatoms;++a) for(int n=1;n<=3;++n) { index=(a-1)*3+n;
       operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index]=new char[5];
       snprintf(operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index],MAXNOFCHARINLINE,"M%i_%i",a,n); }
    int nop = (*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index;
    operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+3*(*clusterpars).cs.nofatoms+1]=NULL;

    // Call myparse to get the sequence of operations for each statement
    cluster_nlines = myparse(sipffilename,numbers,numbernames,strings,stringnames,Iaa,operatornames,
                     cluster_seq,cluster_seqconst,cluster_statements,cluster_seqconstInd,cluster_opassignflags); 
    if(cluster_nlines<1) {
       printf("Error own parsing sipf file %s\n",sipffilename); exit(EXIT_FAILURE); }

    // Classify left (lhs) and right (rhs) hand side variables in each equations
    int lhs[999], rhs[999], isq; memset(lhs,0,999*sizeof(int)); memset(rhs,0,999*sizeof(int));
    for(int i=0; i<cluster_nlines; i++) { 
       lhs[cluster_seq[i][0]] = 1; 
       isq = 1; while(cluster_seq[i][isq]!=0) { if(cluster_seq[i][isq]<0) rhs[-cluster_seq[i][isq]] = 1; isq++; }
       if(cluster_seq[i][1]>2) rhs[cluster_seq[i][0]] = 1; // For assignment operators, += -= *= need to save matrix too.
    }

    // If any operator other than In_m are required on the right hand sides, we need to allocate them now before performing the operations
    for(int i=cluster_Ia_ind0; i<=nop; i++) if(rhs[i]==1) Iaa[i]=new zsMat<double>(dim,dim);//ComplexMatrix(1,dim,1,dim);
    // Allocates the interaction and cluster moment matrices
    Ia= new zsMat<double> * [nofcomponents+1];
    cluster_M= new zsMat<double> * [3+3*(*clusterpars).cs.nofatoms+1];
    for(int n = 1;n<=nofcomponents;++n) Ia[n]=new zsMat<double>(dim,dim);
    for(int n = 1;n<=3;++n) cluster_M[n]=new zsMat<double>(dim,dim);
    for(int a=1;a<=(*clusterpars).cs.nofatoms;++a) for(int n=1;n<=3;++n) {
       int index_M=3*a+n; cluster_M[index_M]=new zsMat<double>(dim,dim); }

    // Runs through the list of operators to be defined
    zsMat<double> *dum = new zsMat<double>(dim,dim);
    for(int iln=0; iln<cluster_nlines; iln++) {
       *dum = 0;
       // If matrix is not required on RHS, it is unallocated and _should_ be zero - set it to a dummy (zero) matrix's address
       if(rhs[cluster_seq[iln][0]]==0) Iaa[cluster_seq[iln][0]] = dum;
       if(myparse_execute(Iaa, operatornames, cluster_seq[iln], cluster_seqconst[iln], 
                          cluster_statements[iln], cluster_seqconstInd[iln], cluster_opassignflags[iln])==false) exit(EXIT_FAILURE);
       if(cluster_seq[iln][0]>=cluster_M_ind0) {
        //for(int i=1;i<=dim;++i) for(int j=1;j<=dim;++j) {
//           if(i<j) {
//              const complex<double> tva = (*Iaa[cluster_seq[iln][0]])(j,i);
//              (*cluster_M[cluster_seq[iln][0]-cluster_M_ind0+1])(i,j) = imag(tva); }
//           else    { 
//              const complex<double> tva = (*Iaa[cluster_seq[iln][0]])(i,j);
//              (*cluster_M[cluster_seq[iln][0]-cluster_M_ind0+1])(i,j) = real(tvz); }
        //   if(i<j) { (*cluster_M[cluster_seq[iln][0]-cluster_M_ind0+1])(i,j) = imag((*Iaa[cluster_seq[iln][0]])[make_pair(j,i)]); }
        //   else    { (*cluster_M[cluster_seq[iln][0]-cluster_M_ind0+1])(i,j) = real((*Iaa[cluster_seq[iln][0]])[make_pair(i,j)]); }
        //}
          *cluster_M[cluster_seq[iln][0]-cluster_M_ind0+1] = *Iaa[cluster_seq[iln][0]];
       }
       else if(cluster_seq[iln][0]>=cluster_Ia_ind0) {
        //for(int i=1;i<=dim;++i) for(int j=1;j<=dim;++j) {
        //   if(i<j) { (*Ia[cluster_seq[iln][0]-cluster_Ia_ind0+1])(i,j) = imag((*Iaa[cluster_seq[iln][0]])[make_pair(j,i)]); }
        //   else    { (*Ia[cluster_seq[iln][0]-cluster_Ia_ind0+1])(i,j) = real((*Iaa[cluster_seq[iln][0]])[make_pair(i,j)]); }
        //}
          *Ia[cluster_seq[iln][0]-cluster_Ia_ind0+1] = *Iaa[cluster_seq[iln][0]]; 
       }
       if(rhs[cluster_seq[iln][0]]==0) Iaa[cluster_seq[iln][0]] = NULL;
    }
    delete dum;

    if(dim<200) {
    char ploutname[MAXNOFCHARINLINE];  snprintf(ploutname,MAXNOFCHARINLINE,"results/_%s.pl.out",sipffilename);
    FILE *PLOUT = fopen(ploutname,"w");
    for(int a=0; a<=nop; a++)
    {
       if(a<cluster_Ia_ind0) dum = Iaa[a]; else if(a<cluster_M_ind0) dum = Ia[a-cluster_Ia_ind0+1]; else dum = cluster_M[a-cluster_M_ind0+1]; 
       fprintf(PLOUT,"%s=\n#Real Part\n",operatornames[a]);
       for(int i=1; i<=dim; i++) { fprintf(PLOUT,"%8.5f",real((*dum)[make_pair(i,1)])); for(int j=2; j<=dim; j++) fprintf(PLOUT," %8.5f",real((*dum)[make_pair(i,j)])); fprintf(PLOUT,"\n"); }
       fprintf(PLOUT,"#Imaginary Part\n");
       for(int i=1; i<=dim; i++) { fprintf(PLOUT,"%+8.5f",imag((*dum)[make_pair(i,1)])); for(int j=2; j<=dim; j++) fprintf(PLOUT," %+8.5f",imag((*dum)[make_pair(i,j)])); fprintf(PLOUT,"\n"); }
    }
    fclose(PLOUT);
    }
    for(int i=0; i<cluster_Ia_ind0; i++) delete Iaa[i];
    for(int i=cluster_Ia_ind0; i<=nop; i++) if(rhs[i]==1) delete Iaa[i];

    // Delete the operatornames arrays
    delete[]operatornames[0];
    for(int a=1;a<=(*clusterpars).cs.nofatoms;++a) for(int i=1;i<=(*clusterpars).cs.nofcomponents;++i) delete[]operatornames[(a-1)*(*clusterpars).cs.nofcomponents+i];
    for(int i=1;i<=nofcomponents;++i) delete[]operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+i];
    for(int i=1;i<=3;++i) delete[]operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+i];
    for(int a=1;a<=(*clusterpars).cs.nofatoms;++a) for(int n=1;n<=3;++n) { int index=(a-1)*3+n; 
       delete[]operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index]; }
    // Delete the sequence arrays
    for(int i=0; i<cluster_nlines; i++) {
       delete[]cluster_seq[i]; delete[]cluster_seqconst[i]; delete[]cluster_statements[i]; }
    delete[]cluster_seqconstInd; delete[]cluster_opassignflags; delete[]cluster_seq; delete[]cluster_seqconst; delete[]cluster_statements;
 }
 else {

 // now initialize the interaction operators
 for(int i=1;i<=nofcomponents;++i)
 {operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+i]=new char[5];
  snprintf(operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+i],MAXNOFCHARINLINE,"I%i",i);
  Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+i]=new zsMat<double>(dim,dim);//ComplexMatrix(1,dim,1,dim);                     
 }
 // now initialize the cluster total magnetic moment operators
 for(int i=1;i<=3;++i) 
 {operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+i]=new char[5];
  snprintf(operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+i],MAXNOFCHARINLINE,"M%i",i);
  Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+ i]=new zsMat<double>(dim,dim);//ComplexMatrix(1,dim,1,dim);                     
 }

 // now initialize the cluster individual magnetic moment operators 
 //(for neutron form factor in dipole approx - function MQ dMQ1)
 for(int a=1;a<=(*clusterpars).cs.nofatoms;++a)for(int n=1;n<=3;++n)
 {int index=(a-1)*3+n;
  operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index]=new char[5];
  snprintf(operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index],MAXNOFCHARINLINE,"M%i_%i",a,n);
  Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index]=new zsMat<double>(dim,dim);//ComplexMatrix(1,dim,1,dim);                     
 } 
  operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+3*(*clusterpars).cs.nofatoms+1]=NULL;

 if(perlparse(sipffilename,numbers,numbernames,strings,stringnames,Iaa,operatornames)==false)
   {printf("Error perl parsing sipf file %s\n",sipffilename);exit(EXIT_FAILURE);}
for(int i=0;i<=(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents;++i)
  {delete Iaa[i];delete operatornames[i];}

// here we fill the interaction operator matrices with values
Ia= new zsMat<double> * [nofcomponents+1];
for(int n = 1;n<=nofcomponents;++n){Ia[n]=new zsMat<double>(dim,dim);
//for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){
//   if(i<j){(*Ia[n])(i,j)=imag((*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+n])[make_pair(j,i)]);}
//      else{(*Ia[n])(i,j)=real((*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+n])[make_pair(i,j)]);}
//}
(*Ia[n])=(*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+n]);
delete           Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+n];
delete operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+n]; 
}

// here we should initialize set the operator matrices cluster_M  ... total magnetic moment
// operator storage and individual magnetic moments: this should then be used
// in the mcalc, dm1calc, MQ and dMQ1 functions of the cluster module ...
// to be done !!!

cluster_M= new zsMat<double> * [3+3*(*clusterpars).cs.nofatoms+1];
// first do the total moment
for(int n = 1;n<=3;++n){cluster_M[n]=new zsMat<double>(dim,dim);
//for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){
// //if(i<j){(*cluster_M[n])(i,j)=imag((*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+n])[make_pair(j,i)]);}
// //   else{(*cluster_M[n])(i,j)=real((*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+n])[make_pair(i,j)]);}
//   (*cluster_M[n])(i,j) = (*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+n])[make_pair(i,j)];
//}
(*cluster_M[n]) = (*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+n]);
delete           Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+n];
delete operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+n]; 
}
// now do the individual atom magnetic moments
 for(int a=1;a<=(*clusterpars).cs.nofatoms;++a)for(int n=1;n<=3;++n)
{int index=(a-1)*3+n; // a .... atom index  n ... xyz components of magnetic moment
 int index_M=3*a+n;
 cluster_M[index_M]=new zsMat<double>(dim,dim);
 (*cluster_M[index_M]) = (*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index]);
 //for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){
 //if(i<j){(*cluster_M[index_M])(i,j)=imag((*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index])[make_pair(j,i)]);}
 //   else{(*cluster_M[index_M])(i,j)=real((*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index])[make_pair(i,j)]);}
 //  (*cluster_M[index_M])(i,j) = (*Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index])[make_pair(i,j)];
 //}
delete           Iaa[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index];
delete operatornames[(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+index]; 
}

}
fdim=dim; 
if (truncate>1e-6 && truncate<1) {
   dim = (int)(ceil(truncate*(double)dim)); is1sttrunc = true; zm = new complexdouble[fdim*fdim]; }
// for(int i=0;i<=(*clusterpars).cs.nofatoms*(*clusterpars).cs.nofcomponents+nofcomponents+3+3*(*clusterpars).cs.nofatoms;++i)
//  {delete Iaa[i];delete operatornames[i];}
printf("#module cluster initialized\n");
}

//------------------------------------------------------------------------------------------------
// Routine to calculate the Iaa (Ia_i) matrix for some particular a,i
//------------------------------------------------------------------------------------------------
void jjjpar::cluster_Iaa(zsMat<double> *Iai, int a, int i)
{
 Vector Hxc(1,(*clusterpars).cs.nofcomponents);Vector Hext(1,3);
 // initialize matrices
//Iaa[0]=new ComplexMatrix(1,dim,1,dim);(*Iaa[0])=1;
//for(int i=1;i<=(*clusterpars).cs.nofcomponents;++i)
  { //Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i]=new ComplexMatrix(1,dim,1,dim);
    //(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])=0;
    (*Iai)=0;
    Matrix Jai((*(*clusterpars).jjj[a]).opmat(i,Hxc,Hext));
    // myPrintMatrix(stdout,Jai);printf("\n");

  //1. determine dimensions
  int da=dnn[a];
  int dx=1,dz=1;
  for(int m=1  ;m<a;++m)dx*=dnn[m];
  for(int m=a+1;m<=(*clusterpars).cs.nofatoms;
                    ++m)dz*=dnn[m];
/*
  int mx=0,ma=1;
  if(a>1)mx=1;
  ma = dz; mx = dz*dnn[a];
//for(int m=a;  m<=(*clusterpars).cs.nofatoms;++m)mx*=dnn[m];
//for(int m=a+1;m<=(*clusterpars).cs.nofatoms;++m)ma*=dnn[m];
  //2. insertion loop
  if(useperl) {
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=da-1;++ai)
  for(int bi=0;bi<=da-1;++bi)
  {int r=ix*mx+ai*ma+iz+1;
   int s=ix*mx+bi*ma+iz+1;
 //printf("hello %i %i %i %i\n",r,s,ai,bi);
  // here we should fill the matrices with values corresponding to the
  // I1 I2 of the individual atoms of the cluster
// if(r<s)        {(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])(r,s)+=complex<double>(Jai(bi+1,ai+1),-Jai(ai+1,bi+1));}
// else if (r==s) {(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])(r,s)+=complex<double>(Jai(ai+1,bi+1),0.0);}   
//            else{(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])(r,s)+=complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1));}   
   if(r<s)        {(*Iai)(r,s)+=complex<double>(Jai(bi+1,ai+1),-Jai(ai+1,bi+1));}
   else if (r==s) {(*Iai)(r,s)+=complex<double>(Jai(ai+1,bi+1),0.0);}   
              else{(*Iai)(r,s)+=complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1));}   
  }
  } else { 
*/
//Matrix outmat(1,dim,1,dim); outmat=0;
/*  for(int ix=0;ix<=dx-1;++ix) for(int iz=0;iz<=dz-1;++iz)
       for(int ai=0;ai<=da-1;++ai) for(int bi=0;bi<=ai;++bi)
       { 
          int r=ix*mx+ai*ma+iz+1; int s=ix*mx+bi*ma+iz+1;
        //if(r<s)        {(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])(r,s)+=complex<double>(Jai(bi+1,ai+1),-Jai(ai+1,bi+1));}
//             if (r==s) {(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])(r,s)+=complex<double>(Jai(ai+1,bi+1),0.0);}   
//        else if(r>s)   {(*Iaa[(a-1)*(*clusterpars).cs.nofcomponents+i])(r,s)+=complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1));}
               if (r==s) {(*Iai)(r,s)+=complex<double>(Jai(ai+1,bi+1),0.0);}   
          else if(r>s)   {(*Iai)(r,s)+=complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1));}
       }*/
/*  for(int ix=0;ix<=dx-1;++ix) for(int iz=0;iz<=dz-1;++iz)
       for(int ai=0;ai<=da-1;++ai) for(int bi=0;bi<=da-1;++bi) 
       {
          int r=ix*mx+ai*ma+iz+1; int s=ix*mx+bi*ma+iz+1; outmat(r,s) = Jai(ai+1,bi+1); 
       }
myPrintMatrix(stdout,outmat);printf("\n");*/
   int mx = dz*dnn[a];
   for(int ix=0; ix<dx; ix++)
    for(int ai=0; ai<da; ai++) { int r=ix*mx+ai*dz+1;
     for(int bi=0; bi<=ai; bi++) { int s=ix*mx+bi*dz+1;
      for(int iz=0; iz<dz; iz++) {
//             if (r==s) {(*Iai)(r+iz,s+iz)+=complex<double>(Jai(ai+1,bi+1),0.0);}   
//        else if(r>s)   {(*Iai)(r+iz,s+iz)+=complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1));}
         if(r==s)     { if(fabs(Jai(ai+1,bi+1))>EPS)                           (*Iai)(r+iz,s+iz) = complex<double>(Jai(ai+1,bi+1),0.); }
         else if(r>s) { if(fabs(Jai(ai+1,bi+1))>EPS||fabs(Jai(bi+1,ai+1))>EPS) (*Iai)(r+iz,s+iz) = complex<double>(Jai(ai+1,bi+1),Jai(bi+1,ai+1)); }
         else if(r<s) { if(fabs(Jai(ai+1,bi+1))>EPS||fabs(Jai(bi+1,ai+1))>EPS) (*Iai)(r+iz,s+iz) = complex<double>(Jai(bi+1,ai+1),-Jai(ai+1,bi+1)); }
      } } }
   (*Iai).to_csc();
//}
 }
}

//------------------------------------------------------------------------------------------------
//routine Icalc for cluster
//------------------------------------------------------------------------------------------------
void jjjpar::cluster_Icalc_mcalc_Micalc (int code,Vector & Jret,double & T, Vector &  Hxc,Vector & Hext, double & lnZ, double & U)
{ /*on input
   code         defining, what should be calculated
         1....   Ia (interaction operators for Icalc)
         2....    m (total magnetic moment for mcalc), or
         3....    Mi (individual magnetic moments (for MQ)
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    J		single ion momentum vector <J>
    Z		single ion partition function
    U		single ion magnetic energy
*/


Vector En(1,dim);
//Matrix zr(1,dim,1,dim);
ComplexMatrix zc(1,dim,1,dim);
cluster_calcH_and_diagonalize(En,zc,Hxc,Hext);

// calculate Z and wn (occupation probability)
     Vector wn(1,dim);double Zs;
     double x,y;
     x=Min(En);
    if(T>0)
     {
     for (int i=1;i<=dim;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y);
      else wn[i]=0.0;
//      printf("%4.4g\n",En(i));
      }
     Zs=Sum(wn);wn/=Zs;
     lnZ=log(Zs)-x/KB/T;
     }
     else
     { printf ("Temperature T<0: please choose probability distribution of states by hand\n");
                         printf ("Number   Energy     Excitation Energy\n");
     for (int i=1;i<=dim;++i) printf ("%i    %4.4g meV   %4.4g meV\n",i,En(i),En(i)-x);
     char instr[MAXNOFCHARINLINE];
     for (int i=1;i<=dim;++i)
      {printf("eigenstate %i: %4.4g meV %4.4g meV  - please enter probability w(%i):",i,En(i),En(i)-x,i);
       fgets(instr, MAXNOFCHARINLINE, stdin);

       wn(i)=strtod(instr,NULL);
      }
       Zs=Sum(wn);wn/=Zs;

       lnZ=log(Zs);
                         printf ("\n\nNumber   Energy     Excitation Energy   Probability\n");
     for (int i=1;i<=dim;++i) printf ("%i    %4.4g meV   %4.4g meV %4.4g  \n",i,En(i),En(i)-x,wn(i));
     }

   // calculate U
     U=En*wn;

 switch(code)

{case 1:for(int a=1;a<=Hxc.Hi();++a)
        {// calculate expecation Value of J
         //printf("Matrix of Operator");
         // myPrintMatrix(stdout,Ja);
         // determine expectation value
         Jret(a)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {Jret(a)+=wn[i]*aMb_real((*Ia[a]),zc,i,i); 
          // if(fabs(aMb_imag((*Ia[a]),zr,zc,i,i))>SMALL){fprintf(stderr,"ERROR module cluster - Icalc: expectation value imaginary\n");exit(EXIT_FAILURE);}
         }
        }break;
 case 2:for(int n=1;n<=3;++n)
        {// calculate expecation Value of m
         Jret(n)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {Jret(n)+=wn[i]*aMb_real((*cluster_M[n]),zc,i,i); 
         }
        }break;
 case 3:for(int a=1;a<=(*clusterpars).cs.nofatoms;++a)for(int n=1;n<=3;++n)
        {int index_M=a*3+n; // a .... atom index  n ... xyz components of magnetic moment
         int index=(a-1)*3+n; // calculate expecation Value of m
         Jret(index)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {Jret(index)+=wn[i]*aMb_real((*cluster_M[index_M]),zc,i,i); 
         }
        }break;
}

//  printf ("Ha=%g Hb=%g Hc=%g Ja=%g Jb=%g Jc=%g \n", 
//     gjmbH[1]/MU_B/gjJ, gjmbH[2]/MU_B/gjJ, gjmbH[3]/MU_B/gjJ, J[1], J[2], J[3]);
}
//------------------------------------------------------------------------------------------------
//routine ducalc for cluster
//------------------------------------------------------------------------------------------------

int jjjpar::cluster_dm(int code,int & tn,double & T,ComplexVector & u1,float & delta,int & ni, int & nf,ComplexMatrix & ests)
{ 
  /*on input
   delta        maxE
   code         defining, what should be calculated
         1....   Ia (interaction operators for du1calc)
         2....    m (total magnetic moment for dm1calc), or
         3....    Mi (individual magnetic moments (for dMQ1)
    transitionnumber ... number of transition to be computed - meaningless for kramers doublet, because there is only 1 transition
    T		temperature[K]
    gjmbH	vector of effective field [meV]
  on output    
    delta	splitting of kramers doublet [meV]
    u1(i)	<-|(Ji-<Ji>)|+> sqrt(tanh(delta/2kT))
    ni,nf        numer of level |-> and |+>
*/
  int pr=0,subtractexpvalue=1;if(T<0){subtractexpvalue=0;T=-T;}
  if (tn<0) {pr=1;tn*=-1;}
   double ninit=u1[1].real();
   double pinit=u1[1].imag();
   double maxE=delta; 
// 1. get i and j from tn
int k=0,ii=1,jj=1;
for(ii=1;ii<=dim;++ii){for(jj=ii;jj<=dim;++jj){++k;if(k==tn)break;}if(k==tn)break;}
ni=ii;nf=jj;
if((delta=real(ests(0,jj))-real(ests(0,ii)))<=maxE)
 {
    Vector En(1,dim);
    Vector wn(1,dim);
    for(int i=1;i<=dim;++i){En(i)=real(ests(0,i));wn(i)=imag(ests(0,i));}

  //Matrix zr(1,dim,1,dim);
  //Matrix zc(1,dim,1,dim);
  //for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){zr(i,j)=real(ests(i,j));zc(i,j)=imag(ests(i,j));}

    // calculate mat and delta for transition number tn
    // 2. set delta
    delta=En(jj)-En(ii);
    if (delta<-0.000001){fprintf(stderr,"ERROR module cluster - du1calc: energy gain delta gets negative\n");exit(EXIT_FAILURE);}
    if(jj==ii)delta=-SMALL_QUASIELASTIC_ENERGY; //if transition within the same level: take negative delta !!- this is needed in routine intcalc

   int maxn=0;
   switch(code)
   {case 1: maxn=u1.Hi();break;
    case 2: maxn=3; break;
    case 3: maxn=(*clusterpars).cs.nofatoms*3;break;
   }

   Vector Jret(1,maxn);Jret=0;
   ComplexVector iJj(1,maxn); //Matrix aMb_tempV(1,6,1,dim); //double rv, iv;
   for(int a=1;a<=maxn;++a)
   {// transition matrix element
    switch(code)
    {case 1: iJj(a) = aMb_complex((*Ia[a]),ests,ii,jj); break;
     case 2: iJj(a) = aMb_complex((*cluster_M[a]),ests,ii,jj); break;
     case 3: iJj(a) = aMb_complex((*cluster_M[a+3]),ests,ii,jj); break;
    }
    // determine expectation value
    Jret(a)=0;
    if (subtractexpvalue==1&&ii==jj)
    { for(int i=1;i<=dim&&wn[i]>0.00001;++i)
     { switch(code)
      {case 1: Jret(a)+=wn[i]*aMb_real((*Ia[a]),ests,i,i);break;
       case 2: Jret(a)+=wn[i]*aMb_real((*cluster_M[a]),ests,i,i); break;
       case 3: Jret(a)+=wn[i]*aMb_real((*cluster_M[a+3]),ests,i,i);break;
      }
     }
    }
   }

   // 3. set u1
   for(int l=1;l<=maxn;++l)
   {if(ii==jj){//take into account thermal expectation values <Jret>
             u1(l)=(iJj(l)-Jret(l));}
    else    {u1(l)=iJj(l);} 
   }

   //printf("code %i: %g %g %g\n",code,Jret(1),Jret(2),Jret(3));
   //myPrintMatrix(stdout,(*cluster_M[1]));
   //myPrintMatrix(stdout,(*cluster_M[2]));
   //myPrintMatrix(stdout,(*cluster_M[3]));

   if (delta>SMALL_QUASIELASTIC_ENERGY)
      { if(pr==1){
         printf("delta(%i->%i)=%4.4gmeV",ii,jj,delta);
         printf(" |<%i|Ja|%i>|^2=%4.4g |<%i|Jb|%i>|^2=%4.4g |<%i|Jc|%i>|^2=%4.4g",ii,jj,abs(u1(1))*abs(u1(1)),ii,jj,abs(u1(2))*abs(u1(2)),ii,jj,abs(u1(3))*abs(u1(3)));
         printf(" n%i-n%i=%4.4g\n",ii,jj,wn(ii)-wn(jj));}
       u1*=sqrt(wn(ii)-wn(jj)); // occupation factor
        }else
      {// quasielastic scattering has not wi-wj but wj*epsilon/kT
        if(pr==1){
         printf("delta(%i->%i)=%4.4gmeV",ii,jj,delta);
         printf(" |<%i|Ja-<Ja>|%i>|^2=%4.4g |<%i|Jb-<Jb>|%i>|^2=%4.4g |<%i|Jc-<Jc>|%i>|^2=%4.4g",ii,jj,abs(u1(1))*abs(u1(1)),ii,jj,abs(u1(2))*abs(u1(2)),ii,jj,abs(u1(3))*abs(u1(3)));
         printf(" n%i=%4.4g\n",ii,wn(ii));}
       u1*=sqrt(wn(ii)/KB/T);
      }


   if (pr==1) printf ("delta=%4.6g meV\n",delta);
}

   // return number of all transitions
      double n=ninit;
      if (n>dim)n=dim;
      double zsum=0,zii,x;
      int noft=0;
      if(T>0)for(int i=1;(i<=n)&((((x=(real(est(0,i))-real(est(0,1)))/KB/T)<200)? zii=exp(-x):zii=0)>=(pinit*zsum));++i)
      {noft+=dim-i+1;zsum+=zii;}
   
      if(T<0)for(int i=1;i<=n;++i)
      {noft+=dim-i+1;}
      //printf("nt=%i",noft);
 return noft;



}

int arpackeig(zsMat<double> &M, Vector &En, complexdouble*zc, int nev, iterwork &workspace)
{
   int n = M.nr();
   int ido=0,ncv=(2*nev>n?n:2*nev),iparam[]={1,0,10000,1,nev,0,1,0,0,0,0};
   //         IPARAM=ISHIFT,,MAXITER,NB,NCONV,,MODE,NP,NUMOP,NUMOPB,NUMREO
   int lworkl=3*ncv*ncv+8*ncv,info=0;  // Must be > 3*NCV**2 + 5*NCV
   char bmat='I';
   double tol=1e-12;

   // Check workspace is large enough. 
   int litwk=n*(n+ncv)+n*4+lworkl+nev+2*ncv;
   if(workspace.zsize<litwk) workspace.realloc_z(litwk); memset(workspace.zwork,0,litwk*sizeof(complexdouble));
   if(workspace.dsize<ncv)   workspace.realloc_d(ncv);   memset(workspace.dwork,0,ncv*sizeof(double));
   if(workspace.isize<ncv)   workspace.realloc_i(ncv);   memset(workspace.iwork,0,ncv*sizeof(int));

   complexdouble zalpha; zalpha.r=1; zalpha.i=0; //complexdouble zbeta; zbeta.r=0; zbeta.i=0;

   char whichp[]="SR";       // LM/SM==Largest/Smallest Magnitude; LR/SR==RealPart; LI/SI==ImaginaryPart
   int ipntr[14];
   complexdouble *z=&workspace.zwork[0], *resid=&workspace.zwork[n*n], *v=&workspace.zwork[n*n+n], *workd=&workspace.zwork[n*(n+ncv)+n];
   complexdouble *workl=&workspace.zwork[n*(n+ncv)+n*4], *d=&workspace.zwork[n*(n+ncv)+n*4+lworkl], *workev=&workspace.zwork[n*(n+ncv)+n*4+lworkl+nev];
   double *rwork = &workspace.dwork[0]; int *select = &workspace.iwork[0];
   while(ido!=99)
   {
      F77NAME(znaupd)(&ido, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam, ipntr, workd, workl, &lworkl, rwork, &info);
      if(ido==1 || ido==-1)
         M.MultMv((std::complex<double>*)&workd[ipntr[0]],(std::complex<double>*)&workd[ipntr[1]]);
   }
   int rvec=1; char howmny='A';
   F77NAME(zneupd)(&rvec, &howmny, select, d, z, &n, &zalpha, workev, &bmat, &n, whichp, &nev, &tol, resid, &ncv, v, &n, iparam,
                   ipntr, workd, workl, &lworkl, rwork, &info);
   int i=1,j=2,ii; double elem; std::vector<int> ind(nev,0);
   for(ii=0; ii<nev; ii++) { En[ii+1] = d[ii].r; ind[ii] = ii; }
   while(i<nev)
   {
      if(En[i]<=En[i+1]) { i=j; j++; }
      else { elem = En[i]; En[i] = En[i+1]; En[i+1] = elem; ii=ind[i-1]; ind[i-1]=ind[i]; ind[i]=ii; i--; if(i==0) i=1; }
   }
   for(i=0; i<nev; i++) {
      for(j=0; j<n; j++) zc[n*j+i]=z[ind[i]*n+j]; //{ zr(j+1,i+1)=z[ind[i]*n+j].r; zc(j+1,i+1)=z[ind[i]*n+j].i; }
   }

   return info;
}


//------------------------------------------------------------------------------------------------
//routine estates for cluster
//------------------------------------------------------------------------------------------------
void jjjpar::cluster_est(ComplexMatrix * eigenstates,Vector &Hxc,Vector &Hext,double & T)
{ /*on input
    gJmbH	vector of effective field [meV]
      on output
    Matrix containing the eigenvalues and eigenfunctions of the crystalfield problem
    eigenvalues ares stored as real part of row zero
    boltzmann population numbers are stored as imaginary part of row zero
*/
 fprintf(stderr,"# calculating eigenstates of cluster ...");
(*eigenstates) = ComplexMatrix(0,dim,1,dim);
 Vector En(1,dim);
 //Matrix zr(1,dim,1,dim);
 //ComplexMatrix zc(1,dim,1,dim);
 cluster_calcH_and_diagonalize(En,(*eigenstates),Hxc,Hext);
     Vector wn(1,dim);double Zs;
     double x,y;
     x=Min(En);
     for (int i=1;i<=dim;++i)
     {if ((y=(En(i)-x)/KB/T)<700) wn[i]=exp(-y);
      else wn[i]=0.0;
      }
     Zs=Sum(wn);wn/=Zs;
 for(int i=1;i<=dim;++i)(*eigenstates)(0,i)=complex<double>(En(i),wn(i));

 //for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j)
 //{(*eigenstates)(i,j)=complex<double>(zr(i,j),zc(i,j));
 //}
 fprintf(stderr," ... done\n");
}

void jjjpar::cluster_calcH_and_diagonalize(Vector & En,ComplexMatrix &zc,Vector & Hxc,Vector & Hext)
{zsMat<double> H(dim,dim); 
 H.zero(); // initialize to zero

 bool isHextsame=true; 
 if(justinit) { isHextsame = false; justinit=false; }
 else
  for(int i=Hext.Lo(); i<=Hext.Hi(); i++) { if(fabs(Hext[i]-(*oldHext)[i])>DBL_EPSILON) { isHextsame=false; break; } }
 if(!isHextsame) {
    (*clusterH).zero(); (*clusterH).to_tri();
    Vector ZeroHxc(1,1);ZeroHxc=0;
// fill H matrix with sum over Hi of individual spins
for (int i=1;i<=(*clusterpars).cs.nofatoms;++i)
{Matrix Hi((*(*clusterpars).jjj[i]).opmat(0,ZeroHxc,Hext)); // here we need ZeroHxc because
                                                            // exchange operators are set by user (Ia)
                                                            // and not in the single ion module
  //1. determine dimensions
  int di=dnn[i];
  int dx=1,dz=1;
  for(int m=1  ;m<i;++m)dx*=dnn[m];
  for(int m=i+1;m<=(*clusterpars).cs.nofatoms;
                    ++m)dz*=dnn[m];
  int mx=0,mi=1;
  if(i>1)mx=1;
  for(int m=i;  m<=(*clusterpars).cs.nofatoms;++m)mx*=dnn[m];
  for(int m=i+1;m<=(*clusterpars).cs.nofatoms;++m)mi*=dnn[m];

  //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int bi=0;bi<=di-1;++bi)
  {int r=ix*mx+ai*mi+iz+1;
   int s=ix*mx+bi*mi+iz+1;
//   printf("hello %i %i %i %i\n",r,s,ai,bi);
   //if(fabs(Hi(ai+1,bi+1))>EPS) { if(r>=s) (*clusterH)(r,s)+=Hi(ai+1,bi+1); else (*clusterH)(s,r)+=complex<double>(0.,Hi(ai+1,bi+1)); }
   if(r>=s) (*clusterH)(r,s)+=Hi(ai+1,bi+1); else (*clusterH)(s,r)+=complex<double>(0.,Hi(ai+1,bi+1));
   //(*clusterH)(r,s)+=Hi(ai+1,bi+1);
  }
}
// myPrintMatrix(stdout,H);
//printf("now enter interactions\n");
// put interactions in H matrix
for (int n=1;n<=(*clusterpars).cs.nofatoms;++n)
for (int nn=1;nn<=(*(*clusterpars).jjj[n]).paranz;++nn)
{// determine interaction operator of atom n with neighbour nn
 // get i and j from n and sublattice[nn]
 int i=n;int j=(*(*clusterpars).jjj[n]).sublattice[nn];
 int ilj=1;
 if(i>j){i=j;j=n;ilj=0;} // if n>nn exchange i and j ... so that we always have i < j
 Matrix SinS(-0.5*(*(*clusterpars).jjj[n]).jij[nn](1,1)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(1,ZeroHxc,Hext),
                          (*(*clusterpars).jjj[j]).opmat(1,ZeroHxc,Hext)
                         )
             );

 for(int a=1;a<=(*(*clusterpars).jjj[i]).nofcomponents;++a)
 for(int b=1;b<=(*(*clusterpars).jjj[j]).nofcomponents;++b)
 {
  if(a==1&&b==1) break;
  if(ilj){ if(fabs((*(*clusterpars).jjj[n]).jij[nn](a,b))>DBL_EPSILON) 
  SinS+=-0.5*(*(*clusterpars).jjj[n]).jij[nn](a,b)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(a,ZeroHxc,Hext),
                          (*(*clusterpars).jjj[j]).opmat(b,ZeroHxc,Hext)
                    );
          } else { // if n>nn then take transpose of exchange parameter jij !! (J(ij)=JT(ji))
           if(fabs((*(*clusterpars).jjj[n]).jij[nn](b,a))>DBL_EPSILON) 
  SinS+=-0.5*(*(*clusterpars).jjj[n]).jij[nn](b,a)*
             herm_dirprod((*(*clusterpars).jjj[i]).opmat(a,ZeroHxc,Hext),
                          (*(*clusterpars).jjj[j]).opmat(b,ZeroHxc,Hext)
                    );
          }
 }
//printf("hello n=%i nn=%i i=%i j=%i ...\n",n,nn,i,j);  myPrintMatrix(stdout,SinS);
 // insert SinS into H
 //1. determine dimensions
  int di=dnn[i];
  int dj=dnn[j];
 int dx=1,dy=1,dz=1;
 for(int m=1  ;m<i;++m)dx*=dnn[m];
 for(int m=i+1;m<j;++m)dy*=dnn[m];
 for(int m=j+1;m<=(*clusterpars).cs.nofatoms;
                   ++m)dz*=dnn[m];

 int mx=0,my=0,mi=1,mj=1;
 if(i>1)mx=1;
 if(i+1<=(*clusterpars).cs.nofatoms&&i+1<j)my=1;
 
 for(int m=i;  m<=(*clusterpars).cs.nofatoms;++m)mx*=dnn[m];
 for(int m=i+1;m<=(*clusterpars).cs.nofatoms;++m)mi*=dnn[m];
 for(int m=j  ;m<=(*clusterpars).cs.nofatoms;++m)my*=dnn[m];
 for(int m=j+1;m<=(*clusterpars).cs.nofatoms;++m)mj*=dnn[m];

 //2. insertion loop
  for(int ix=0;ix<=dx-1;++ix)
  for(int iy=0;iy<=dy-1;++iy)
  for(int iz=0;iz<=dz-1;++iz)
  for(int ai=0;ai<=di-1;++ai)
  for(int aj=0;aj<=dj-1;++aj)
  for(int bi=0;bi<=di-1;++bi)
  for(int bj=0;bj<=dj-1;++bj)
  {int r=ix*mx+ai*mi+iy*my+aj*mj+iz+1;
   int s=ix*mx+bi*mi+iy*my+bj*mj+iz+1;
   int k=ai*dj+aj+1;// ??
   int l=bi*dj+bj+1;// ??
   //if(fabs(SinS(k,l))>EPS) { if(r>=s) (*clusterH)(r,s)+=SinS(k,l); else (*clusterH)(s,r)+=complex<double>(0.,SinS(k,l)); }
   if(r>=s) (*clusterH)(r,s)+=SinS(k,l); else (*clusterH)(s,r)+=complex<double>(0.,SinS(k,l));
   //(*clusterH)(r,s)+=SinS(k,l);
//printf("hello %i %i %i %i %g %g\n",r,s,k,l,H(r,s),SinS(k,l));
  } 
}

if (truncate>1e-6 && truncate<1) 
{
   char jobz = 'V', uplo = 'U'; int lda=fdim, n=fdim, info=0, lwork=4*n, lrwork = 3*n-2;
   int zsz = lwork + fdim*dim + dim*dim, dsz = lrwork+fdim;

   if(workspace->zsize<zsz) workspace->realloc_z(zsz); memset(workspace->zwork,0,zsz*sizeof(complexdouble));
   if(workspace->dsize<dsz) workspace->realloc_d(dsz); memset(workspace->dwork,0,dsz*sizeof(double));
   complexdouble *zwork=&workspace->zwork[0], *zmt=&workspace->zwork[lwork]; 
   complexdouble *opr=&workspace->zwork[lwork+fdim*dim], zv;
   double *rwork=&workspace->dwork[0], *eigv=&workspace->dwork[lrwork];

   // Diagonalise the full Hamiltonian
   if(is1sttrunc)
   {
      clock_t start,end; start = clock();
      printf("Truncate cluster: full dimension=%d\ttruncated dimension=%d\n",fdim,dim); fflush(stdout);
      (*clusterH).h_array((std::complex<double>*)zm);
      if(fdim>200) { printf("Diagonalising Full Hamiltonian..."); fflush(stdout); }
      F77NAME(zheev)(&jobz, &uplo, &n, zm, &lda, eigv, zwork, &lwork, rwork, &info);
      if(info!=0) { 
         fprintf(stderr,"cluster_module:zheev return error %d\n",info); exit(-1); }
    //memset(zm,0,fdim*dim*sizeof(complexdouble)); arpackeig(*clusterH,En,zm,dim,*workspace); if(workspace->dsize!=dsz) { rwork=&workspace->dwork[0]; }
    //if(workspace->zsize!=zsz) { zwork=&workspace->zwork[0]; zmt=&workspace->zwork[lwork]; opr=&workspace->zwork[lwork+fdim*dim]; }
      end = clock(); if(fdim>200) printf("done. In %f s\n",(double)(end-start)/CLOCKS_PER_SEC); fflush(stdout);
      if(fdim>200) { printf("eigval="); for(int ii=1; ii<10; ii++) printf("%g ",eigv[ii]); printf(" ... "); for(int ii=fdim-10; ii<fdim; ii++) printf("%g ",eigv[ii]); printf("\n"); }
    //if(fdim>200) { printf("eigval="); for(int ii=1; ii<10; ii++) printf("%g ",En[ii]); printf(" ... "); for(int ii=dim-10; ii<dim; ii++) printf("%g ",En[ii]); printf("\n"); }
   }
   // Use the eigenvectors to transform the full Hamiltonian
   char notranspose='N',transpose='C'; complexdouble zalpha; zalpha.r=1; zalpha.i=0; complexdouble zbeta; zbeta.r=0; zbeta.i=0;
   (*clusterH).MultMMH((std::complex<double>*)zmt,(std::complex<double>*)zm,dim);
   F77NAME(zgemm)(&transpose,&notranspose,&dim,&dim,&fdim,&zalpha,zm,&fdim,zmt,&fdim,&zbeta,opr,&dim);
   (*clusterH).zero(dim,dim); (*clusterH).to_tri();
   for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) { 
      zv=opr[dim*j+i]; if(fabs(zv.r)>DBL_EPSILON || fabs(zv.i)>DBL_EPSILON) (*clusterH)(i+1,j+1) = std::complex<double>(zv.r,zv.i); }
   // Transforms all the operator matrices
   if (is1sttrunc)
   {
      is1sttrunc=false;
      clock_t start,end; start = clock();
      if(fdim>200) { printf("Rotating operator matrices..."); fflush(stdout); }
      for(int a=1; a<=nofcomponents; a++)
      {
         (*Ia[a]).MultMMH((std::complex<double>*)zmt,(std::complex<double>*)zm,dim);
         F77NAME(zgemm)(&transpose,&notranspose,&dim,&dim,&fdim,&zalpha,zm,&fdim,zmt,&fdim,&zbeta,opr,&dim);
         (*Ia[a]).zero(dim,dim); (*Ia[a]).to_tri();
         for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) { 
            zv=opr[dim*j+i]; if(fabs(zv.r)>DBL_EPSILON || fabs(zv.i)>DBL_EPSILON) (*Ia[a])(i+1,j+1) = std::complex<double>(zv.r,zv.i); }
      }
      for(int a=1; a<=(*clusterpars).cs.nofatoms*3+3; a++)
      {
       //(*cluster_M[a]).to_csc();
         (*cluster_M[a]).MultMMH((std::complex<double>*)zmt,(std::complex<double>*)zm,dim);
         F77NAME(zgemm)(&transpose,&notranspose,&dim,&dim,&fdim,&zalpha,zm,&fdim,zmt,&fdim,&zbeta,opr,&dim);
         (*cluster_M[a]).zero(dim,dim); (*cluster_M[a]).to_tri();
         for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) { 
            zv=opr[dim*j+i]; if(fabs(zv.r)>DBL_EPSILON || fabs(zv.i)>DBL_EPSILON) (*cluster_M[a])(i+1,j+1) = std::complex<double>(zv.r,zv.i); }
      }
      end = clock(); if(fdim>200) { printf("done. In %f s\n",(double)(end-start)/CLOCKS_PER_SEC); fflush(stdout); }
   }
}
else
{
   (*clusterH).to_csc();
   (*clusterH).tridupl(); for(int i=1; i<=Hxc.Hi(); i++) (*clusterH) += 0.*(*Ia[i]); (*clusterH).to_csc(); 
}

*oldHext = Hext; } H = *clusterH;

if(oldeig) {
// insert exchange field
for(int i =1;i<=Hxc.Hi();++i)H-=Hxc(i)*(*Ia[i]);
// diagonalize H
int sort=1;int maxiter=1000000;
Matrix Hp = H.fp_matrix();
// myPrintMatrix(stdout,H);
//printf("now diagonalise\n");
Matrix zr(1,dim,1,dim);
Matrix zi(1,dim,1,dim);
EigenSystemHermitean (Hp,En,zr,zi,sort,maxiter);
for(int i=1; i<=dim; i++) { for(int j=1; j<=dim; j++) { zc(i,j) = std::complex<double>(zr(i,j),zi(i,j)); } }
//printf("Eigenvector real\n");
// myPrintMatrix(stdout,zr);
//printf("Eigenvector imag\n"); 
// myPrintMatrix(stdout,zc);
}
else if (arpack>1e-6 && arpack!=1 && !(truncate>1e-6 && truncate<1) ) {
   int nev = ceil(arpack*dim);
   memset(&zc[1][1],0,dim*dim*sizeof(complexdouble)); 
   arpackeig(*clusterH,En,(complexdouble*)&zc[1][1],nev,*workspace);
   double emax = fabs(En[nev-1]*100); if(fabs(En[1]*100)>emax) emax = fabs(En[1]*100);
   for (int ii=nev; ii<dim; ii++) { En[ii+1] = emax; }
}
else
{
   clock_t start,end; start = clock();
   char jobz = 'V', uplo = 'U', range = 'A'; int n=dim, lda=n, ldz=n, info=0, lwork=4*n, il, iu, numfnd, lrwork = 24*n, liwork=10*n;
   double vl=-9e9, vu=9e9, abstol = 0.00001;
   int zsz = lwork + dim*dim*2, dsz = lrwork, isz = 2*n + liwork;
   if(workspace->isize<isz) workspace->realloc_i(isz); memset(workspace->iwork,0,isz*sizeof(int));
   if(workspace->zsize<zsz) workspace->realloc_z(zsz); memset(workspace->zwork,0,zsz*sizeof(complexdouble));
   if(workspace->dsize<dsz) workspace->realloc_d(dsz); memset(workspace->dwork,0,dsz*sizeof(double));
   complexdouble *zwork = &workspace->zwork[0], *zmt = &workspace->zwork[lwork], *zmo = &workspace->zwork[lwork+dim*dim];  
   double *rwork = &workspace->dwork[0];
   int *isuppz = &workspace->iwork[0], *iwork = &workspace->iwork[2*n];  
   for(int i =1;i<=Hxc.Hi();++i)H-=Hxc(i)*(*Ia[i]); H.h_array((std::complex<double>*)zmt);
   F77NAME(zheevr)(&jobz, &range, &uplo, &n, zmt, &lda, &vl, &vu, &il, &iu, &abstol, &numfnd, &En[1],
          zmo, &ldz, isuppz, zwork, &lwork, rwork, &lrwork, iwork, &liwork, &info);
   for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) { 
      zc(i+1,j+1)=std::complex<double>(zmo[dim*j+i].r,zmo[dim*j+i].i); }
   if(info!=0) { 
      fprintf(stderr,"cluster_module:zheevr return error %d\n",info); exit(-1); }
   end = clock(); if(dim>200) fprintf(stderr,"zheev took %f s\n",(double)(end-start)/CLOCKS_PER_SEC); fflush(stderr);
}
/*{
   clock_t start,end; start = clock();
   char jobz = 'V', uplo = 'U'; int n=dim, lda=n, info=0, lwork=2*n;
   int zsz = lwork + dim*dim*2, dsz = 3*n-2;
   if(workspace->zsize<zsz) workspace->realloc_z(zsz); memset(workspace->zwork,0,zsz*sizeof(complexdouble));
   if(workspace->dsize<dsz) workspace->realloc_d(dsz); memset(workspace->dwork,0,dsz*sizeof(double));
   complexdouble *zwork = &workspace->zwork[0], *zmt = &workspace->zwork[lwork];
   double *rwork = &workspace->dwork[0];
   for(int i =1;i<=Hxc.Hi();++i)H-=Hxc(i)*(*Ia[i]); H.h_array((std::complex<double>*)zmt);
   F77NAME(zheev)(&jobz, &uplo, &n, zmt, &lda, &En[1], zwork, &lwork, rwork, &info);
   for (int i=0; i<dim; i++) for (int j=0; j<dim; j++) { 
      zc(i+1,j+1)=std::complex<double>(zmt[dim*j+i].r,zmt[dim*j+i].i); }
   if(info!=0) { 
      fprintf(stderr,"cluster_module:zheev return error %d\n",info); exit(-1); }
   end = clock(); if(dim>200) fprintf(stderr,"zheev took %f s\n",(double)(end-start)/CLOCKS_PER_SEC); fflush(stderr);
}*/

}

/**************************************************************************/
//                   OBSERVABLES
/**************************************************************************/

void jjjpar::cluster_Micalc (Vector &mom,ComplexMatrix & ests)
{
Vector En(1,dim);
Vector wn(1,dim);
//Matrix zr(1,dim,1,dim);
//Matrix zc(1,dim,1,dim);
 for(int i=1;i<=dim;++i){En(i)=real(ests(0,i));wn(i)=imag(ests(0,i));}

   //for(int i=1;i<=dim;++i)for(int j=1;j<=dim;++j){zr(i,j)=real(ests(i,j));zc(i,j)=imag(ests(i,j));}

for(int a=1;a<=(*clusterpars).cs.nofatoms;++a)for(int n=1;n<=3;++n)
        {int index_M=a*3+n; // a .... atom index  n ... xyz components of magnetic moment
         int index=(a-1)*3+n; // calculate expecation Value of m
         mom(index)=0;
         for(int i=1;i<=dim&&wn[i]>0.00001;++i)
         {mom(index)+=wn[i]*aMb_real((*cluster_M[index_M]),ests,i,i); 
         }
        }
}
