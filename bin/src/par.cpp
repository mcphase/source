#include "par.hpp"
#include "../../version"
#include <martin.h>
#include <cstring>

#define MAXNOFNUMBERSINLINE 20
#define MAXNOFCHARINLINE 7024
#define SMALL_MATCH_LATTICEVECTOR 1e-4


 // *************************************************************************
 // ************************ parameters *************************************
 // *************************************************************************


//constructor 
par::par(float ai,float bi,float ci,float alphai,float betai,float gammai,int nofci)
{r=Matrix(1,3,1,3);rez=Matrix(1,3,1,3);Cel=Matrix(1,6,1,6);CelInv=Matrix(1,6,1,6);CelInv=0;
 a=ai;b=bi;c=ci;alpha=alphai;beta=betai;gamma=gammai;nofcomponents=nofci;
  rems[1]=new char[40];rems[1][0]='\0';
  rems[2]=new char[40];rems[2][0]='\0';
  rems[3]=new char[40];rems[3][0]='\0';
  nofatoms=0;
  jjj=new jjjpar * [nofatoms+1];
}

par::par (const char *filejjj)
{ int i,j,n,l;
  FILE *fin_coq;
  char instr[MAXNOFCHARINLINE];
  Vector hkl(1,3),hkl_rint(1,3);
  r=Matrix(1,3,1,3);rez=Matrix(1,3,1,3);
    Cel=Matrix(1,6,1,6);Cel=0;CelInv=Matrix(1,6,1,6);CelInv=0;
  fin_coq = fopen_errchk (filejjj, "rb");

 // input file header ------------------------------------------------------------------
  fgets (instr, MAXNOFCHARINLINE, fin_coq);
  extract(instr,"a",a);extract(instr,"b",b); extract(instr,"c",c); 
                 extract(instr,"alpha",alpha);  extract(instr,"beta",beta);extract(instr,"gamma",gamma); 
  instr[0]='#';
   // inserted 12.11.07 in order to format output correctly (characterstring 13 spoiled output string)
   for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
   rems[1]=new char[strlen(instr)+2];strcpy(rems[1],instr);
   rems[2]=new char[40];strcpy(rems[2],"#!<--mcphas.mcphas.j-->");
  nofatoms=0;
  instr[0]='#';a=0;b=0;c=0;
 while (nofatoms==0||(strstr(instr,"*******")==NULL&&instr[strspn(instr," \t")]=='#')) 
  {fgets(instr,MAXNOFCHARINLINE,fin_coq);
   if(a==0&&b==0&&c==0){extract(instr,"a",a);extract(instr,"b",b); extract(instr,"c",c); 
                 extract(instr,"alpha",alpha);  extract(instr,"beta",beta);extract(instr,"gamma",gamma); 
   }
   extract(instr,"r1x",r[1][1]);extract(instr,"r2x",r[1][2]); extract(instr,"r3x",r[1][3]); 
   extract(instr,"r1y",r[2][1]); extract(instr,"r2y",r[2][2]); extract(instr,"r3y",r[2][3]);
   extract(instr,"r1z",r[3][1]); extract(instr,"r2z",r[3][2]); extract(instr,"r3z",r[3][3]);
   extract(instr,"r1a",r[1][1]);extract(instr,"r2a",r[1][2]); extract(instr,"r3a",r[1][3]); 
   extract(instr,"r1b",r[2][1]); extract(instr,"r2b",r[2][2]); extract(instr,"r3b",r[2][3]);
   extract(instr,"r1c",r[3][1]); extract(instr,"r2c",r[3][2]); extract(instr,"r3c",r[3][3]);

    // read optional elastic constants        
  char Celstr[6];
   for(i=1;i<=6;++i)for(j=1;j<=6;++j){
  sprintf(Celstr,"Cel%i%i",i,j);// printf("%s\n",Celstr);
  extract(instr,Celstr,Cel(i,j));Cel(j,i)=Cel(i,j);}

   extract(instr,"nofatoms",nofatoms);extract(instr,"nofcomponents",nofcomponents); 
		  if(feof(fin_coq)!=0)
                    {fprintf(stderr,"ERROR reading header of file %s: line '#! nofatoms=...' not found\n",filejjj);exit(EXIT_FAILURE);}
  }
  if(nofatoms>MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL)
  {fprintf(stderr,"ERROR reading mcphas.j: maximum number of atoms in unit cell exceeded - enlarge it in par.hpp and recompile\n");exit(EXIT_FAILURE);}
  

  rez=r.Inverse();
  rems[3]=new char[strlen(instr)+2];strcpy(rems[3],instr);
  
  //read parameter sets for every atom 
  jjj=new jjjpar * [nofatoms+1];
  //gJ=Vector(1,nofatoms);
  for(i=1;i<=nofatoms;++i)  
  {//printf("creating atom %i (of %i)...\n",i,nofatoms);
   jjj[i]=new jjjpar(fin_coq,nofcomponents);
   if(jjj[i]==NULL){ fprintf (stderr, "Out of memory creating atoms jjjpar by constructor\n");exit (EXIT_FAILURE);}
   //gJ(i)=(*jjj[i]).gJ;
   
   if(nofcomponents!=(*jjj[i]).nofcomponents)
   {fprintf(stderr,"ERROR reading mcphas.j: nofcomponents (%i) not consistent for atom %i (%i read in fileheader)\n",(*jjj[i]).nofcomponents,i,nofcomponents);exit(EXIT_FAILURE);}
  }
  //determine sublattices
  for(i=1;i<=nofatoms;++i)
  {for(n=1;n<=(*jjj[i]).paranz;++n)
   {(*jjj[i]).sublattice[n]=0;
    for(j=1;j<=nofatoms;++j)
    {//try if neighbour n is on sublattice j
     hkl=rez*((*jjj[i]).xyz+(*jjj[i]).dn[n]-(*jjj[j]).xyz);
     // check if hkl is integer - if yes then the neighbour n is on sublattice j
     for(l=1;l<=3;++l){hkl_rint(l)=rint(hkl(l));}
     if(Norm(hkl_rint-hkl)<0.001){(*jjj[i]).sublattice[n]=j;}
    }
    if((*jjj[i]).sublattice[n]==0){fprintf(stderr,"Warning mcphas - par.cpp: file %s inconsistent:  neighbour %i of atom %i at %g %g %g is not on any sublattice. Continuing putting it onto sublattice 1 ...\n",filejjj,n,i,(*jjj[i]).dn[n](1),(*jjj[i]).dn[n](2),(*jjj[i]).dn[n](3));
                                   (*jjj[i]).sublattice[n]=1;}
   }
  //check consistency of mcphas.j
  if(nofcomponents!=(*jjj[i]).nofcomponents)
   {fprintf(stderr,"Error loading mcphas.j: the number of spin components for different ions have to be equal\n");
    exit(EXIT_FAILURE);
   }
  
  }
 if  (ferror(fin_coq)==1)
  {fprintf(stderr,"ERROR Reading file %s\n",filejjj);exit(1);}
  fclose (fin_coq);
  fprintf(stderr,"#Finished Reading file %s\n",filejjj);
//myPrintMatrix(stdout,(*(*jjj[1]).G));

}

//kopier-konstruktor
par::par(const par & p)
{ int i;
  a=p.a;b=p.b;c=p.c;alpha=p.alpha;beta=p.beta;gamma=p.gamma;
  r=p.r;rez=p.rez;
  nofatoms=p.nofatoms;
  nofcomponents=p.nofcomponents;
  Cel=p.Cel;CelInv=p.CelInv;
//dimension arrays
  for (i=1;i<=3;++i)
  {rems[i] = new char[strlen(p.rems[i])+2];
   strcpy(rems[i],p.rems[i]);}

 jjj=new jjjpar * [nofatoms+1]; for (i=1;i<=nofatoms;++i){ jjj[i] = new jjjpar(*p.jjj[i]);
  for(int n=1;n<=(*jjj[i]).paranz;++n)
   {(*jjj[i]).sublattice[n]=(*p.jjj[i]).sublattice[n];
 }}
}

//destruktor
par::~par ()
{// printf("hello destruktor par\n"); 
  int i;
  for(i=1;i<=3;++i)
  {delete []rems[i];}
  for(i=1;i<=nofatoms;++i){delete jjj[i];}delete []jjj;
// printf("hello end of destruktor par\n");   
}

int par::newatom(jjjpar * p) //creates new atom from an existing and returns its index
{ jjjpar ** nnn;
  int j;
                  ++nofatoms; // the number of atoms has to be increased
  if(nofatoms>MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL)
  {fprintf(stderr,"ERROR par.cpp: maximum number of atoms in unit cell exceeded - enlarge it in par.hpp and recompile\n");exit(EXIT_FAILURE);}
                  nnn=new jjjpar * [nofatoms+1];
                  for (j=1;j<nofatoms;++j){nnn[j]=jjj[j];} 
                  nnn[nofatoms]=new jjjpar((*p));// use copy constructor to create new atom parameter set    
		  delete []jjj;
		  jjj=nnn;  
return nofatoms;                 
}

int par::delatom(int n) // removes atom number n 
{jjjpar ** nnn;
 int j,again=0;
 FILE * out;
 if(n<1||n>nofatoms){fprintf(stderr,"ERROR par.cpp:delatom n=%i out of range [1:nofatoms=%i]\n",n,nofatoms);exit(EXIT_FAILURE);}
  --nofatoms; // the number of atoms has to be decreased
 nnn=new jjjpar * [nofatoms+1];
 for (j=1;j<n;++j){nnn[j]=jjj[j];
if(0==strcmp((*jjj[n]).sipffilename,(*jjj[j]).sipffilename)){again=1;}} 
 for (j=n+1;j<=nofatoms+1;++j){nnn[j-1]=jjj[j];
if(0==strcmp((*jjj[n]).sipffilename,(*jjj[j]).sipffilename)){again=1;}}

if(again==0){out=fopen("reduce_unitcell_sipf.del","a");
fprintf(out,"%s\n",(*jjj[n]).sipffilename);fclose(out);} 
 delete []jjj;
 jjj=nnn;           
return nofatoms; 
}

void par::reduce_unitcell()
{//checks every atom in the unit cell and removes
// any atom, which is connected to another by a lattice vector
 int i,j,nold=nofatoms;
 Vector d(1,3),n(1,3);
FILE * out;out=fopen("reduce_unitcell_sipf.del","w");fclose(out);

 for(i=1;i<nofatoms;++i){
  for(j=i+1;j<=nofatoms;++j){//printf("nofatoms=%i %i %i\n",nofatoms,i,j);
  d=(*jjj[j]).xyz-(*jjj[i]).xyz;
  n=rez*d;
  if(fabs(rint(n(1))-n(1))<SMALL_MATCH_LATTICEVECTOR&&
     fabs(rint(n(2))-n(2))<SMALL_MATCH_LATTICEVECTOR&&
     fabs(rint(n(3))-n(3))<SMALL_MATCH_LATTICEVECTOR){//printf("del %i\n",j);
                                                      delatom(j);--j;
                                                      }
 }}
 Cel*=(double)nofatoms/nold; // renormalise elastic constants to reduced unit cell dimension
}

void par::add (par & p1)
{int i;
    Vector abc(1,3),abc1(1,3);
    abc(1)=a;abc(2)=b;abc(3)=c;
    abc1(1)=p1.a;abc1(2)=p1.b;abc1(3)=p1.c;
    if (Norm(abc-abc1)>0.0001){fprintf(stderr,"ERROR adding parameter sets: lattice parameters abc not equal\n");exit(EXIT_FAILURE);}
//    if (nofcomponents!=p1.nofcomponents)
//    {fprintf(stderr,"ERROR adding parameter sets: number of spin components not equal\n");exit(EXIT_FAILURE);}
    if (nofcomponents>p1.nofcomponents)
{p1.increase_nofcomponents(nofcomponents-p1.nofcomponents);}
    if (nofcomponents<p1.nofcomponents)
{increase_nofcomponents(p1.nofcomponents-nofcomponents);}

if(p1.nofatoms<nofatoms)
{fprintf(stderr,"# Warning program addj: nofatoms=%i of 1. parameter set greater than %i - continuing, check result with care ! \n",nofatoms,p1.nofatoms);
}
 for(i=1;i<=p1.nofatoms;++i)
 {if (i>nofatoms)
   {newatom(p1.jjj[i]);
   fprintf(stderr,"# Warning program addj: nofatoms=%i not equal %i - adding atom number %i with sipffilename=%s\n",nofatoms,p1.nofatoms,i,(*jjj[i]).sipffilename);

   }
   else
   {if (strcmp((*jjj[i]).sipffilename,(*p1.jjj[i]).sipffilename)!=0){
    fprintf(stderr,"# Warning program addj: adding parameter sets atom %i sipffilename %s does not match %s - taking %s for output\n",i,(*jjj[i]).sipffilename,(*p1.jjj[i]).sipffilename,(*jjj[i]).sipffilename);}
     // add the parameters of p1 to the parameters of this
    (*jjj[i]).add((*p1.jjj[i]),abc);
   }

 }
 Cel+=p1.Cel; // add elastic constants
}

void par::scale(double scalefactor) // scale all interaction parameters by scalefactor
{int i;
 for(i=1;i<=nofatoms;++i)
 {
    (*jjj[i]).scalepars(scalefactor);
 }
}

void par::set_nofcomponents (int n)
{// sets the numberofcomponents to n
 if(n<nofcomponents){decrease_nofcomponents(nofcomponents-n);}
 if(n>nofcomponents){increase_nofcomponents(n-nofcomponents);}
}
void par::increase_nofcomponents (int n)
{//increases the number of components in the interaction vector by n

 int i;
 if (n<1) {fprintf(stderr,"ERROR increasing number of compoments in parameter set: n negative - number cannot be decreased\n");exit(EXIT_FAILURE);}
  fprintf(stderr,"Warning: increasing nofcomponents not tested  yet ... addition of parameter sets may be erroneous\n");

 for(i=1;i<=nofatoms;++i)
 {
    (*jjj[i]).increase_nofcomponents(n);
 }
 nofcomponents+=n;
}

void par::decrease_nofcomponents (int n)
{//decreases the number of components in the interaction vector by n

 int i;
 if (n<1) {fprintf(stderr,"ERROR decreasing number of compoments in parameter set: n negative - number cannot be decreased\n");exit(EXIT_FAILURE);}
 if (nofcomponents-1<n) {fprintf(stderr,"ERROR decreasing number of compoments in parameter set: n = %i must be smaller than nofcomponents = %i\n",n,nofcomponents);exit(EXIT_FAILURE);}
  fprintf(stderr,"Warning: decreasing nofcomponents not tested  yet ... addition of parameter sets may be erroneous\n");

 for(i=1;i<=nofatoms;++i)
 {
    (*jjj[i]).decrease_nofcomponents(n);
 }
 nofcomponents-=n;
}

//save to file
void par::save (const char * filename,int noindexchange)
{ FILE * fout;
  fout = fopen_errchk (filename, "w");
  fprintf(fout,"%s",rems[1]);
  fprintf(fout,"#<!--mcphase.mcphas.j-->\n");
  fprintf(fout,"#***************************************************************\n");
  fprintf(fout,"# Lattice and Exchange Parameter file for\n");
  fprintf(fout,"# %s\n",MCPHASVERSION);
  fprintf(fout,"# - program to calculate static magnetic properties\n");
  fprintf(fout,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
  fprintf(fout,"# %s\n",MCDISPVERSION);
  fprintf(fout,"# - program to calculate the dispersion of magnetic excitations\n");
  fprintf(fout,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
  fprintf(fout,"#***************************************************************\n");
  save (fout,noindexchange);
  fclose(fout);
}

void par::save (FILE * file,int noindexchange)
{ int i;
  errno = 0;
  savelattice(file);
  for (i=1;i<=nofatoms;++i)
  {
  (*jjj[i]).save(file,noindexchange);
//  fprintf(file,"%s",rems[3+i]); // changed 3.03 - I believe here should be only line with stars
    fprintf(file,"#*********************************************************************\n");
  }
 
}
void par::savelattice (FILE *file)
{ 
  errno = 0;
  fprintf(file,"#\n# Lattice Constants (A)\n");
  fprintf(file,"#! a=%4.6g b=%4.6g c=%4.6g alpha=%4.6g beta=%4.6g gamma=%4.6g\n",a,b,c,alpha,beta,gamma);
  fprintf(file,"#! r1a=%4.12g r2a=%4.12g r3a=%4.12g\n",r[1][1],r[1][2],r[1][3]);
  fprintf(file,"#! r1b=%4.12g r2b=%4.12g r3b=%4.12g   primitive lattice vectors [a][b][c]\n",r[2][1],r[2][2],r[2][3]);
  fprintf(file,"#! r1c=%4.12g r2c=%4.12g r3c=%4.12g\n",r[3][1],r[3][2],r[3][3]);

 // save elastic constants
fprintf(file,"#\n# Nonzero Elastic constants   in meV per primitive unit cell in Voigt notation only first index<=second index has to be given\n");
fprintf(file,"# because the constants are symmetric Celij=Celji\n");
fprintf(file,"# Elastic constants refer to the Euclidean coordinate system ijk defined\n");
fprintf(file,"# with respect to abc as j||b, k||(a x b) and i normal to k and j\n");
int i1=0,i,j;fprintf(file,"#! ");
 for(i=1;i<=6;++i){if(i1>0){i1=0;fprintf(file,"#! ");}
  for(j=i;j<=6;++j){if(fabs(Cel(i,j))>SMALL){++i1;fprintf(file," Cel%i%i=%+10.9g",i,j,Cel(i,j));}}
if(i1>0)fprintf(file,"\n"); }
if(i1==0)fprintf(file,"\n");

  fprintf(file,"#\n#! nofatoms=%i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n",nofatoms,nofcomponents);
//  fprintf(file,"%s",rems[3]);// changed 8.09 - I believe here should be only line with stars
  fprintf(file,"#*********************************************************************\n");
}

void par::saveatoms (FILE * file)
{ int i;errno = 0;
  for (i=1;i<=nofatoms;++i)
    {(*jjj[i]).saveatom(file);}
}

void par::save_sipfs(const char *path)   //save single ion parameter files filename to path*
{int i;
 for(i=1;i<=nofatoms;++i)
 {
    (*jjj[i]).save_sipf(path);
 }
}
