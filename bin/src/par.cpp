#include "par.hpp"
#include "../../version"
#include <martin.h>
#include <cstring>

#define MAXNOFNUMBERSINLINE 20
#define MAXNOFCHARINLINE 7024
#define SMALL_MATCH_LATTICEVECTOR 1e-4


 // *************************************************************************
 // ************************ Class for Crystal Structure Parameters *************************************
 // *************************************************************************


//constructor 
par::par(Vector abc,int nofci)
{rez=Matrix(1,3,1,3);
 Cel=Matrix(1,6,1,6);CelInv=Matrix(1,6,1,6);CelInv=0;
 cs.nofcomponents=nofci;
 cs.abc=abc;
  rems[1]=new char[40];rems[1][0]='\0';
  rems[2]=new char[40];rems[2][0]='\0';
  rems[3]=new char[40];rems[3][0]='\0';
  cs.nofatoms=0;
  jjj=new jjjpar * [cs.nofatoms+1];
}

par::par (const char *filejjj,int verbose)
{ int i,j,n,l;
  FILE *fin_coq;
  char instr[MAXNOFCHARINLINE];
  Vector hkl(1,3),hkl_rint(1,3);
  rez=Matrix(1,3,1,3);
  Cel=Matrix(1,6,1,6);Cel=0;CelInv=Matrix(1,6,1,6);CelInv=0;
  fin_coq = fopen_errchk (filejjj, "rb");

 // input file header ------------------------------------------------------------------
  fgets (instr, MAXNOFCHARINLINE, fin_coq);
  extract(instr,"a",cs.abc(1));extract(instr,"b",cs.abc(2)); extract(instr,"c",cs.abc(3)); 
                 extract(instr,"alpha",cs.abc(4));  extract(instr,"beta",cs.abc(5));extract(instr,"gamma",cs.abc(6)); 
  instr[0]='#';
   // inserted 12.11.07 in order to format output correctly (characterstring 13 spoiled output string)
   for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
   rems[1]=new char[strlen(instr)+2];strcpy(rems[1],instr);
   rems[2]=new char[40];strcpy(rems[2],"#!<--mcphas.mcphas.j-->");
  cs.nofatoms=0;
  instr[0]='#';cs.abc=0;
 while (cs.nofatoms==0||(strstr(instr,"*******")==NULL&&instr[strspn(instr," \t")]=='#')) 
  {fgets(instr,MAXNOFCHARINLINE,fin_coq);
   if(Norm(cs.abc)<1e-10){extract(instr,"a",cs.abc(1));extract(instr,"b",cs.abc(2)); extract(instr,"c",cs.abc(3)); 
                 extract(instr,"alpha",cs.abc(4));  extract(instr,"beta",cs.abc(5));extract(instr,"gamma",cs.abc(6)); 
   }
   extract(instr,"r1x",cs.r[1][1]);extract(instr,"r2x",cs.r[1][2]); extract(instr,"r3x",cs.r[1][3]); 
   extract(instr,"r1y",cs.r[2][1]); extract(instr,"r2y",cs.r[2][2]); extract(instr,"r3y",cs.r[2][3]);
   extract(instr,"r1z",cs.r[3][1]); extract(instr,"r2z",cs.r[3][2]); extract(instr,"r3z",cs.r[3][3]);
   extract(instr,"r1a",cs.r[1][1]);extract(instr,"r2a",cs.r[1][2]); extract(instr,"r3a",cs.r[1][3]); 
   extract(instr,"r1b",cs.r[2][1]); extract(instr,"r2b",cs.r[2][2]); extract(instr,"r3b",cs.r[2][3]);
   extract(instr,"r1c",cs.r[3][1]); extract(instr,"r2c",cs.r[3][2]); extract(instr,"r3c",cs.r[3][3]);

    // read optional elastic constants        
  char Celstr[6];
   for(i=1;i<=6;++i)for(j=1;j<=6;++j){
  snprintf(Celstr,sizeof(Celstr),"Cel%i%i",i,j);// printf("%s\n",Celstr);
  extract(instr,Celstr,Cel(i,j));Cel(j,i)=Cel(i,j);}

   extract(instr,"nofatoms",cs.nofatoms);extract(instr,"nofcomponents",cs.nofcomponents); 
		  if(feof(fin_coq)!=0)
                    {fprintf(stderr,"ERROR reading header of file %s: line '#! nofatoms=...' not found\n",filejjj);exit(EXIT_FAILURE);}
  }
  if(cs.nofatoms>MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL)
  {fprintf(stderr,"ERROR reading mcphas.j: maximum number of atoms in unit cell exceeded - enlarge it in par.hpp and recompile\n");exit(EXIT_FAILURE);}
  // check if primitive lattice is right handed
 if(Det(cs.r)<0){fprintf(stderr,"ERROR reading mcphas.j: primitive lattice r1 r2 r3 not right handed\n");exit(EXIT_FAILURE);}
 
  rez=cs.r.Inverse();
  rems[3]=new char[strlen(instr)+2];strcpy(rems[3],instr);
  
  //read parameter sets for every atom 
  jjj=new jjjpar * [cs.nofatoms+1];
  //gJ=Vector(1,cs.nofatoms);
  for(i=1;i<=cs.nofatoms;++i)  
  {//printf("creating atom %i (of %i)...\n",i,cs.nofatoms);
   jjj[i]=new jjjpar(fin_coq,cs.nofcomponents,verbose);
   if(jjj[i]==NULL){ fprintf (stderr, "Out of memory creating atoms jjjpar by constructor\n");exit (EXIT_FAILURE);}
   //gJ(i)=(*jjj[i]).gJ;
   cs.sipffilenames[i]=(*jjj[i]).sipffilename;
   cs.x[i]=(*jjj[i]).xyz[1];
   cs.y[i]=(*jjj[i]).xyz[2];
   cs.z[i]=(*jjj[i]).xyz[3];
   if(cs.nofcomponents!=(*jjj[i]).nofcomponents)
   {fprintf(stderr,"ERROR reading mcphas.j: nofcomponents (%i) not consistent for atom %i (%i read in fileheader)\n",(*jjj[i]).nofcomponents,i,cs.nofcomponents);exit(EXIT_FAILURE);}
  }
  //determine sublattices and rij 
  for(i=1;i<=cs.nofatoms;++i)
  {for(n=1;n<=(*jjj[i]).paranz;++n)
   {(*jjj[i]).sublattice[n]=0;
    for(j=1;j<=cs.nofatoms;++j)
    {//try if neighbour n is on sublattice j
     hkl=rez*((*jjj[i]).xyz+(*jjj[i]).dn[n]-(*jjj[j]).xyz);
     // check if hkl is integer - if yes then the neighbour n is on sublattice j
     for(l=1;l<=3;++l){hkl_rint(l)=rint(hkl(l));}
     if(Norm(hkl_rint-hkl)<0.001){(*jjj[i]).sublattice[n]=j;}
    }
    if((*jjj[i]).sublattice[n]==0){fprintf(stderr,"Warning mcphas - par.cpp: file %s inconsistent:  neighbour %i of atom %i at %g %g %g is not on any sublattice. Continuing putting it onto sublattice 1 ...\n",filejjj,n,i,(*jjj[i]).dn[n](1),(*jjj[i]).dn[n](2),(*jjj[i]).dn[n](3));
                                   (*jjj[i]).sublattice[n]=1;}
   // transform (*jjj[i]).dn[n] to Euclidean frame and store in dr[n]
  dadbdc2ijk((*jjj[i]).dr[n],(*jjj[i]).dn[n], cs.abc);

   }
  //check consistency of mcphas.j
  if(cs.nofcomponents!=(*jjj[i]).nofcomponents)
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

//kopier-konstruktor copy constructor
par::par(const par & p)
{ int i;
  cs.abc=p.cs.abc;
  cs.r=p.cs.r;rez=p.rez;
  cs.nofatoms=p.cs.nofatoms;
  cs.nofcomponents=p.cs.nofcomponents;
  Cel=p.Cel;CelInv=p.CelInv;
//dimension arrays
  for (i=1;i<=3;++i)
  {rems[i] = new char[strlen(p.rems[i])+2];
   strcpy(rems[i],p.rems[i]);}
 jjj=new jjjpar * [cs.nofatoms+1]; for (i=1;i<=cs.nofatoms;++i){ jjj[i] = new jjjpar(*p.jjj[i]);
 
  cs.sipffilenames[i]=(*jjj[i]).sipffilename;
  cs.x[i]=(*jjj[i]).xyz[1];
  cs.y[i]=(*jjj[i]).xyz[2];
  cs.z[i]=(*jjj[i]).xyz[3];

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
  for(i=1;i<=cs.nofatoms;++i){delete jjj[i];}delete []jjj;
// printf("hello end of destruktor par\n");   
}

int par::newatom(jjjpar * p) //creates new atom from an existing and returns its index
{ jjjpar ** nnn;
  int j;
                  ++cs.nofatoms; // the number of atoms has to be increased
  if(cs.nofatoms>MAX_NOF_ATOMS_IN_PRIMITIVE_CRYST_UNITCELL)
  {fprintf(stderr,"ERROR par.cpp: maximum number of atoms in unit cell exceeded - enlarge it in par.hpp and recompile\n");exit(EXIT_FAILURE);}
                  nnn=new jjjpar * [cs.nofatoms+1];
                  for (j=1;j<cs.nofatoms;++j){nnn[j]=jjj[j];} 
                  nnn[cs.nofatoms]=new jjjpar((*p));// use copy constructor to create new atom parameter set 
                  cs.sipffilenames[cs.nofatoms]=(*nnn[cs.nofatoms]).sipffilename;
                  cs.x[cs.nofatoms]=(*nnn[cs.nofatoms]).xyz[1];
                  cs.y[cs.nofatoms]=(*nnn[cs.nofatoms]).xyz[2];
                  cs.z[cs.nofatoms]=(*nnn[cs.nofatoms]).xyz[3];

		  delete []jjj;
		  jjj=nnn;  
return cs.nofatoms;                 
}

int par::delatom(int nn, Matrix & distribute,int verbose) // removes atom number n 
// if n<0 then atom number |n| is removed and also all interactions of other atoms
// with this atom are removed from the interaction table 
// if n>0 interactions with the other atoms are kept and transferred to 
// a group of atoms (numbers given in column 1 of distribute) with 
// coefficients given in col 2 of distribute. Only interactions
// with atoms given in column 1 of distribute are removed completely.
// Attention: sublattice index is not changed by this function !
//            -- thus after running it sublattice[s] refers still to
//            original numbering with all atoms in the parameters set !
{jjjpar ** nnn;
 int j,s,again=0,n=nn; if(nn<0)n=-nn;
 FILE * out;
 if(n<1||n>cs.nofatoms){fprintf(stderr,"ERROR par.cpp:delatom n=%i out of range [1:nofatoms=%i]\n",n,cs.nofatoms);exit(EXIT_FAILURE);}
  --cs.nofatoms; // the number of atoms has to be decreased
 nnn=new jjjpar * [cs.nofatoms+1];
if(verbose)fprintf(stderr,"Deleting atom %i\n",n);
 for (j=1;j<=cs.nofatoms+1;++j){if(j==n)++j;
if(verbose)fprintf(stderr,"caring about interactions of atom %i\n",j);
// take care of all interactions to be removed because atom is removed
 for(s=1;s<=(*jjj[j]).paranz;++s){
     if((*jjj[j]).sublattice[s]==n){if(nn>0){// transfer interactions to the remaining atoms
                                             // unless atom j is to be distributed on
                                             int j_in_distribute=0;
                                             for(int i=1;i<=distribute.Rhi();++i)
                                               {if(j==(int)distribute(i,1))j_in_distribute=1;
                                               }
                                             if(j_in_distribute==0){// ok j is not in list .. transfer interaction to those in the list
                                                   for(int i=1;i<=distribute.Rhi();++i)
                                               {int sl=(int)distribute(i,1);int mult=0;double Rmax=1e10;
                                                double coeff=distribute(i,2);
						for(int ss=1;ss<=(*jjj[j]).paranz;++ss){
                                                if((*jjj[j]).sublattice[ss]==sl){double R=Norm((*jjj[j]).dr[ss]-(*jjj[j]).dr[s]);
                                                         if(R<Rmax+SMALL_MATCH_LATTICEVECTOR){++mult;}
                                                         if(R<Rmax-SMALL_MATCH_LATTICEVECTOR){mult=1;Rmax=R;}
                                                                                 }
                                                                                       }
if(mult==0){fprintf(stderr,"Error reduce_unitcell on distributing interaction %i of atom %i onto atom %i: atom %i has no neighbour on sublattice %i\n",s,j,sl,j,sl);exit(1);}
if(verbose)fprintf(stderr,"distributing interaction %i of atom %i onto %i atoms on sublattice %i with coefficient %g \n",s,j,mult,sl,coeff);
						for(int ss=1;ss<=(*jjj[j]).paranz;++ss){
                                                if((*jjj[j]).sublattice[ss]==sl){double R=Norm((*jjj[j]).dr[ss]-(*jjj[j]).dr[s]);
                                                    if(R<Rmax+SMALL_MATCH_LATTICEVECTOR){
                                                       (*jjj[j]).jij[ss]+=(*jjj[j]).jij[s]*(coeff/mult);int found=0;
for(int si=1;si<=(*jjj[sl]).paranz;++si){// look for corresponding neighbour in the interaction sl and also change to preserve Symmetry of interaction table
if(Norm((*jjj[sl]).dr[si]+(*jjj[j]).dr[ss])<SMALL_MATCH_LATTICEVECTOR){(*jjj[sl]).jij[si]+=(*jjj[j]).jij[s].Transpose()*(coeff/mult);
++found;}                                       }
if(found!=1){fprintf(stderr,"Error reduce_unitcell on distributing interaction %i of atom %i onto atom %i: no corresponding bond found from atom %i \n",s,j,sl,sl);exit(1);}

                                                                                        }
                                                                                 }
                                                                                       }

                                               }
                                                                   }
                                            }
                                    (*jjj[j]).delpar(s);--s;}
          //fprintf(stderr,"%i neighbour of ion %i paranz=%i sublattice=%i \n",s,j,(*nnn[j]).paranz,(*nnn[j]).sublattice[s]);
          }
   if(j>n){                
 cs.sipffilenames[j-1]=(*jjj[j]).sipffilename;
 cs.x[j-1]=(*jjj[j]).xyz[1];
 cs.y[j-1]=(*jjj[j]).xyz[2];
 cs.z[j-1]=(*jjj[j]).xyz[3];
          }
if(0==strcmp((*jjj[n]).sipffilename,(*jjj[j]).sipffilename)){again=1;}
              if(j+1==n)++j;}

for (j=1;j<n;++j){nnn[j]=jjj[j];}
for (j=n+1;j<=cs.nofatoms+1;++j){nnn[j-1]=jjj[j];}
// correct the sublattice numbering 
for (j=1;j<=cs.nofatoms;++j){
for(s=1;s<=(*nnn[j]).paranz;++s){if((*nnn[j]).sublattice[s]>n)--(*nnn[j]).sublattice[s];}
                             }

if(again==0){out=fopen("reduce_unitcell_sipf.del","a");
fprintf(out,"%s\n",(*jjj[n]).sipffilename);fclose(out);} 
 delete []jjj;
 jjj=nnn;           
return cs.nofatoms; 
}


void par::reduce_unitcell(int verbose)
{//checks every atom in the unit cell and removes
// any atom, which is connected to another by a lattice vector
 int i,j,nold=cs.nofatoms;
 Vector d(1,3),n(1,3);Matrix dis(1,1,1,1);dis=0;
FILE * out;out=fopen("reduce_unitcell_sipf.del","w");fclose(out);

 for(i=1;i<cs.nofatoms;++i){int ct=0;
  for(j=i+1;j<=cs.nofatoms;++j){//printf("nofatoms=%i %i %i\n",cs.nofatoms,i,j);
  d=(*jjj[j]).xyz-(*jjj[i]).xyz;
  n=rez*d;
  if(fabs(rint(n(1))-n(1))<SMALL_MATCH_LATTICEVECTOR&&
     fabs(rint(n(2))-n(2))<SMALL_MATCH_LATTICEVECTOR&&
     fabs(rint(n(3))-n(3))<SMALL_MATCH_LATTICEVECTOR){//printf("del %i\n",j);
                                                      delatom(-j,dis,verbose);--j;++ct;
                                                      }
    }
  if(verbose){fprintf(stderr,"For atom %i there have been deleted %i equivalent atoms\n",i,ct);}
    // check if atom is indeed in primitive unit cell - if not move it there
   n=rez*(*jjj[i]).xyz;
   while(n(1)>=1.0)--n(1);
   while(n(2)>=1.0)--n(2);
   while(n(3)>=1.0)--n(3);
   while(n(1)<0.0)++n(1);
   while(n(2)<0.0)++n(2);
   while(n(3)<0.0)++n(3);
   (*jjj[i]).xyz=cs.r*n;
 }
 Cel*=(double)cs.nofatoms/nold; // renormalise elastic constants to reduced unit cell dimension
}

void par::add (par & p1)
{int i,i1;
    if (Norm(cs.abc-p1.cs.abc)>0.0001){fprintf(stderr,"ERROR adding parameter sets: lattice parameters abc not equal\n");exit(EXIT_FAILURE);}
//    if (cs.nofcomponents!=p1.cs.nofcomponents)
//    {fprintf(stderr,"ERROR adding parameter sets: number of spin components not equal\n");exit(EXIT_FAILURE);}
    if (cs.nofcomponents>p1.cs.nofcomponents)
{p1.increase_nofcomponents(cs.nofcomponents-p1.cs.nofcomponents);}
    if (cs.nofcomponents<p1.cs.nofcomponents)
{increase_nofcomponents(p1.cs.nofcomponents-cs.nofcomponents);}

if(p1.cs.nofatoms<cs.nofatoms)
{fprintf(stderr,"# Warning program addj: nofatoms=%i of 1. parameter set greater than %i - continuing, check result with care ! \n",cs.nofatoms,p1.cs.nofatoms);
}
 for(i=1;i<=p1.cs.nofatoms;++i)
 {if (i>cs.nofatoms)
   {newatom(p1.jjj[i]);
   fprintf(stderr,"# Warning program addj: nofatoms=%i not equal %i - adding atom number %i with sipffilename=%s\n",cs.nofatoms-1,p1.cs.nofatoms,i,(*jjj[i]).sipffilename);

   }
   else
   {
i1=i;
while(Norm((*jjj[i1]).xyz-(*p1.jjj[i]).xyz)>0.0001)
{if(i1==i)fprintf(stderr,"# Problem program addj: atomic positions da db dc of atom number %i  at %g %g %g do not match atom %i at %g %g %g \n# ... trying to find a matching atom\n",
i1,(*jjj[i1]).xyz(1),(*jjj[i1]).xyz(2),(*jjj[i1]).xyz(3),i,(*p1.jjj[i]).xyz(1),(*p1.jjj[i]).xyz(2),(*p1.jjj[i]).xyz(3));
if(i1==cs.nofatoms)i1=0;
++i1;
if(i1==i){fprintf(stderr,"# ... no matching atom found: adding new atom number %i with sipffilename=%s\n",cs.nofatoms+1,(*p1.jjj[i]).sipffilename);
newatom(p1.jjj[i]);i1=cs.nofatoms;
(*jjj[i1]).scalepars(0.0); // set jjjpars to zero - because there comes the add command in line 242
 }
}
    if (strcmp((*jjj[i1]).sipffilename,(*p1.jjj[i]).sipffilename)!=0){
    fprintf(stderr,"# Warning program addj: adding parameter sets atom %i sipffilename %s does not match %s - taking %s for output\n",
i,(*jjj[i]).sipffilename,(*p1.jjj[i]).sipffilename,(*jjj[i]).sipffilename);}
     // add the parameters of p1 to the parameters of this
    (*jjj[i1]).add((*p1.jjj[i]),cs.abc);
     }
  

 }
 Cel+=p1.Cel; // add elastic constants
}

void par::scale(double scalefactor) // scale all interaction parameters by scalefactor
{int i;
 for(i=1;i<=cs.nofatoms;++i)
 {
    (*jjj[i]).scalepars(scalefactor);
 }
}

void par::set_nofcomponents (int n)
{// sets the numberofcomponents to n
 if(n<cs.nofcomponents){decrease_nofcomponents(cs.nofcomponents-n);}
 if(n>cs.nofcomponents){increase_nofcomponents(n-cs.nofcomponents);}
}
void par::increase_nofcomponents (int n)
{//increases the number of components in the interaction vector by n

 int i;
 if (n<1) {fprintf(stderr,"ERROR increasing number of compoments in parameter set: n negative - number cannot be decreased\n");exit(EXIT_FAILURE);}
  fprintf(stderr,"Warning: increasing nofcomponents not tested  yet ... addition of parameter sets may be erroneous\n");

 for(i=1;i<=cs.nofatoms;++i)
 {
    (*jjj[i]).increase_nofcomponents(n);
 }
 cs.nofcomponents+=n;
}

void par::decrease_nofcomponents (int n)
{//decreases the number of components in the interaction vector by n

 int i;
 if (n<1) {fprintf(stderr,"ERROR decreasing number of compoments in parameter set: n negative - number cannot be decreased\n");exit(EXIT_FAILURE);}
 if (cs.nofcomponents-1<n) {fprintf(stderr,"ERROR decreasing number of compoments in parameter set: n = %i must be smaller than nofcomponents = %i\n",n,cs.nofcomponents);exit(EXIT_FAILURE);}
  fprintf(stderr,"Warning: decreasing nofcomponents not tested  yet ... addition of parameter sets may be erroneous\n");

 for(i=1;i<=cs.nofatoms;++i)
 {
    (*jjj[i]).decrease_nofcomponents(n);
 }
 cs.nofcomponents-=n;
}

//save to file
void par::save (const char * filename,int noindexchange)
{ FILE * fout;
  fout = fopen_errchk (filename, "w");
  save (fout,noindexchange);
  fclose(fout);
}



void par::save (FILE * file,int noindexchange)
{ int i;
  errno = 0;
  fprintf(file,"%s",rems[1]);
//  fprintf(file,"#<!--mcphase.mcphas.j-->\n");
  fprintf(file,"#***************************************************************\n");
  fprintf(file,"# Lattice and Exchange Parameter file for\n");
  fprintf(file,"# %s\n",MCPHASVERSION);
  fprintf(file,"# - program to calculate static magnetic properties\n");
  fprintf(file,"# reference: M. Rotter JMMM 272-276 (2004) 481\n");
  fprintf(file,"# %s\n",MCDISPVERSION);
  fprintf(file,"# - program to calculate the dispersion of magnetic excitations\n");
  fprintf(file,"# reference: M. Rotter et al. J. Appl. Phys. A74 (2002) 5751\n");
  fprintf(file,"#***************************************************************\n");
  savelattice(file);
  for (i=1;i<=cs.nofatoms;++i)
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
  fprintf(file,"#! a=%4.6g b=%4.6g c=%4.6g alpha=%4.6g beta=%4.6g gamma=%4.6g\n",cs.abc(1),cs.abc(2),cs.abc(3),cs.abc(4),cs.abc(5),cs.abc(6));
  fprintf(file,"#! r1a=%4.12g r2a=%4.12g r3a=%4.12g\n",cs.r[1][1],cs.r[1][2],cs.r[1][3]);
  fprintf(file,"#! r1b=%4.12g r2b=%4.12g r3b=%4.12g   primitive lattice vectors [a][b][c]\n",cs.r[2][1],cs.r[2][2],cs.r[2][3]);
  fprintf(file,"#! r1c=%4.12g r2c=%4.12g r3c=%4.12g\n",cs.r[3][1],cs.r[3][2],cs.r[3][3]);

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

  fprintf(file,"#\n#! nofatoms=%i  nofcomponents=%i  number of atoms in primitive unit cell/number of components of each spin\n",cs.nofatoms,cs.nofcomponents);
//  fprintf(file,"%s",rems[3]);// changed 8.09 - I believe here should be only line with stars
  fprintf(file,"#*********************************************************************\n");
}

void par::saveatoms (FILE * file)
{ int i;errno = 0;
  for (i=1;i<=cs.nofatoms;++i)
    {(*jjj[i]).saveatom(file);}
}

void par::save_sipfs(const char *path)   //save single ion parameter files filename to path*
{int i;
 for(i=1;i<=cs.nofatoms;++i)
 {
    (*jjj[i]).save_sipf(path);
 }
}

void par::save_mcdiff_in (const char * program)
{FILE * fout;Vector r(1,3),d(1,3);
  fout = fopen_errchk ("mcdiff.in", "w");
// if atom is not magnetic print it out here
int natcryst=0;for(int i=1;i<=cs.nofatoms;++i)if((*jjj[i]).magnetic==0)++natcryst;
  cs.print_mcdiff_in_header(fout,program,natcryst);

for(int i=1;i<=cs.nofatoms;++i)
{if((*jjj[i]).magnetic==0)
 {d(1)=cs.x[i];d(2)=cs.y[i];d(3)=cs.z[i];
  r=rez*d;
  fprintf(fout,"%+10.8f  %+10.8f         %+10.8f %+10.8f %+10.8f    %+10.8f %+10.8f %+10.8f  %+10.8f  # %s\n",
  (*jjj[i]).SLR,(*jjj[i]).SLI,cs.x[i],cs.y[i],cs.z[i],r(1),r(2),r(3),(*jjj[i]).DWF,cs.sipffilenames[i]);
 }

}

fprintf(fout,"\
#\n\
#\n\
# %%SECTION 3%% DESCRIPTION OF THE LATTICE\n\
#\n\
# -----------------------------------------------------------------------------\n");
savelattice(fout);
int nat=0;for(int i=1;i<=cs.nofatoms;++i)if((*jjj[i]).magnetic!=0)++nat;
fprintf(fout,"\
#\n\
# %%SECTION 4%% DESCRIPTION OF MAGNETIC UNIT CELL AND LIST OF MAGNETIC ATOMS\n\
#\n\
#\n\
# here follows the description of the magnetic unit cell with respect\n\
# to the primitive crystallographic unit cell:\n\
# 'nr1', 'nr2', 'nr3' ...the crystallographic unit cell has to be taken\n\
#                        nr1 nr2 and nr3 times along r1 r2 and r3,\n\
#                        respectively to get magnetic unit cell\n\
# 'nat' denotes the number of magnetic atoms in magnetic unit cell\n\
#\n\
# Temperature,  External Magnetic Field: Magnetic Unit Cell\n\
#! T=1 K Ha=0 T Hb= 0 T Hc= 0 T: nr1=1 nr2=1 nr3=1 nat=%i \n\
#\n\
#\n\
# It follows a list of nat lines with to describe the magnetic moment configuration\n\
# Notes:\n\
# 'atom-filename' means the single ion property filename of this magnetic atom:\n\
#                 -it must contain the Formfactor Coefficients (e.g. see international tables)\n\
#                                      Lande factor\n\
#                                      Neutron Scattering Length (10^-12 cm) \n\
#                 -it may contain a    Debey Waller Factor\n\
# 'da' 'db' and 'dc' are not used by the program (unless you enter a line #! use_dadbdc=1)\n\
# 'dr1','dr2' and 'dr3' refer to the primitive lattice given below\n\
# 'Ma','Mb','Mc' denote the magnetic moment components in Bohr magnetons\n\
#                in case of non orthogonal lattices instead of Ma Mb Mc the components Mx My Mz\n\
#                have to be given, which refer to an right handed orthogonal coordinate system \n\
#                defined by y||b, z||(a x b) and x normal to y and z\n\
#  <Sa> <Sb> <Sc>  <La> <Lb > <Lc>  (optional) denote the spin and orbital angular momentum components \n\
# 'mf1' 'mf2' 'mf3' (optional line, used to go beyond dipole approx for formfactor)\n\
#                                     denote the corresponding exchange fields in meV\n\
#\n\
#{atom-file} da[a]  db[b]    dc[c]     dr1[r1]  dr2[r2]  dr3[r3]   <Ma>     <Mb>     <Mc> [mb] [optional <Sa> <Sb> <Sc> <La> <Lb> <Lc> ]\n\
#{corresponding exchange fields [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)}\n",nat);

fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [Moment created by program %s]\n",program);
for(int i=1;i<=cs.nofatoms;++i)
{if((*jjj[i]).magnetic!=0)
 {d(1)=cs.x[i];d(2)=cs.y[i];d(3)=cs.z[i];
  r=rez.Transpose()*d;
  fprintf(fout,"{%s}         %+10.8f %+10.8f %+10.8f    %+10.8f %+10.8f %+10.8f 0 0 0 \n",
  cs.sipffilenames[i],cs.x[i],cs.y[i],cs.z[i],r(1),r(2),r(3));
 }

}


  fclose(fout);
 
}

// operator !=  returns 8 7 6 5 4 3 2 1 0depending on agreement of
 //  8 abc 7 nofatoms 6 atomic positions 5 sipffilenames 4 nofcomponents 3 nofneighbours disagreement
 //  2 neighbour position 1 interaction parmeter disagreement i.e. 0 is perfect match
int par::operator!= (par & op2) // match
{if(Norm(cs.abc-op2.cs.abc)>0.0001){fprintf(stderr,"# Lattice does not match\n");return 8;}
 if(cs.nofatoms!=op2.cs.nofatoms){fprintf(stderr,"# Number of atoms nofatoms do not match\n");return 7;}
 for(int i=1;i<=cs.nofatoms;++i)
  {if(Norm((*op2.jjj[i]).xyz-(*jjj[i]).xyz)>0.0001){fprintf(stderr,"# Atomic position of atom %i do not match\n",i);return 6;}
 }
for(int i=1;i<=cs.nofatoms;++i)
{
if (strcmp((*jjj[i]).sipffilename,(*op2.jjj[i]).sipffilename)!=0)
  {fprintf(stderr,"# Sipffilename %s and %s of atom %i do not match\n",(*jjj[i]).sipffilename,(*op2.jjj[i]).sipffilename,i);return 5;}

  }
 if(cs.nofcomponents!=op2.cs.nofcomponents){fprintf(stderr,"# nofcomponents do not match\n");return 4;}

for(int i=1;i<=cs.nofatoms;++i)
  {if((*op2.jjj[i]).paranz!=(*jjj[i]).paranz){fprintf(stderr,"# nofneighbours of atom %i do not match\n",i);return 3;}
 }

for(int i=1;i<=cs.nofatoms;++i)
 {for(int n=1;n<=(*op2.jjj[i]).paranz;++n)
  if(Norm((*jjj[i]).dn[n]-(*op2.jjj[i]).dn[n])>0.001){fprintf(stderr,"# neighbour position %i of atom %i do not match\n",n,i);return 2;}
 }
for(int i=1;i<=cs.nofatoms;++i)
 {for(int n=1;n<=(*op2.jjj[i]).paranz;++n)
  if(NormFro((*jjj[i]).jij[n]-(*op2.jjj[i]).jij[n])>0.001){fprintf(stderr,"# interactions %i of atom %i do not match\n",n,i);return 1;}
 }
return 0;
}