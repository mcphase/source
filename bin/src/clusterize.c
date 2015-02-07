/**************************************************************************
 *
 * clusterize - program to read mcphas.j and create clusterized input files
 *
 **************************************************************************/
#include "par.hpp"
#include <complex>
#include "martin.h"

#include <cerrno>
#include <sys/types.h>
#include <sys/stat.h>
#include <vector.h>

// main program
int main (int argc, char **argv)
{printf("# Program Clusterize \n");
 if(argc<2){printf ("\n\
Creates groups of atoms to form clusters. Each cluster is treated \n\
as a single subsystem and exactly diagonalised. The subsystems are coupled \n\
and treated in the standard mean field and dynamical matrix diagonalisation (DMD) \n\
procedures in mcphas and mcdisp. Such coupled cluster calculations are \n\
appropriate for solids with groups of magnetic ions which interact strongly and \n\
different groups are coupled weakly. For an example see Cu2Te2O5Cl2, Jens Jensen \n\
PRB 2008 \n\
\n\n\
- use as:\n\n\
 clusterize mcphas_in.j 1 2 3 0 4 5 6 0 7 8 9 \n \
 ... creates mcphas.j cluster1.sipf cluster1.j cluster2.sipf cluster2.j ... \n \
    from input file mcphas_in.j. Groups of atoms have to be separated by a zero.\n \
		\n");exit(1);}
 // 1st load parameters as they are
 par inp(argv[1]);
 printf("# program clusterize: finished reading input file\n");
 // determine the number of components (=orig nofcomp x number of atoms in biggest cluster)
 int i=2,cnr,noc=inp.nofcomponents;
 while(i<argc)
  {int n=0;
   while(i<argc&&argv[i][0]!='0')
   {n+=inp.nofcomponents;++i;
   }++i;
   if(n>noc)noc=n;
  }
 printf("# ... clusterizing and testing output\n");
 par out(inp.a,inp.b,inp.c,inp.alpha,inp.beta,inp.gamma,noc);
 out.r=inp.r;out.rez=inp.rez;
 char filename[1000];
 // create cluster.j files
 i=2;cnr=0;
 while(i<argc)
 {// for each cluster create a par and store it to cluster.j
  par clust(inp.a,inp.b,inp.c,inp.alpha,inp.beta,inp.gamma,inp.nofcomponents);
  clust.r=2.0*inp.r;clust.rez=0.5*inp.rez;
  ++cnr;int iatom1=i,nofcat=0;Vector cog(1,3);cog=0;
  while(i<argc&&argv[i][0]!='0')
  {int n=(int)strtod(argv[i],NULL);
   if(n>inp.nofatoms||n<1){fprintf(stderr,"Error clusterize: atom number out of range 1-%i\n",inp.nofatoms);exit(1);}
   int newat=clust.newatom(inp.jjj[n]);cog+=(*inp.jjj[n]).xyz;++nofcat;
   ++i;
   // here delete neighbours which do not belong to cluster
   for(int nn=1;nn<=(*clust.jjj[newat]).paranz;++nn)// loop through neighbours of the new atom
   {int i1=iatom1,keep=0;
    while(i1<argc&&argv[i1][0]!='0')// loop through atoms of cluster
    {int n1=(int)strtod(argv[i1],NULL);
     // check if neighbor is one of the atoms
     if((*clust.jjj[newat]).sublattice[nn]==n1)
      {// test also if distance is really that of atom n1 in first unit cell 
       if(Norm((*clust.jjj[newat]).dn[nn]-(*inp.jjj[n1]).xyz+(*clust.jjj[newat]).xyz)<0.001)
       {keep=1;}
      }++i1;
    }
    if(keep==0){// delete neighbour from neighbour list 
                (*clust.jjj[newat]).delpar(nn);--nn;
               }
   }
  }++i;
  cog/=nofcat; // now center of gravity of cluster is determined
  for(int n=1;n<=clust.nofatoms;++n)(*clust.jjj[n]).xyz-=cog;

  // store cluster'cnr'.j
  sprintf(filename,"cluster%i.j",cnr);
  clust.save(filename);
  sprintf(filename,"cluster%i.sipf",cnr);
  FILE * fout;
  fout =fopen_errchk (filename, "w");
  fprintf(fout,"\
#!MODULE=cluster\n\
#<!--mcphase.sipf-->\n");
  fprintf(fout,"#!structurefile=cluster%i.j\n",cnr);
  fprintf(fout,"\
# ... the noperl option allows only simple matrix operations\n\
#     but saves a lot of memory ...\n\
#!noperl\n\
# the individual moments\n");
  for(int n=1;n<=clust.nofatoms;++n)
  {fprintf(fout,"\
$M%i_1=2.0*$I%i_1;\n\
$M%i_2=2.0*$I%i_2;\n\
$M%i_3=2.0*$I%i_3;\n",n,n,n,n,n,n);
  }
  fprintf(fout,"\
# the total moment\n\
$M1=$M1_1;\n\
$M2=$M1_2;\n\
$M3=$M1_3;\n");
  for(int n=2;n<=clust.nofatoms;++n)
  {fprintf(fout,"\
$M1=$M1+$M%i_1;\n\
$M2=$M2+$M%i_2;\n\
$M3=$M3+$M%i_3;\n",n,n,n);
   }
   fprintf(fout,"\
# for running mcdisp it is needed to define a interaction operator\n");
   for(int n=1;n<=clust.nofatoms;++n)for(int nn=1;nn<=inp.nofcomponents;++nn)
   {fprintf(fout,"$I%i=$I%i_%i;\n",(n-1)*inp.nofcomponents+nn,n,nn);}
   for(int n=clust.nofatoms*inp.nofcomponents+1;n<=noc;++n)
   {fprintf(fout,"$I%i=$I1_1;\n",n);}    
  fclose(fout);
  jjjpar jjjnew(cog(1),cog(2),cog(3),filename,noc);
  jjjnew.diagonalexchange=0;
  out.newatom(&jjjnew);
 }

    // first: increase the number of components in the in file so that we can store
     // these atoms correctly
     int nofcomponents=inp.nofcomponents;
     inp.increase_nofcomponents(noc-nofcomponents);
 

  // here we need to go through the atoms of all of the clusters and check all their neighbours
  i=2;cnr=0;while(i<argc)
  {int iatom1=i,catom1=0;++cnr;
   while(i<argc&&argv[i][0]!='0')
   {int n=(int)strtod(argv[i],NULL);++i;++catom1;
    for(int nnn=1;nnn<=(*inp.jjj[n]).paranz;++nnn)
    {
  // ... if
  // 1) a neighbour is within the cluster, do not store it in the nn list
     int i1=iatom1,keep=1;
     while(i1<argc&&argv[i1][0]!='0')// loop through atoms of cluster
     {int n1=(int)strtod(argv[i1],NULL);
      // check if neighbor is one of the atoms
      if((*inp.jjj[n]).sublattice[nnn]==n1)
      {// test also if distance is really that of atom n1 in first unit cell 
       if(Norm((*inp.jjj[n]).dn[nnn]-(*inp.jjj[n1]).xyz+(*inp.jjj[n]).xyz)<0.001)
       {keep=0;}
      }++i1;
     } //printf("atom %i neighbour %i keep %i\n",n,nnn,keep);
    if(keep==1) // hm atom is not within the cluster 
     { // check if it is within any cluster
       i1=2;int in_a_cluster=0,cnr1=0,catom2=0;
       while(i1<argc)
       {++cnr1;int catom2i=0;
        while(i1<argc&&argv[i1][0]!='0')// loop through atoms of cluster
        {int n1=(int)strtod(argv[i1],NULL);++catom2i;
         // check if neighbor is one of the atoms
         if((*inp.jjj[n]).sublattice[nnn]==n1){in_a_cluster=cnr1;catom2=catom2i;}
         ++i1;
        }++i1;
       }
       if(in_a_cluster==0)
       {// 2) a neighbour is not within any cluster - add interaction to interaction list taking
       //    care of atomic position - cog of cluster 
                     // determine vector of atom to center of cluster
                       Vector dd(1,3);
                         dd=(*out.jjj[cnr]).xyz-(*inp.jjj[n]).xyz;
                     // subtract this vector from the neighbour coordinate
                         (*inp.jjj[n]).dn[nnn]-=dd;
       //and mind to put nonzero the appropriate off diagonal component
       // ... modify (*inp.jjj[n]).jij[nnn] 
                     // put the interaction constants to the appropriate degree of freedom
                     // clear interaction constants which have been automatically set
                      Matrix jij(1,noc,1,noc);jij=0;
                     // interaction constants are in (*inp.jjj[n]).jij[nnn] ... put them into the correct places
                      for(int nx=1;nx<=nofcomponents;++nx)for(int ny=1;ny<=nofcomponents;++ny)
                       jij(nx+nofcomponents*(catom1-1),ny)=
                            (*inp.jjj[n]).jij[nnn](nx,ny);
                     // put the correctly distributed matrix to inp.jjj
                       (*inp.jjj[n]).jij[nnn]=jij;                 
       //printf("not in any cluster");
       }else{//if(in_a_cluster==cnr){fprintf(stderr,"Error: atom %i neighbour %i is in same cluster but keep=1\n",n,nnn);exit(EXIT_FAILURE);}
       // 3) a neighbour is within another cluster: take distance to cog of that cluster 
                   // determine vector of atom to center of cluster
                       Vector dd(1,3);
                         dd=(*out.jjj[cnr]).xyz-(*inp.jjj[n]).xyz;
                     // subtract this vector from the neighbour coordinate
                         (*inp.jjj[n]).dn[nnn]-=dd;
                   // determine vector of atom to center of cluster of the other atom
                         dd=(*out.jjj[in_a_cluster]).xyz-(*inp.jjj[(*inp.jjj[n]).sublattice[nnn]]).xyz;
                     // add this vector to the neighbour coordinate
                         (*inp.jjj[n]).dn[nnn]+=dd;
        // ...and store
       //    interaction minding to put nonzero the appropriate off diagonal component 
       // ... modify (*inp.jjj[n]).jij[nnn] !!! 
              // put the interaction constants to the appropriate degree of freedom
                     // clear interaction constants which have been automatically set
                      Matrix jij(1,noc,1,noc);jij=0;
                     // interaction constants are in (*inp.jjj[n]).jij[nnn] ... put them into the correct places
                      for(int nx=1;nx<=nofcomponents;++nx)for(int ny=1;ny<=nofcomponents;++ny)
                       jij(nx+nofcomponents*(catom1-1),ny+nofcomponents*(catom2-1))=
                            (*inp.jjj[n]).jij[nnn](nx,ny);
                     // put the correctly distributed matrix to inp.jjj
                       (*inp.jjj[n]).jij[nnn]=jij;  
       //     printf("in  cluster %i\n",in_a_cluster);
       }
     } // fi keep ==1
     else
     {(*inp.jjj[n]).delpar(nnn);--nnn;}
    }
    Vector abc(1,3);
    abc(1)=out.a;abc(2)=out.b;abc(3)=out.c;
    (*out.jjj[cnr]).add((*inp.jjj[n]),abc);
   }++i;
  }
  
   
 // here add the remaining atoms (which are not in clusters) to out !
 for(int nn=1;nn<=inp.nofatoms;++nn)
 {// check if atom is in cluster
  int inclust=0;
  i=2;while(i<argc)
  {while(i<argc&&argv[i][0]!='0')
   {int n=(int)strtod(argv[i],NULL);++i;
    if (n==nn)inclust=1;
   }++i;
  }
  if(inclust==0){// add atom to out
                 jjjpar jjjnew((*inp.jjj[nn]).xyz(1),(*inp.jjj[nn]).xyz(2),(*inp.jjj[nn]).xyz(3),(*inp.jjj[nn]).sipffilename,noc);
                 int nt=out.newatom(&jjjnew);
 
                 //int nt=out.newatom(inp.jjj[nn]); // !!! still needs to take care of neihgbours !!!
                      // loop through neighbours of newatom and check if some are in clusters
                      // if not - leave them as they are
                      // if yes, set newatoms diagonalexchange =0 and determine distance to coq of corresponding cluster and
                      // put interaction to (off diagonal) component corresponding to the appropriate 
                      // degree of freedom of that cluster
                      for(int nnn=1;nnn<=(*inp.jjj[nn]).paranz;++nnn)
                   { inclust=0;cnr=1;int catom=0,catomi=0;
                     i=2;while(i<argc)
                     {while(i<argc&&argv[i][0]!='0')
                      {int n=(int)strtod(argv[i],NULL);++i;++catomi;
                       if (n==(*inp.jjj[nn]).sublattice[nnn]){inclust=cnr;catom=catomi;}
                      }++i;++cnr;catomi=0;
                     }
                     if(inclust>0){//printf("hello incluster\n");
                          (*out.jjj[nt]).diagonalexchange=0;
                   // determine distance to cluster coq 
                     // determine vector of atom to center of cluster
                       Vector dd(1,3);
                         dd=(*out.jjj[inclust]).xyz-(*inp.jjj[(*inp.jjj[nn]).sublattice[nnn]]).xyz;
                     // add this vector to the neighbour coordinate
                         (*inp.jjj[nn]).dn[nnn]+=dd;
                  // put the interaction constants to the appropriate degree of freedom
                     // clear interaction constants which have been automatically set
                      Matrix jij(1,noc,1,noc);jij=0;
                     // interaction constants are in (*inp.jjj[nn]).jij[nnn] ... put them into the correct places
                      for(int nx=1;nx<=nofcomponents;++nx)for(int ny=1;ny<=nofcomponents;++ny)
                       jij(nx,ny+nofcomponents*(catom-1))=
                            (*inp.jjj[nn]).jij[nnn](nx,ny);
                     // put the correctly distributed matrix to inp.jjj
                       (*inp.jjj[nn]).jij[nnn]=jij;
                                   }
                   }
                  Vector abc(1,3);
                  abc(1)=out.a;abc(2)=out.b;abc(3)=out.c;
                  (*out.jjj[nt]).add((*inp.jjj[nn]),abc);
                 }
 }
 out.save("mcphas.j");
}