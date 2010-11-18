/***********************************************************************
 *
 * jjj2j.c - program to transform the mcphas.jjj format to the mcphas.j format
 *
 ***********************************************************************/

#include <jjj2j.h>
/**********************************************************************/
/**********************************************************************/
// main program
int main (int argc, char **argv)
{
 
// check command line
  if (argc < 2)
    { printf ("use as: jjj2j *.jjj - parameters from file *.jjj\n");
      printf ("are reformatted in .j format to stdout\n");
      printf ("with 'jjj2j mcphas.jjj > mcphas.j' you can create the\n");
      printf ("new input file mcphas.j\n");
      exit (1);
    }

// load  parameters from file  (as object of class parameters)
   par inputpars (argv[1]);
  if (inputpars.nofcomponents!=3)
  {fprintf(stderr,"ERROR jjj2j: more than 3 moment components not supported\n");exit(EXIT_FAILURE);}
  
int i,ii;
for (i=1;i<=inputpars.nofatoms;++i)
 {// calculate additional neighbours for atom number i
  jjjpar addjjj(1,(*inputpars.jjj[i]).diagonalexchange);
  //go throug the list of neighbours of atom number i
  int k,j,jmax,sa,sb,sc,signa,signb,signc;
  Vector d(1,3),d_rint(1,3),xyz(1,3),xyz_rint(1,3);// some vector
  k=1;
  jmax=(*inputpars.jjj[i]).paranz;
  for(j=1;j<=jmax;++j)
   {
    // 1. transform dn(k) to primitive lattice
    xyz=(*inputpars.jjj[i]).dn[k];//new
    d=inputpars.rez*(const Vector&)xyz;
    for (ii=1;ii<=3;++ii)d_rint(ii)=rint(d(ii));
        
       xyz_rint=inputpars.r*(const Vector&)d_rint;
       signa=1;signb=1;signc=1;
       // determination of multiplicity of parameter l
       if (fabs(xyz(1))>0.001&&fabs(fabs(xyz(1))-fabs(xyz_rint(1)))<0.001) 
                               {signa=-1;} // create aditional neighbours at -x(a)
       if (fabs(xyz(2))>0.001&&fabs(fabs(xyz(2))-fabs(xyz_rint(2)))<0.001)
                               {signb=-1;} // create additional neighbours at -y(b)
       if (fabs(xyz(3))>0.001&&fabs(fabs(xyz(3))-fabs(xyz_rint(3)))<0.001) 
                               {signc=-1;} // create additional neighbours at -z(c)
       // go through the different equivalent neighbours
       for (sa=signa;sa<=1;sa+=2){for(sb=signb;sb<=1;sb+=2){for(sc=signc;sc<=1;sc+=2)   
	{if(sa>0&&sb>0&&sc>0){++k;}else
	                               {// create additional neighbour
				        // coordinates
					addjjj.dn[1](1)=(double)sa*(*inputpars.jjj[i]).dn[k](1);//new
				        addjjj.dn[1](2)=(double)sb*(*inputpars.jjj[i]).dn[k](2);//new
				        addjjj.dn[1](3)=(double)sc*(*inputpars.jjj[i]).dn[k](3);//new
	                                // diagonal couplings
					addjjj.jij[1]=(*inputpars.jjj[i]).jij[k];
					// offdiagonal couplings if required
					if((*inputpars.jjj[i]).diagonalexchange==0)
					  { 
					   addjjj.jij[1](1,2)=(double)sa*sb*(*inputpars.jjj[i]).jij[k](1,2);//ab
					   addjjj.jij[1](2,1)=(double)sa*sb*(*inputpars.jjj[i]).jij[k](2,1);//ba
					   addjjj.jij[1](1,3)=(double)sc*sa*(*inputpars.jjj[i]).jij[k](1,3);//ac
					   addjjj.jij[1](3,1)=(double)sc*sa*(*inputpars.jjj[i]).jij[k](3,1);//ca
					   addjjj.jij[1](2,3)=(double)sc*sb*(*inputpars.jjj[i]).jij[k](2,3);//bc
					   addjjj.jij[1](3,2)=(double)sc*sb*(*inputpars.jjj[i]).jij[k](3,2);//cb
                                          }
					(*inputpars.jjj[i]).addpars(k,addjjj);++k;
				       }
        }}}
   }
 }

// save parameters
   inputpars.save(stdout);

   return 0;
}



