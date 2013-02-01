/***********************************************************
 *
 * chrgplt - program to display charge density of 
 *           an ion given its CF pars, T and effective H
 * Author: M. Rotter
 **********************************************************/


#include "../../version"
#include "density.hpp"
#include "spincf.hpp"
#include "jjjpar.hpp"
#include "martin.h"
#include "plt_func.c"
#include "cryststruct.hpp"

/**********************************************************************/
// main program
int main (int argc, char **argv)
{
printf("***********************************************************\n");
printf("*\n");
printf("* chrgplt - program to display charge density of \n");
printf("*           an ion given its CF pars, T and effective H\n");
printf("* Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("* %s\n",MCPHASVERSION);
printf("***********************************************************\n");
// check command line
  if (argc < 7)
    { printf ("\nprogram chrgplt - calculate chargdensity of a single ion\n\n\
                use as: chrgplt threshhold T Ha Hb Hc  mcphas.sipf \n\n\
                - given is temperature T[K] and magnetic effective field H[T]\n \
		- crystal field  parameters Blm should be read from a \n\
	           standard mcphas single ion property file mcphas.sipf \n\
                options: if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered\n\
                by hand.\n\
		\n");
      printf ("\n\
                ... to view chargeplots type: javaview results/chrgplt.jvx\n\
                                              displaycontour 2 3 4 results/chrgplti.grid\n\
                                              displaycontour 1 3 4 results/chrgpltj.grid\n\
                                              displaycontour 1 2 4 results/chrgpltk.grid\n");
      exit (1);
    }

 graphic_parameters gp;
gp.show_abc_unitcell=0;
gp.show_primitive_crystal_unitcell=0;
gp.show_magnetic_unitcell=0;
gp.threshhold=strtod(argv[1],NULL);
gp.scale_pointcharges=1;
gp.show_pointcharges=1;
  gp.show_atoms=1;gp.showprim=1;
  gp.spins_scale_moment=1;
 gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
cryststruct cs;

  double T,ha,hb,hc;
  T=strtod(argv[2],NULL);
  ha=strtod(argv[3],NULL);
  hb=strtod(argv[4],NULL);
  hc=strtod(argv[5],NULL);

 FILE * fout;
 // read cf-parameters into class object jjjpar
jjjpar jjjps(0.0,0.0,0.0,argv[6]);
  int dim;
  double lnz,u;
  if(jjjps.gJ==0)
  {dim=51;} // external module ic1ion
  else
  {dim=48;}

  if(jjjps.module_type!=0&&jjjps.module_type!=2&&jjjps.module_type!=4)
  {fprintf(stderr,"ERROR chrgplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}

  if (jjjps.module_type==0&&jjjps.gJ!=0)
  {fprintf(stderr,"************** WARNING **********************\n reading external single ion module with gJ not zero: gJ=%g - please check if calculation of density is supported !!!\n*********************************************\n",jjjps.gJ);}

Vector h(1,dim);
  Vector moments(1,dim);
  h=0;

fout = fopen_errchk ("results/chrgplt.coeff", "w");
fprintf(fout,"# coefficients for density calculation\n#T=%g K field H=(%g,%g,%g) Tesla\n",T,ha,hb,hc);

printf("# T=%g K field H=(%g,%g,%g) Tesla\n",T,ha,hb,hc);
  if(jjjps.gJ==0){ h(1)=2.0*MU_B*ha;h(3)=2.0*MU_B*hb;h(5)=2.0*MU_B*hc;h(2)=MU_B*ha;h(4)=MU_B*hb;h(6)=MU_B*hc;}
  else { h(1)=jjjps.gJ*MU_B*ha;h(2)=jjjps.gJ*MU_B*hb;h(3)=jjjps.gJ*MU_B*hc;}
  //int dj=(int)(2.0*(*iops).J+1);

  jjjps.mcalc_parameter_storage_init(h,T);
  printf("calculating expectation values ....\n");
  jjjps.mcalc(moments,T,h,lnz,u,jjjps.mcalc_parstorage);
//  cfield  has to be used to calculate all the <Olm>.
//  printf("Stevens factors: alpha beta gamma = %4g %4g %4g \n",(*iops).alpha,(*iops).beta,(*iops).gamma);
//  printf("Lande Factor: gJ = %4g\n",(*iops).gJ);

int i,nofpc=0;
const char lm[]=             "Jy  Jz  Jx  O22SO21SO20 O21 O22 O33SO32SO31SO30 O31 O32 O33 O44SO43SO42SO41SO40 O41 O42 O43 O44 O55SO54SO53SO52SO51SO50 O51 O52 O53 O54 O55 O66SO65SO64SO63SO62SO61SO60 O61 O62 O63 O64 O65 O66 ";
const char lme[]="Sx  Lx  Sy  Ly  Sz  Lz  O22SO21SO20 O21 O22 O33SO32SO31SO30 O31 O32 O33 O44SO43SO42SO41SO40 O41 O42 O43 O44 O55SO54SO53SO52SO51SO50 O51 O52 O53 O54 O55 O66SO65SO64SO63SO62SO61SO60 O61 O62 O63 O64 O65 O66 ";
char lm4[5];lm4[4]='\0';
printf("#chargedensity is expanded in tesseral harmonics as\n#   ro(r)= -|e||R4f(r)|^2 sum_lm clm Zlm(Omega)\n#\n");
fprintf(fout,"#chargedensity is expanded in tesseral harmonics as\n#   ro(r)= -|e||R4f(r)|^2 sum_lm clm Zlm(Omega)\n#\n");
   for(i=1;i<=dim;++i){int l,m;double factor;
                       if (dim==48) strncpy(lm4,lm+(i-1)*4,4);
                       if (dim==51) strncpy(lm4,lme+(i-1)*4,4);
                       factor=0;
                  // the following prefactors should correspond exactly to jjjpar_basmodfunc.cpp lines 933 ff !!!
                       if(jjjps.module_type==2||jjjps.module_type==4){if(i>3){l=lm4[1]-48;m=lm4[2]-48;
                                                        factor=jjjps.tetan()(l)*jjjps.cnst(l,m);
                                                        }
                                               } // these are prefactors in case of module cfield (stevens parameters tetan and zlm prefactors)
                                      else     {if(i>6){l=lm4[1]-48;m=lm4[2]-48;
                                                       if(m!=0){factor=sqrt((2.0*l+1)/8/PI);}
                                                           else{factor=sqrt((2.0*l+1)/4/PI);}
                                                       }
                                               } 
                       printf(" <J%c> = <%s> =%12.6f   clm=%12.6f\n",'a'-1+i,lm4,moments(i),moments(i)*factor);
                       fprintf(fout," <J%c> = <%s> =%12.6f   clm=%12.6f\n",'a'-1+i,lm4,myround(moments(i)),myround(moments(i)*factor));}
printf("\n");
  fclose (fout);

 char text[1000];
 if(jjjps.module_type==0||jjjps.module_type==4){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=abc</title>\n", T, ha, hb, hc);}
 if(jjjps.module_type==2){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=cab</title>\n", T, ha, hb, hc);}

 cs.cffilenames[1]=new char[MAXNOFCHARINLINE];
 cs.abc(1)=6.0;cs.abc(2)=6.0;cs.abc(3)=6.0;
 cs.r=0;cs.r(1,1)=1.0;cs.r(2,2)=1.0;cs.r(3,3)=1.0;
 strcpy(cs.cffilenames[1],argv[6]);
 cs.x[1]=0.5*gp.scale_view_1;cs.y[1]=0.5*gp.scale_view_2;cs.z[1]=0.5*gp.scale_view_3; // put atom in middle of cell

// read pointcharge-parameters
if(gp.show_pointcharges>0) nofpc=read_pointcharge_parameters(gp,cs.cffilenames,argv[6],cs.x,cs.y,cs.z,jjjps,cs.abc);
  spincf s(1,1,1,nofpc+1,dim);
 // Vector gJJ(1,nofpc+1);gJJ=0;gJJ(1)=jjjps.gJ;
  cs.gJ[1]=jjjps.gJ;
  Vector hkl(1,3);hkl=0;s=s*0;
  spincf ev_real(s),ev_imag(s);

  for(i=1;i<=dim;++i)s.m(1,1,1)(i)=moments(i);
  fout = fopen_errchk ("results/chrgplt.jvx", "w");
   s.jvx_cd(fout,text,cs,gp,0.0,ev_real,ev_imag,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/chrgplt.grid", "w");
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/chrgplti.grid", "w");
  gp.gridi=1;gp.gridj=200;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/chrgpltj.grid", "w");
  gp.gridi=200;gp.gridj=1;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/chrgpltk.grid", "w");
  gp.gridi=200;gp.gridj=200;gp.gridk=1;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program chrgplt\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/chrgplt.jvx\n");
fprintf(stderr,"# * displaycontour 2 3 4 results/chrgplti.grid\n");
fprintf(stderr,"# * displaycontour 1 3 4 results/chrgpltj.grid\n");
fprintf(stderr,"# * displaycontour 1 2 4 results/chrgpltk.grid\n");
fprintf(stderr,"# * saved density mesh in results/chrgplt.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  int dj;
  for(dj=1;dj<=nofpc+1;++dj){delete cs.cffilenames[dj];}
  return 0;
}


