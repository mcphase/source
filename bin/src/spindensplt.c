/***********************************************************
 *
 * spindensplt - program to display spin density of
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
printf("* spindensplt - program to display spin density of \n");
printf("*               an ion given its CF pars, T and effective H\n");
printf("* Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("* %s\n",MCPHASVERSION);
printf("***********************************************************\n");
// check command line
  if (argc < 7)
    { printf ("\nprogram spindensplt - calculate spindensity of a single ion\n\n\
                use as: spindensplt threshhold T Ha Hb Hc i j k mcphas.sipf \n\n\
                use or: spindensplt threshhold T Ha Hb Hc mcphas.sipf \n\n\
                use or: spindensplt threshhold T Ha Hb Hc -div mcphas.sipf \n\n\
                - given is temperature T[K] and magnetic effective field H[T]\n\
		- the spindensity vector component along direction (i,j,k) is\n\
                  calculated (if omitted absolute value of spind. is calc.)\n\
                - crystal field  parameters Blm should be read from a \n\
	          standard mcphas single ion property file mcphas.sipf \n\
                options: \n\
                -div triggers calculation of divergence of the vector field\n\
                - if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered \n\
                by hand.\n\
		\n");
      printf ("\n\
...to view spindensityplots type: javaview results/spindensplt.jvx\n\
                                  displaycontour 2 3 4 results/spindenplti.grid\n\
                                  displaycontour 1 3 4 results/spindenpltj.grid\n\
                                  displaycontour 1 2 4 results/spindenpltk.grid\n\
                \n");
      exit (1);
    }

  double T,ha,hb,hc,xx=0,yy=0,zz=0;


  T=strtod(argv[2],NULL);
  ha=strtod(argv[3],NULL);
  hb=strtod(argv[4],NULL);
  hc=strtod(argv[5],NULL);
int doijk=0;
if (argc>9){
  xx=strtod(argv[6],NULL);
  yy=strtod(argv[7],NULL);
  zz=strtod(argv[8],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(xx*xx+yy*yy+zz*zz);
  xx/=rr;yy/=rr;zz/=rr;
  doijk=3;
 }

if (argc>7&&strncmp(argv[6],"-div",4)==0)
{doijk=1;
}

 jjjpar jjjps(0.0,0.0,0.0,argv[6+doijk]);


  if(jjjps.module_type!=0){fprintf(stderr,"ERROR spindensplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}
    int dim=49; if(doijk<3){dim=3*49;}

  if (jjjps.module_type==0&&jjjps.gJ!=0)
  {fprintf(stderr,"************** WARNING **********************\n reading external single ion module with gJ not zero: gJ=%g - please check if calculation of density is supported !!!\n*********************************************\n",jjjps.gJ);}



// Indices for spindensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};



  Vector h(1,6);
  Vector moms(1,6);
  Vector moments(1,dim);
  Vector momentsx(1,49);
  Vector momentsy(1,49);
  Vector momentsz(1,49);
  h=0;
 FILE * fout;


 if(jjjps.gJ==0){ h(1)=2.0*MU_B*ha;h(3)=2.0*MU_B*hb;h(5)=2.0*MU_B*hc;h(2)=MU_B*ha;h(4)=MU_B*hb;h(6)=MU_B*hc;}
  else { h(1)=jjjps.gJ*MU_B*ha;h(2)=jjjps.gJ*MU_B*hb;h(3)=jjjps.gJ*MU_B*hc;}
   //int dj=(int)(2.0*(*iops).J+1);

  jjjps.mcalc_parameter_storage_init(h,T);


  printf("calculating expectation values ....\n");
//double lnz,u;
  //jjjps.mcalc(moms,T,h,lnz,u,jjjps.mcalc_parstorage);
if(xx!=0||doijk<3)jjjps.spindensity_mcalc (momentsx,1, T, h, jjjps.mcalc_parstorage);
if(yy!=0||doijk<3)jjjps.spindensity_mcalc (momentsy,2, T, h, jjjps.mcalc_parstorage);
if(zz!=0||doijk<3)jjjps.spindensity_mcalc (momentsz,3, T, h, jjjps.mcalc_parstorage);
if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}

int i,nofpc=0;
fout = fopen_errchk ("results/spindensplt.coeff", "w");
fprintf(fout,"# coefficients for density calculation\n#T=%g K field H=(%g,%g,%g) Tesla\n",T,ha,hb,hc);
printf("# T=%g K field H=(%g,%g,%g) Tesla\n",T,ha,hb,hc);
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   Ms(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
fprintf(fout,"#spindensity is expanded in tesseral harmonics Zlm\n\
#   Ms(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
   for(i=1;i<=49;++i){
  if(doijk==3){             printf(" aS(%i,%i) =%12.6f\n",k[i],q[i],myround(moments(i)));}
  else{printf(" aSx(%i,%i) =%12.6f",k[i],q[i],myround(momentsx(i)));
       printf(" aSy(%i,%i) =%12.6f",k[i],q[i],myround(momentsy(i)));
       printf(" aSz(%i,%i) =%12.6f\n",k[i],q[i],myround(momentsz(i)));
       fprintf(fout," aSx(%i,%i) =%12.6f",k[i],q[i],myround(momentsx(i)));
       fprintf(fout," aSy(%i,%i) =%12.6f",k[i],q[i],myround(momentsy(i)));
       fprintf(fout," aSz(%i,%i) =%12.6f\n",k[i],q[i],myround(momentsz(i)));
       moments(i)=momentsx(i);moments(i+49)=momentsy(i);moments(i+2*49)=momentsz(i);
      }

}
printf("\n");
fclose(fout);
graphic_parameters gp;
gp.show_abc_unitcell=0;
gp.show_primitive_crystal_unitcell=0;
gp.show_magnetic_unitcell=0;
gp.threshhold=strtod(argv[1],NULL);
gp.scale_pointcharges=1;
gp.show_pointcharges=1;
gp.scale_density_vectors=1;
 gp.show_atoms=1;gp.showprim=1;
  gp.spins_scale_moment=0;
  gp.spins_show_static_moment_direction=0;
if(gp.read())printf("#reading graphic parameters from results/graphic_parameters.set\n");// read graphic parameters which are set by user in file results/graphic_parameters.set
 gp.spins_scale_moment=0;
 gp.spins_show_static_moment_direction=0;
 gp.spins_wave_amplitude=0;
 gp.spins_show_oscillation=0;
 gp.spins_show_ellipses=0;

// read cf-parameters into class object jjjpar
cryststruct cs;

 char text[1000];
 if(jjjps.module_type==0||jjjps.module_type==4){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=abc, spindensity Ms(r).(%g,%g,%g)</title>\n", T, ha, hb, hc,xx,yy,zz);}
 if(jjjps.module_type==2){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=cab, spindensity Ms(r).(%g,%g,%g)</title>\n", T, ha, hb, hc,xx,yy,zz);}

 cs.cffilenames[1]=new char[MAXNOFCHARINLINE];
 cs.abc(1)=6.0;cs.abc(2)=6.0;cs.abc(3)=6.0;
 cs.r=0;cs.r(1,1)=1.0;cs.r(2,2)=1.0;cs.r(3,3)=1.0;
 strcpy(cs.cffilenames[1],argv[6+doijk]);
 cs.x[1]=0.5*gp.scale_view_1;cs.y[1]=0.5*gp.scale_view_2;cs.z[1]=0.5*gp.scale_view_3; // put atom in middle of cell


// read pointcharge-parameters
if(gp.show_pointcharges>0) nofpc=read_pointcharge_parameters(gp,cs.cffilenames,argv[6+doijk],cs.x,cs.y,cs.z,jjjps,cs.abc);

  spincf s(1,1,1,nofpc+1,dim);
  Vector hkl(1,3);hkl=0;s=s*0;
  spincf ev_real(s),ev_imag(s);

  for(i=1;i<=dim;++i)s.m(1,1,1)(i)=moments(i);
  cs.gJ[1]=jjjps.gJ;
  if(doijk==3) sprintf(gp.title,"projection of spindensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of spindensity div Ms(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of spindensity |Ms(r)|");
  printf("%s\n",gp.title);
  fout = fopen_errchk ("results/spindensplt.jvx", "w");
   s.jvx_cd(fout,text,cs,gp,0.0,ev_real,ev_imag,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/spindensplt.grid", "w");
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/spindensplti.grid", "w");
  gp.gridi=1;gp.gridj=200;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/spindenspltj.grid", "w");
  gp.gridi=200;gp.gridj=1;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
  fout = fopen_errchk ("results/spindenspltk.grid", "w");
  gp.gridi=200;gp.gridj=200;gp.gridk=1;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h);
  fclose (fout);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spindensplt\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/spindensplt.jvx\n");
fprintf(stderr,"# * displaycontour 2 3 4 results/spindensplti.grid\n");
fprintf(stderr,"# * displaycontour 1 3 4 results/spindenspltj.grid\n");
fprintf(stderr,"# * displaycontour 1 2 4 results/spindenspltk.grid\n");
fprintf(stderr,"# * saved density mesh in results/spindensplt.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  int dj;
  for(dj=1;dj<=nofpc+1;++dj){delete cs.cffilenames[dj];}
  return 0;
}


