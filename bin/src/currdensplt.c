/***********************************************************
 *
 * currdensplt - program to display current density of
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
printf("* currdensplt - program to display current density of \n");
printf("*               an ion given its CF pars, T and effective H\n");
printf("* Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("* %s\n",MCPHASVERSION);
printf("***********************************************************\n");
// check command line
  if (argc < 7)
    { printf ("\nprogram currdensplt - calculate current density of a single ion\n\n\
                use as: currdensplt threshhold T Ha Hb Hc mcphas.sipf \n\n\
                    or: currdensplt threshhold T Ha Hb Hc i j k mcphas.sipf \n\n\
                    or: currdensplt threshhold T Ha Hb Hc -div mcphas.sipf \n\n\
                - given is temperature T[K] and magnetic effective field H[T]\n\
		- the current density vector component along direction (i,j,k)\n\
                  is calculated (if omitted abs value of currdens is calc.)\n\
                - crystal field  parameters Blm should be read from a \n\
	          standard mcphas single ion property file mcphas.sipf \n\
                options: \n\
                -div triggers calculation of divergence of the vector field\n\
                - if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered \n\
                by hand.\n\
		\n");
      printf ("\n\
...to view currdensityplots type: javaview results/currdensplt.jvx\n\
                                  displaycontour 2 3 4 results/currdensplti.grid\n\
                                  displaycontour 1 3 4 results/currdenspltj.grid\n\
                                  displaycontour 1 2 4 results/currdenspltk.grid\n\
                \n");
      exit (1);
    }


  double T,xx=0,yy=0,zz=0;
  Vector Hext(1,3);
  T=strtod(argv[2],NULL);
  Hext(1)=strtod(argv[3],NULL);
  Hext(2)=strtod(argv[4],NULL);
  Hext(3)=strtod(argv[5],NULL);
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
 FILE * fout;

 // read cf-parameters into class object jjjpar
 jjjpar jjjps(0.0,0.0,0.0,argv[6+doijk]);
 if(jjjps.module_type!=0){fprintf(stderr,"ERROR currdensplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}
  if (jjjps.module_type==0&&jjjps.gJ!=0)
  {fprintf(stderr,"************** WARNING **********************\n reading external single ion module with gJ not zero: gJ=%g - please check if calculation of density is supported !!!\n*********************************************\n",jjjps.gJ);}

  
  int dim=3*49;
// Indices for currdensity
//          0 not used
//          0 1  2 3 4  5  6 7 8  9 10 11 1213141516 17 18 192021222324 25 26 27 28 29303132333435 36 37 38 39 40 414243444546474849
int k[] = {-1,0, 1,1,1, 2, 2,2,2,2, 3, 3, 3,3,3,3,3, 4, 4, 4, 4,4,4,4,4,4, 5, 5, 5, 5, 5,5,5,5,5,5,5, 6, 6, 6, 6, 6, 6,6,6,6,6,6,6,6};
int q[] = {-1,0,-1,0,1,-2,-1,0,1,2,-3,-2,-1,0,1,2,3,-4,-3,-2,-1,0,1,2,3,4,-5,-4,-3,-2,-1,0,1,2,3,4,5,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6};

  Vector h(1,6);
  Vector momS(1,6);
  Vector moments(1,2*dim);


  
  Vector momentlx(1,49);
  Vector momently(1,49);
  Vector momentlz(1,49);
  h=0;

  //int dj=(int)(2.0*(*iops).J+1);

  jjjps.Icalc_parameter_storage_init(h,Hext,T);
  printf("calculating expectation values ....\n");
//double lnz,u;
 // jjjps.Icalc(momS,T,h,lnz,u,jjjps.Icalc_parstorage);
jjjps.orbmomdensity_coeff (momentlx,1, T, h,Hext, jjjps.Icalc_parstorage);
jjjps.orbmomdensity_coeff (momently,2, T, h,Hext, jjjps.Icalc_parstorage);
jjjps.orbmomdensity_coeff (momentlz,3, T, h,Hext, jjjps.Icalc_parstorage);
int i,nofpc=0;
fout = fopen_errchk ("results/currdensplt.coeff", "w");
fprintf(fout,"# coefficients for density calculation\n#T=%g K field H=(%g,%g,%g) Tesla\n",T,Hext(1),Hext(2),Hext(3));
printf("# T=%g K field H=(%g,%g,%g) Tesla\n",T,Hext(1),Hext(2),Hext(3));
printf("#currdensity is expanded in tesseral harmonics Zlm\n\
#   j(r).(%g,%g,%g)= sum_lm (b(l,m) R^2(r)+ d(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
fprintf(fout,"# T=%g K field H=(%g,%g,%g) Tesla\n",T,Hext(1),Hext(2),Hext(3));
fprintf(fout,"#currdensity is expanded in tesseral harmonics Zlm\n\
#   j(r).(%g,%g,%g)= sum_lm (b(l,m) R^2(r)+ d(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
   for(i=1;i<=49;++i){
                       printf(" aLx(%i,%i) =%12.6f  ",k[i],q[i],myround(momentlx(i)));
                       printf(" aLy(%i,%i) =%12.6f  ",k[i],q[i],myround(momently(i)));
                       printf(" aLz(%i,%i) =%12.6f\n",k[i],q[i],myround(momentlz(i)));
                       fprintf(fout," aLx(%i,%i) =%12.6f  ",k[i],q[i],myround(momentlx(i)));
                       fprintf(fout," aLy(%i,%i) =%12.6f  ",k[i],q[i],myround(momently(i)));
                       fprintf(fout," aLz(%i,%i) =%12.6f\n",k[i],q[i],myround(momentlz(i)));
             moments(i)=momentlx(i);moments(i+49)=momently(i);moments(i+2*49)=momentlz(i);
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
 gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
 gp.spins_scale_moment=0;
 gp.spins_show_static_moment_direction=0;
 gp.spins_wave_amplitude=0;
 gp.spins_show_oscillation=0;
 gp.spins_show_ellipses=0;
cryststruct cs;

 char text[1000];
 if(jjjps.module_type==0||jjjps.module_type==4){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=abc, currdensity j(r).(%g,%g,%g)(milliAmp/A^2)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);}
 if(jjjps.module_type==2){sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates xyz=cab, currdensity j(r).(%g,%g,%g)(milliAmp/A^2)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);}

 cs.cffilenames[1]=new char[MAXNOFCHARINLINE];
 cs.abc(1)=6.0;cs.abc(2)=6.0;cs.abc(3)=6.0;
 cs.r=0;cs.r(1,1)=1.0;cs.r(2,2)=1.0;cs.r(3,3)=1.0;
 strcpy(cs.cffilenames[1],argv[6+doijk]);
 cs.x[1]=0.5*gp.scale_view_1;cs.y[1]=0.5*gp.scale_view_2;cs.z[1]=0.5*gp.scale_view_3; // put atom in middle of cell


// read pointcharge-parameters
if(gp.show_pointcharges>0) nofpc=read_pointcharge_parameters(gp,cs.cffilenames,argv[6+doijk],cs.x,cs.y,cs.z,jjjps,cs.abc);

  spincf s(1,1,1,nofpc+1,2*dim);
  cs.gJ[1]=jjjps.gJ;
  Vector hkl(1,3);hkl=0;s=s*0;
  spincf ev_real(s),ev_imag(s);

  for(i=1;i<=2*dim;++i)s.m(1,1,1)(i)=moments(i);
  gp.show_atoms=1;gp.showprim=1;
  gp.spins_scale_moment=0;
  gp.spins_show_static_moment_direction=0;
  if(doijk==3) sprintf(gp.title,"projection of currdensity j(r).(i=%g,j=%g,k=%g)(milliAmp/A^2)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of currdensity div j(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of currdensity |j(r)|(milliAmp/A^2)");
  printf("%s\n",gp.title);
fout = fopen_errchk ("results/currdensplt.jvx", "w");
   s.jvx_cd(fout,text,cs,gp,0.0,ev_real,ev_imag,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/currdensplt.grid", "w");
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/currdensplti.grid", "w");
  gp.showprim=1;
  gp.gridi=1;gp.gridj=200;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/currdenspltj.grid", "w");
  gp.gridi=200;gp.gridj=1;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/currdenspltk.grid", "w");
  gp.gridi=200;gp.gridj=200;gp.gridk=1;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program currdensplt\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/currdensplt.jvx\n");
fprintf(stderr,"# * displaycontour 2 3 4 results/currdensplti.grid\n");
fprintf(stderr,"# * displaycontour 1 3 4 results/currdenspltj.grid\n");
fprintf(stderr,"# * displaycontour 1 2 4 results/currdenspltk.grid\n");
fprintf(stderr,"# * saved density mesh in results/currdensplt.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  int dj;
  for(dj=1;dj<=nofpc+1;++dj){delete cs.cffilenames[dj];}
  return 0;
}


