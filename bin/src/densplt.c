/***********************************************************
 *
 * densplt - program to display charge density of 
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

void help_and_exit()
    { printf ("\nprogram densplt - calculate density of a single ion\n\n\
                use as: densplt c|s|o|m|j [-p i j k|-div] [-S|-L|-M] mcphas.sipf T Ha Hb Hc\n\n\
                 c ... calculate chargedensity\n\
                 s ... calculate spindensity\n\
                 o ... calculate angular orbital momentum density\n\
                 m ... calculate magnetic moment density\n\
                 j ... calculate currentdensity\n\
        optional p i j k ... calculate projection of spin/orbital/current/magnetic moment density along direction i j k, e.g. 0 0 1\n\
        optional -div    ... calculate divergence of spin/orbital/current/magnetic moment density  \n\
        optional -S  ... show arrow indicating spin\n\
        optional -L  ... show arrow indicating orbital angular momentum\n\
        optional -M  ... show arrow indicating magnetic moment\n\
                 - crystal field  parameters Blm should be read from a \n\
	           standard mcphas single ion property file mcphas.sipf \n\
                 - given is temperature T[K] and magnetic effective field H[T]\n \
		options: if T<0 then no thermal boltzmann distribution is taken\n\
		the statistical probability of each CF state has to be entered\n\
                by hand.\n\
     example:\n\
               densplt c Pr3p.sipf 2 0 0 1\n\
             ...calculates the charge density using crystal field from Pr3p.sipf\n\
                at T=2K and H=(0,0,1) Tesla\n \
		\n");
      printf ("\n\
                ... to view chargeplots type: javaview results/densplt.jvx\n\
                                              displaycontour 2 3 4 results/densplti.grid\n\
                                              displaycontour 1 3 4 results/denspltj.grid\n\
                                              displaycontour 1 2 4 results/denspltk.grid\n");
      exit (1);
    }
/**********************************************************************/
// main program
int main (int argc, char **argv)
{
printf("***********************************************************\n");
printf("*\n");
printf("* densplt - program to display charge density of \n");
printf("*           an ion given its CF pars, T and effective H\n");
printf("* Reference: M. Rotter PRB 79 (2009) 140405R\n");
printf("* %s\n",MCPHASVERSION);
printf("***********************************************************\n");
// check command line
if (argc < 7){help_and_exit();}

graphic_parameters gp;
gp.show_abc_unitcell=0;
gp.show_primitive_crystal_unitcell=0;
gp.show_magnetic_unitcell=0;
gp.show_pointcharges=1;
gp.show_atoms=1;gp.showprim=1;
gp.scale_pointcharges=1;
gp.scale_density_vectors=1;
gp.spins_scale_moment=1;
gp.spins_wave_amplitude=0;
gp.spins_show_static_moment_direction=0;
gp.spins_show_oscillation=0;
gp.spins_show_ellipses=0;

cryststruct cs;
int os=0; int doijk=0,arrow=0,arrowdim=3;
 double T,xx=0,yy=0,zz=0;
 
if(strcmp(argv[2],"-div")==0){os=1;doijk=1;}
if(strcmp(argv[2],"-p")==0){os=4;
  xx=strtod(argv[3],NULL);
  yy=strtod(argv[4],NULL);
  zz=strtod(argv[5],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(xx*xx+yy*yy+zz*zz);
  xx/=rr;yy/=rr;zz/=rr;
  doijk=3;
                        }

if(strcmp(argv[2+os],"-S")==0){os+=1;arrow=1;arrowdim=SPIN_EV_DIM;gp.spins_colour=3;}
if(strcmp(argv[2+os],"-L")==0){os+=1;arrow=2;arrowdim=ORBMOM_EV_DIM;gp.spins_colour=2;}
if(strcmp(argv[2+os],"-M")==0){os+=1;arrow=3;arrowdim=MAGMOM_EV_DIM;gp.spins_colour=1;}

  // read cf-parameters into class object jjjpar
  jjjpar jjjps(0.0,0.0,0.0,argv[2+os],1);
  Vector Hext(1,3);
  T=strtod(argv[3+os],NULL);
  Hext(1)=strtod(argv[4+os],NULL);
  Hext(2)=strtod(argv[5+os],NULL);
  Hext(3)=strtod(argv[6+os],NULL);

 FILE * fout;
  int dim;
  dim=28;
 char text[1000];
 char coord[]="xyz=abc";
if(jjjps.module_type==2){strcpy(coord,"xyz=cba");}

switch(argv[1][0]) // dimension definition from jjjpar.hpp
{case 'c': dim=CHARGEDENS_EV_DIM;
printf("#chargedensity is expanded in tesseral harmonics Zlm\n\
#   ro(r) sum_lm (a(l,m) R^2(r) Zlm(Omega)\n\
#   M. Rotter et al. J Phys: Conf Ser. 325 (2011) 012005\n#\n ");
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates %s, chargedensity ro(r)</title>\n", T,Hext(1),Hext(2),Hext(3),coord);
 sprintf(gp.title,"chargedensity ro(r)");
 gp.threshhold=-0.05;
           break;
 case 's': dim=SPINDENS_EV_DIM;
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates %s, spindensity S(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),coord,xx,yy,zz);
  if(doijk==3) sprintf(gp.title,"projection of spindensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of spindensity div Ms(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of spindensity |Ms(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'o': dim=ORBMOMDENS_EV_DIM;
printf("#orbital momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm  aL(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates %s, orbital momdensity L(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),coord,xx,yy,zz);
  if(doijk==3) sprintf(gp.title,"projection of orbmomdensity Ms(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of orbmomdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of orbmomdensity |ML(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'm': dim=SPINDENS_EV_DIM+ORBMOMDENS_EV_DIM;
printf("#magnetic momdensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm (aS(l,m) R^2(r)+ aL(l,m) F(r)) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates %s, magnetic momdensity M(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),coord,xx,yy,zz);
  if(doijk==3) sprintf(gp.title,"projection of momdensity M(r).(%g,%g,%g)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of momdensity div ML(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of momdensity |ML(r)|");
if(doijk<3){dim*=3;}
gp.threshhold=0.05;
break;
 case 'j': dim=ORBMOMDENS_EV_DIM;
printf("#currdensity is expanded in tesseral harmonics Zlm\n\
#   j(r).(%g,%g,%g)= sum_lm (b(l,m) R^2(r)+ d(l,m) F(r) Zlm(Omega)\n\
#   with F(r)=1/r int_r^inf R^2(x) dx\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT with coordinates %s, currentdensity j(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),coord,xx,yy,zz);
  if(doijk==3) sprintf(gp.title,"projection of currdensity j(r).(i=%g,j=%g,k=%g)(milliAmp/A^2)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of currdensity div j(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of currdensity |j(r)|(milliAmp/A^2)");
  dim*=6;
gp.threshhold=0.05;
break;
 default: help_and_exit();break;
}

if(gp.read())printf("#reading graphic parameters from results/graphic_parameters.set\n");
// read graphic parameters which are set by user in file results/graphic_parameters.set
//  if(jjjps.module_type!=0&&jjjps.module_type!=2&&jjjps.module_type!=4)
//  {fprintf(stderr,"ERROR densplt: calculation not possible for this single ion module\n");exit(EXIT_FAILURE);}
//  if (jjjps.module_type==0&&jjjps.gJ!=0)
//  {fprintf(stderr,"************** WARNING **********************\n reading external single ion module with gJ not zero: gJ=%g - please check if calculation of density is supported !!!\n*********************************************\n",jjjps.gJ);}

  Vector moments(1,dim),mom(1,arrowdim);mom=0;

  Vector momS(1,SPINDENS_EV_DIM);
  Vector momL(1,ORBMOMDENS_EV_DIM);
  Vector momentsx(1,SPINDENS_EV_DIM);
  Vector momentsy(1,SPINDENS_EV_DIM);
  Vector momentsz(1,SPINDENS_EV_DIM);
  Vector momentlx(1,ORBMOMDENS_EV_DIM);
  Vector momently(1,ORBMOMDENS_EV_DIM);
  Vector momentlz(1,ORBMOMDENS_EV_DIM);
  Vector h(1,jjjps.nofcomponents);h=0; // exchange field =0 ... dimension 6 ok for every module ?

  jjjps.Icalc_parameter_storage_init(h,Hext,T);

  printf("calculating expectation values of density coefficients ....\n");

switch(arrow)
{case 1: jjjps.Scalc(mom,T,h,Hext,jjjps.Icalc_parstorage);break;
 case 2: jjjps.Lcalc(mom,T,h,Hext,jjjps.Icalc_parstorage);break;
 case 3: jjjps.mcalc(mom,T,h,Hext,jjjps.Icalc_parstorage);break;
}

switch(argv[1][0]) // dimension definition from jjjpar.hpp
{case 'c':  jjjps.chargedensity_coeff (moments, T, h, Hext, jjjps.Icalc_parstorage); break;
 case 's':  if(xx!=0||doijk<3)jjjps.spindensity_coeff (momentsx,1, T, h,Hext, jjjps.Icalc_parstorage);
            if(yy!=0||doijk<3)jjjps.spindensity_coeff (momentsy,2, T, h,Hext, jjjps.Icalc_parstorage);
            if(zz!=0||doijk<3)jjjps.spindensity_coeff (momentsz,3, T, h,Hext, jjjps.Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=SPINDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'o':  if(xx!=0||doijk<3)jjjps.orbmomdensity_coeff (momentsx,1, T, h,Hext, jjjps.Icalc_parstorage);
            if(yy!=0||doijk<3)jjjps.orbmomdensity_coeff (momentsy,2, T, h,Hext, jjjps.Icalc_parstorage);
            if(zz!=0||doijk<3)jjjps.orbmomdensity_coeff (momentsz,3, T, h,Hext, jjjps.Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+ORBMOMDENS_EV_DIM)=momentsy(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'm':  if(xx!=0||doijk<3)jjjps.spindensity_coeff (momentsx,1, T, h,Hext, jjjps.Icalc_parstorage);
            if(yy!=0||doijk<3)jjjps.spindensity_coeff (momentsy,2, T, h,Hext, jjjps.Icalc_parstorage);
            if(zz!=0||doijk<3)jjjps.spindensity_coeff (momentsz,3, T, h,Hext, jjjps.Icalc_parstorage);
            momS=xx*momentsx+yy*momentsy+zz*momentsz;
            if(xx!=0||doijk<3)jjjps.orbmomdensity_coeff (momentlx,1, T, h,Hext, jjjps.Icalc_parstorage);
            if(yy!=0||doijk<3)jjjps.orbmomdensity_coeff (momently,2, T, h,Hext, jjjps.Icalc_parstorage);
            if(zz!=0||doijk<3)jjjps.orbmomdensity_coeff (momentlz,3, T, h,Hext, jjjps.Icalc_parstorage);
            momL=xx*momentlx+yy*momently+zz*momentlz;
            for(int i=1;i<=SPINDENS_EV_DIM;++i){
            if(doijk==3){moments(i)=momS(i);moments(i+SPINDENS_EV_DIM)=momL(i);}
            else{moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);
                 moments(i+3*SPINDENS_EV_DIM)=momentlx(i);moments(i+4*SPINDENS_EV_DIM)=momently(i);moments(i+5*SPINDENS_EV_DIM)=momentlz(i);
                }
                                               }
            break;
 case 'j':  jjjps.orbmomdensity_coeff (momentlx,1, T, h,Hext, jjjps.Icalc_parstorage);
            jjjps.orbmomdensity_coeff (momently,2, T, h,Hext, jjjps.Icalc_parstorage);
            jjjps.orbmomdensity_coeff (momentlz,3, T, h,Hext, jjjps.Icalc_parstorage);
            for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){
             moments(i)=momentlx(i);moments(i+ORBMOMDENS_EV_DIM)=momently(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentlz(i);
             }
            break;
 default: help_and_exit();
}

int i,nofpc=0;


 cs.sipffilenames[1]=new char[MAXNOFCHARINLINE];
 cs.abc(1)=6.0;cs.abc(2)=6.0;cs.abc(3)=6.0;cs.nofatoms=nofpc+1;
 cs.r=0;cs.r(1,1)=1.0;cs.r(2,2)=1.0;cs.r(3,3)=1.0;
 strcpy(cs.sipffilenames[1],argv[2+os]);
 cs.x[1]=0.5*gp.scale_view_1;cs.y[1]=0.5*gp.scale_view_2;cs.z[1]=0.5*gp.scale_view_3; // put atom in middle of cell

// read pointcharge-parameters
if(gp.show_pointcharges>0) nofpc=read_pointcharge_parameters(gp,cs.sipffilenames,argv[2+os],cs.x,cs.y,cs.z,jjjps,cs.abc);
  spincf s(1,1,1,nofpc+1,dim);
  spincf magmom(1,1,1,nofpc+1,3);
  Vector hkl(1,3);hkl=0;s=s*0;magmom=magmom*0;
  spincf ev_real(s),ev_imag(s);
  spincf magmomev_real(magmom),magmomev_imag(magmom);

  for(i=1;i<=dim;++i)s.m(1,1,1)(i)=moments(i);
  for(i=1;i<=3;++i)magmom.m(1,1,1)(i)=mom(i);

  printf("# T=%g K field H=(%g,%g,%g) Tesla\n",T,Hext(1),Hext(2),Hext(3));
  printf("#%s\n",gp.title);

  fout = fopen_errchk ("results/densplt.jvx", "w");
   s.jvx_cd(fout,text,cs,gp,0.0,ev_real,ev_imag,hkl,T,h,Hext,magmom,magmomev_real,magmomev_imag);
  fclose (fout);

  fout = fopen_errchk ("results/densplt.grid", "w");
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/densplti.grid", "w");
  gp.gridi=1;gp.gridj=200;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/denspltj.grid", "w");
  gp.gridi=200;gp.gridj=1;gp.gridk=200;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
  fout = fopen_errchk ("results/denspltk.grid", "w");
  gp.gridi=200;gp.gridj=200;gp.gridk=1;
  s.cd(fout,cs,gp,ev_real,ev_imag,0.0,hkl,T,h,Hext);
  fclose (fout);
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program densplt\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/densplt.jvx\n");
fprintf(stderr,"# * displaycontour 2 3 4 results/densplti.grid\n");
fprintf(stderr,"# * displaycontour 1 3 4 results/denspltj.grid\n");
fprintf(stderr,"# * displaycontour 1 2 4 results/denspltk.grid\n");
fprintf(stderr,"# * saved density mesh in results/densplt.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  int dj;
  for(dj=1;dj<=nofpc+1;++dj){delete cs.sipffilenames[dj];}
  return 0;
}


