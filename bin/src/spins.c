/****************************************************
 * spins - display spinconfiguration at given htpoint
 * Author: Martin Rotter
 ****************************************************/
#include "../../version"
#include <par.hpp>
#include "spincf.hpp"
#include "martin.h"
#include "graphic_parameters.hpp"
#include "cryststruct.hpp"
#include "densities_func.c"

void help_and_exit()
    { printf ("\n\
program spins - popout spin/exchange field configuration\n\
              - and/or display 3d animation of spin/moment/densities and animations\n\n\
use as: spins -f mcphas.sps T Ha Hb Hc\n\
    or: spins [-c|-s|-o|-m|-j] [-p i j k|-div] [-S|-L|-M] [-P] T Ha Hb Hc [h k l E]\n\
                    \n\
1) if used with -f filename this file has to be a mcphas.mf or mcphas.sps file, the spin configuration\n\
   at given temperature T[K] and magnetic effective field H[T]\n \
   is read and extracted from this file and printed on screen (stdout), nothing else is done\n\
  \n\
2) if used without a filename, the information is read from results/mcphas.* results/mcdisp.*\n\
   output files and 3d graphical animations are created.\n\
   options are:\n\
         -c ... calculate chargedensity\n\
         -s ... calculate spindensity\n\
         -o ... calculate angular orbital momentum density\n\
         -m ... calculate magnetic moment density\n\
         -j ... calculate currentdensity\n\
         -p i j k ... calculate projection of spin/orbital/current/magnetic moment density\n\
                  along direction i j k, e.g. 0 0 1\n\
         -div    ... calculate divergence of spin/orbital/current/magnetic moment density  \n\
         -S  ... show arrow indicating spin\n\
         -L  ... show arrow indicating orbital angular momentum\n\
         -M  ... show arrow indicating magnetic moment (for cluster show total moment)\n\
         -Mi ... show arrow indicating magnetic moment (for cluster show individual moments)\n\
         -P  ... calculate phononic displacement\n\
         note, that in order to animate changes in the above quantities, the corresponding\n\
         switch has to be enabled in the mcdisp calculation (mcdisp.par) and the single ion\n\
         modules have to be capable of calculating the corresponding observables. \n\
 \n\
     example:\n\
        spins -c 2 0 0 1\n\
        ...calculates the charge density at T=2K and H=(0,0,1) Tesla\n \
\n\
 This program outputs a magnetic structure (and magnetic excitation)\n \
 graphic/movie in the output files of different format:\n \
 results/spins*.eps (postscript), results/spins*.fst (fp_studio), \n \
 results/spins.out (ascii) and results/spins*.jvx (javaview)\n\n \
 the graphics output format can be fine tuned in .sps and .qev input files\n"
" or results/graphics_parameters.set by show_abc_unitcell,\n"
" show_primitive_crystal_unitcell, spins_scale_moment, spins_wave_amplitude\n \
 show_magnetic_unitcell, show_atoms, scale_view_1,scale_view_2, scale_view_3 ...\n\n \
 jvx files can be viewed by:\n\
 java javaview results/spins.jvx \n \
 java javaview \"model=results/spins.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n"
" gif images stored by javaview can be connected to movies by ImageMagick \n"
" convert  -delay 1 -size 100x100 -loop 1 geomAnim.*.gif output.gif\n");
 exit (1);
    }


/**********************************************************************/
// hauptprogramm
int main (int argc, char **argv)
{ 
printf("# **********************************************************\n");
printf("# * spins - display 3d graphics of spins,moments,densities,*\n");
printf("# * at given H and T                                       *\n");
printf("# * Reference: M. Rotter PRB 79 (2009) 140405R             *\n");
printf("# * %s                                     *\n",MCPHASVERSION);
printf("# **********************************************************\n");

 FILE * fin, * fout;
double T; Vector Hext(1,3);
 int i,n=0;//,dophon=0;
 cryststruct cs,cs4;
 spincf savmf;
 //float numbers[13];numbers[9]=1;numbers[10]=3;
 //numbers[0]=13;
 char outstr[MAXNOFCHARINLINE];
  int dim=28;
 char text[1000];
 int os=0; int doijk=0,arrow=0,density=0;//,arrowdim=3;
 double xx=0,yy=0,zz=0;
graphic_parameters gp;
gp.show_abc_unitcell=1.0;
gp.show_primitive_crystal_unitcell=1.0;
gp.show_magnetic_unitcell=1.0;
gp.show_atoms=1.0;
gp.scale_view_1=1.0;
gp.scale_view_2=1.0;
gp.scale_view_3=1.0;
gp.spins_scale_moment=0;
gp.show_density=0;
sprintf(gp.title,"output of program spins");

 // check command line
 if (argc < 5){help_and_exit();}
// first: option without graphics just screendump <I> or exchange field configuration at given HT
 if (strcmp(argv[1],"-f")==0)
 { fin = fopen_errchk (argv[2], "rb");os=2;printf("# reading from file %s\n",argv[2]);}
 else  // second ... other options with graphics !!
 { fin = fopen_errchk ("./results/mcphas.mf", "rb");printf("# reading from file ./results/mcphas.mf\n");
  if(strcmp(argv[1],"-c")==0){os=1;}
  if(strcmp(argv[1],"-s")==0){os=1;}
  if(strcmp(argv[1],"-o")==0){os=1;} 
  if(strcmp(argv[1],"-m")==0){os=1;}
  if(strcmp(argv[1],"-j")==0){os=1;}
  if(os==1)
  {gp.show_density=1;density=1;

switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c': dim=CHARGEDENS_EV_DIM;
printf("#chargedensity is expanded in tesseral harmonics Zlm\n\
#   ro(r) sum_lm (a(l,m) R^2(r) Zlm(Omega)\n\
#   M. Rotter et al. J Phys: Conf Ser. 325 (2011) 012005\n#\n ");
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , chargedensity ro(r)</title>\n", T,Hext(1),Hext(2),Hext(3));
 sprintf(gp.title,"chargedensity ro(r)");
 gp.threshhold=-0.05;
           break;
 case 's': dim=SPINDENS_EV_DIM;
printf("#spindensity is expanded in tesseral harmonics Zlm\n\
#   M(r).(%g,%g,%g)= sum_lm aS(l,m) R^2(r) Zlm(Omega)\n\
#   E. Balcar J. Phys. C. 8 (1975) 1581\n#\n ",xx,yy,zz);
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , spindensity S(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
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
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , orbital momdensity L(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
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
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , magnetic momdensity M(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
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
 sprintf(text,"<title>T=%4gK h||a=%4gT h||b=%4gT h||c=%4gT , currentdensity j(r).(%g,%g,%g)</title>\n", T,Hext(1),Hext(2),Hext(3),xx,yy,zz);
  if(doijk==3) sprintf(gp.title,"projection of currdensity j(r).(i=%g,j=%g,k=%g)(milliAmp/A^2)",xx,yy,zz);
  if(doijk==1){sprintf(gp.title,"divergence of currdensity div j(r)");gp.scale_density_vectors=0;}
  if(doijk==0) sprintf(gp.title,"abs value  of currdensity |j(r)|(milliAmp/A^2)");
  dim*=6;
gp.threshhold=0.05;
break;
 default: help_and_exit();break;
}
  }
  if(strcmp(argv[os+1],"-div")==0){os+=1;doijk=1;}
  else if(strcmp(argv[os+1],"-p")==0){os+=4;
  xx=strtod(argv[3],NULL);
  yy=strtod(argv[4],NULL);
  zz=strtod(argv[5],NULL);
  double rr;
  // normalize direction vector
  rr=sqrt(xx*xx+yy*yy+zz*zz);
  xx/=rr;yy/=rr;zz/=rr;
  doijk=3;
                                     }

if(strcmp(argv[1+os],"-S")==0){os+=1;arrow=1;gp.spins_colour=3; gp.spins_scale_moment=1;//arrowdim=SPIN_EV_DIM;
                              sprintf(gp.title,"%s arrows correspond to the spins",gp.title);}
else if(strcmp(argv[1+os],"-L")==0){os+=1;arrow=2;gp.spins_colour=2; gp.spins_scale_moment=1;//arrowdim=ORBMOM_EV_DIM;
                                   sprintf(gp.title,"%s arrows correspond to the orbital angular momenta",gp.title);}
else if(strncmp(argv[1+os],"-M",2)==0){os+=1;arrow=3;gp.spins_colour=1; gp.spins_scale_moment=1;//arrowdim=MAGMOM_EV_DIM;
                                   sprintf(gp.title,"%s arrows correspond to the magnetic moments",gp.title);
                                   if(strcmp(argv[os],"-Mi")==0){arrow=4;}
                                   }

if(strcmp(argv[1+os],"-P")==0){os+=1;}//dophon=1;}

 }

   fout = fopen_errchk ("./results/spins.out", "w");
// input file header and conf------------------------------------------------------------------
   n=headerinput(fin,fout,gp,cs);cs4.abc=cs.abc;cs4.r=cs.r;cs4.nofatoms=cs.nofatoms;cs4.nofcomponents=cs.nofcomponents;
// load spinsconfigurations and check which one is nearest -------------------------------   
check_for_best(fin,strtod(argv[1+os],NULL),strtod(argv[2+os],NULL),strtod(argv[3+os],NULL),strtod(argv[4+os],NULL),savmf,T,Hext,outstr);
fclose (fin);

  printf("#! %s - configuration\n",outstr);
  savmf.print(stdout);
  if (strcmp(argv[1],"-f")==0) {fclose(fout);exit(0);}

// FROM HERE ON IT IS ONLY EXECUTED IF GRAPHICS ARE DESIRED ... 

gp.read();

  int ii,nt,k,j;
  par inputpars("./mcphas.j");
  Vector hh(1,savmf.nofcomponents*savmf.nofatoms);
  spincf densitycf(savmf.na(),savmf.nb(),savmf.nc(),savmf.nofatoms,dim);
  ii=0; 
  // if individual ions of the cluster are to be shown make nofatoms in spincf larger !
  for (j=1;j<=inputpars.nofatoms;++j){
         if(arrow==4&&(*inputpars.jjj[j]).module_type==5){for(k=1;k<=(*(*inputpars.jjj[j]).clusterpars).nofatoms;++k)
                                                 {++ii;par inputpars4((*(*inputpars.jjj[j]).clusterpars));
                                                 cs4.x[ii]=cs.x[j]+(*inputpars4.jjj[k]).xyz(1);
                                                 cs4.y[ii]=cs.y[j]+(*inputpars4.jjj[k]).xyz(2);
                                                 cs4.z[ii]=cs.z[j]+(*inputpars4.jjj[k]).xyz(3);cs4.sipffilenames[ii]=new char[MAXNOFCHARINLINE];
                                                 strcpy(cs4.sipffilenames[ii],(*inputpars4.jjj[k]).sipffilename);
                                                 //check if cluster abc  are the same as inputpars
           if(fabs(cs4.abc(1)-inputpars4.a)>1e-6)
              {fprintf(stderr,"Error program spins - a=%g in cluster %s not the same as a=%g in mcphas.j\n",cs4.abc(1),(*inputpars.jjj[j]).sipffilename,inputpars4.a);exit(1);}
           if(fabs(cs4.abc(2)-inputpars4.b)>1e-6)
              {fprintf(stderr,"Error program spins - b=%g in cluster %s not the same as b=%g in mcphas.j\n",cs4.abc(2),(*inputpars.jjj[j]).sipffilename,inputpars4.b);exit(1);}
           if(fabs(cs4.abc(3)-inputpars4.c)>1e-6)
              {fprintf(stderr,"Error program spins - c=%g in cluster %s not the same as c=%g in mcphas.j\n",cs4.abc(3),(*inputpars.jjj[j]).sipffilename,inputpars4.c);exit(1);}
           if(fabs(cs4.abc(4)-inputpars4.alpha)>1e-6)
              {fprintf(stderr,"Error program spins - alpha=%g in cluster %s not the same as alpha=%g in mcphas.j\n",cs4.abc(4),(*inputpars.jjj[j]).sipffilename,inputpars4.alpha);exit(1);}
           if(fabs(cs4.abc(5)-inputpars4.beta)>1e-6)
              {fprintf(stderr,"Error program spins - beta=%g in cluster %s not the same as beta=%g in mcphas.j\n",cs4.abc(5),(*inputpars.jjj[j]).sipffilename,inputpars4.beta);exit(1);}
           if(fabs(cs4.abc(6)-inputpars4.gamma)>1e-6)
              {fprintf(stderr,"Error program spins - gamma=%g in cluster %s not the same as gamma=%g in mcphas.j\n",cs4.abc(6),(*inputpars.jjj[j]).sipffilename,inputpars4.gamma);exit(1);}                                                 
                                                 }
                                               }
                                               else
                                              {++ii;cs4.x[ii]=cs.x[j];cs4.y[ii]=cs.y[j];cs4.z[ii]=cs.z[j];cs4.sipffilenames[ii]=cs.sipffilenames[j];}
              } cs4.nofatoms=ii;

  spincf spinconf(savmf.na(),savmf.nb(),savmf.nc(),ii,3);


// the following is for the printout of spins.out ...........................
fprintf(fout,"#!T=%g K Ha=%g T Hb= %g T Hc= %g T: nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",T,Hext(1),Hext(2),Hext(3),savmf.na(),savmf.nb(),savmf.nc(),cs4.nofatoms*savmf.na()*savmf.nb()*savmf.nc());
fprintf(fout,"#{sipf-file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3] <Ma> <Mb> <Mc> [mb] [optional <Sa> <La> <Sb> <Lb> <Sc> <Lc>\n");
fprintf(fout,"#          corresponding exchange fields hxc [meV]- if passed to mcdiff only these are used for calculation (not the magnetic moments)\n");
// determine primitive magnetic unit cell
Matrix p(1,3,1,3);Vector xyz(1,3),dd0(1,3);
savmf.calc_prim_mag_unitcell(p,cs.abc,cs.r);
// .............................................................................                                
	       
//  1. from the meanfieldconfiguration (savmf) the <Olm> have to be calculated for all l=2,4,6
// 1.a: the mcphas.j has to be used to determine the structure + single ione properties (copy something from singleion.c)
// 1.b: Icalc has to be used to calculate all the <Olm>.
Vector h(1,inputpars.nofcomponents);
h=0;for(ii=1;ii<=inputpars.nofatoms;++ii)
{(*inputpars.jjj[ii]).Icalc_parameter_storage_init(h,Hext,T);} // initialize Icalc module parameter storage

 for (i=1;i<=savmf.na();++i){for(j=1;j<=savmf.nb();++j){for(k=1;k<=savmf.nc();++k)
 {
    hh=savmf.m(i,j,k);
  densitycf.m(i,j,k)=0;int i4=1,ii4=1;
  for(ii=1;ii<=inputpars.nofatoms;++ii)
 {  
    Vector magmom(1,3),mom(1,3);
    Vector Lmom(1,3);
    Vector Smom(1,3);
    Vector moments(1,dim);
      Vector momS(1,SPINDENS_EV_DIM);
  Vector momL(1,ORBMOMDENS_EV_DIM);
  Vector momentsx(1,SPINDENS_EV_DIM);
  Vector momentsy(1,SPINDENS_EV_DIM);
  Vector momentsz(1,SPINDENS_EV_DIM);
  Vector momentlx(1,ORBMOMDENS_EV_DIM);
  Vector momently(1,ORBMOMDENS_EV_DIM);
  Vector momentlz(1,ORBMOMDENS_EV_DIM);
    h=0;
   for(nt=1;nt<=inputpars.nofcomponents;++nt){h(nt)=hh(nt+inputpars.nofcomponents*(ii-1));}

switch(arrow)
{case 1: (*inputpars.jjj[ii]).Scalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 2: (*inputpars.jjj[ii]).Lcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
 case 3: (*inputpars.jjj[ii]).mcalc(mom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);break;
}
switch(arrow)
{case 1:
 case 2:
 case 3:
         for(nt=1;nt<=3;++nt)
		        {spinconf.m(i,j,k)(nt+3*(ii-1))=mom(nt);
                    };break;
 case 4: int dim4;
         dim4=3;if((*inputpars.jjj[ii]).module_type==5)
                           dim4=(*(*inputpars.jjj[ii]).clusterpars).nofatoms*3;
         Vector momi(1,dim4);
         (*inputpars.jjj[ii]).micalc(momi,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage);
         // now put the components of momi to the spinconf                                               
         for(nt=1;nt<=dim4;++nt){spinconf.m(i,j,k)(nt+i4-1)=momi(nt);};i4+=dim4;
         break;
}

// output atoms and moments in primitive unit cell to fout  ------------------------------------
Vector dd3(1,3);
if(arrow==4&&(*inputpars.jjj[ii]).module_type==5){
    for(nt=1;nt<=(*(*inputpars.jjj[ii]).clusterpars).nofatoms;++nt)
     {dd3=spinconf.pos(i,j,k,ii4, cs4);
      dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
      fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            cs4.sipffilenames[ii4],dd3(1)/cs.abc(1),dd3(2)/cs.abc(2),dd3(3)/cs.abc(3),dd0(1),dd0(2),dd0(3));
            fprintf(fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(1+3*(ii4-1))));
            fprintf(fout," %4.4f",myround(1e-5,spinconf.m(i,j,k)(2+3*(ii4-1))));
            fprintf(fout," %4.4f\n",myround(1e-5,spinconf.m(i,j,k)(3+3*(ii4-1))));
       // do not print out L S or exchange fields for cluster module
      ++ii4;
     }
    }
else{++ii4;
    // output the positions
    dd3=savmf.pos(i,j,k,ii, cs);
    dd0=p.Inverse()*dd3;dd0(1)*=savmf.na();dd0(2)*=savmf.nb();dd0(3)*=savmf.nc();
    fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
            cs.sipffilenames[ii],dd3(1)/cs.abc(1),dd3(2)/cs.abc(2),dd3(3)/cs.abc(3),dd0(1),dd0(2),dd0(3));
    //ouput the magnetic moment if possible
    if((*inputpars.jjj[ii]).mcalc(magmom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
   {           for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f",myround(1e-5,magmom(nt)));}
     // and output the orbital and spin momentum if possible 
     if((*inputpars.jjj[ii]).Lcalc(Lmom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage)&&
        (*inputpars.jjj[ii]).Scalc(Smom,T,h,Hext,(*inputpars.jjj[ii]).Icalc_parstorage))
      {        for(nt=1;nt<=3;++nt){fprintf(fout," %4.4f %4.4f",myround(1e-5,Smom(nt)),myround(1e-5,Lmom(nt)));}
      }
   }
    // finally output a line with the exchange fields 
  fprintf(fout,"\n                 corresponding exchange fields hxc [meV]-->          ");
                      for(nt=1;nt<=savmf.nofcomponents;++nt)  // printout exchangefields
                        {fprintf(fout," %4.4f",myround(1e-5,h(nt)));}
                         fprintf(fout,"\n");
  }
// -----------------------------------------------------------------------------------------------


if(density){
switch(argv[1][1]) // dimension definition from jjjpar.hpp
{case 'c':  (*inputpars.jjj[ii]).chargedensity_coeff (moments, T, h, Hext, (*inputpars.jjj[ii]).Icalc_parstorage); break;
 case 's':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsy,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=SPINDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'o':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsy,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentsz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(doijk==3){ moments=xx*momentsx+yy*momentsy+zz*momentsz;}
            else{for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){moments(i)=momentsx(i);moments(i+ORBMOMDENS_EV_DIM)=momentsy(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentsz(i);}
                }
            break;
 case 'm':  if(xx!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsy,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).spindensity_coeff (momentsz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            momS=xx*momentsx+yy*momentsy+zz*momentsz;
            if(xx!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentlx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(yy!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momently,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            if(zz!=0||doijk<3)(*inputpars.jjj[ii]).orbmomdensity_coeff (momentlz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            momL=xx*momentlx+yy*momently+zz*momentlz;
            for(int i=1;i<=SPINDENS_EV_DIM;++i){
            if(doijk==3){moments(i)=momS(i);moments(i+SPINDENS_EV_DIM)=momL(i);}
            else{moments(i)=momentsx(i);moments(i+SPINDENS_EV_DIM)=momentsy(i);moments(i+2*SPINDENS_EV_DIM)=momentsz(i);
                 moments(i+3*SPINDENS_EV_DIM)=momentlx(i);moments(i+4*SPINDENS_EV_DIM)=momently(i);moments(i+5*SPINDENS_EV_DIM)=momentlz(i);
                }
                                               }
            break;
 case 'j':  (*inputpars.jjj[ii]).orbmomdensity_coeff (momentlx,1, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            (*inputpars.jjj[ii]).orbmomdensity_coeff (momently,2, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            (*inputpars.jjj[ii]).orbmomdensity_coeff (momentlz,3, T, h,Hext, (*inputpars.jjj[ii]).Icalc_parstorage);
            for(int i=1;i<=ORBMOMDENS_EV_DIM;++i){
             moments(i)=momentlx(i);moments(i+ORBMOMDENS_EV_DIM)=momently(i);moments(i+2*ORBMOMDENS_EV_DIM)=momentlz(i);
             }
            break;
 default: help_and_exit();
}
                   for(nt=1;nt<=dim;++nt)
		        {densitycf.m(i,j,k)(nt+dim*(ii-1))=moments(nt);
                    }
} // gp.show_density

  }}
}}
  fclose (fout);
  
// create plot of spinconfiguration -----------------------------------------------------------
printf("# ************************************************************************\n");
printf("#%s\n",gp.title);
printf("# ************************************************************************\n");
              if(arrow==0)gp.spins_scale_moment=0;
              if(density==0)gp.show_density=0;

    fin = fopen_errchk ("./results/spins.eps", "w");
     spinconf.eps(fin,outstr);
    fclose (fin);

// here the 3d file should be created
    fin = fopen_errchk ("./results/spinsab.eps", "w");
Vector gJJ(1,spinconf.nofatoms); for (i=1;i<=spinconf.nofatoms;++i){gJJ(i)=1;}
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,1,gJJ,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spinsac.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,2,gJJ,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spinsbc.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,3,gJJ,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spins3dab.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,4,gJJ,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spins3dac.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,5,gJJ,spinconf);
    fclose (fin);
    fin = fopen_errchk ("./results/spins3dbc.eps", "w");
     spinconf.eps3d(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,6,gJJ,spinconf);
    fclose (fin);

    fin = fopen_errchk ("./results/spins.fst", "w");
     spinconf.fst(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,gJJ,spinconf);
    fclose (fin);

    
   fin = fopen_errchk ("./results/spins_prim.fst", "w");
     spinconf.fstprim(fin,outstr,cs4.abc,cs4.r,cs4.x,cs4.y,cs4.z,gJJ,spinconf);
    fclose (fin);
             Vector hkl(1,3);hkl=0;
             Vector gjmbHxc(1,3);gjmbHxc=0;
             spincf densityev_real(densitycf*0.0);
             spincf densityev_imag(densitycf*0.0);
             spincf spinconfev_real(spinconf*0.0);
             spincf spinconfev_imag(spinconf*0.0);
            // to do jvx output of static structure put zeros into these spinconfigurations

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin = fopen_errchk ("./results/spins.jvx", "w");
    gp.showprim=0;gp.spins_wave_amplitude=0;
     densitycf.jvx_cd(fin,outstr,cs,gp,0.0,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag);
    fclose (fin);

// create jvx file of spinconfiguration - checkout polytope/goldfarb3.jvx  primitive/cubewithedges.jvx
   fin = fopen_errchk ("./results/spins_prim.jvx", "w");
     gp.showprim=1;
     densitycf.jvx_cd(fin,outstr,cs,gp,0.0,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag);
    fclose (fin);

//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
// try a spinwave picture !!!  ... include phonons and spindensity changes ...
//***************************************************************************************************************
//***************************************************************************************************************

if (argc-os>=6){
              // double E;
             long int pos=0;
             int extended_eigenvector_dimension;
              char instr[MAXNOFCHARINLINE];
              float numbers[13];numbers[9]=1;numbers[10]=3;
              numbers[0]=13;
             gp.spins_wave_amplitude=1.0;gp.spins_show_ellipses=1.0;gp.spins_show_oscillation=1.0;
//----------------------------------------------------------------------------------------------------------
            if(arrow>0){double checkdd=1e7;
             switch(arrow)
             {case 1:  fin = fopen_errchk ("./results/mcdisp.qes", "rb");break;
              case 2:  fin = fopen_errchk ("./results/mcdisp.qeo", "rb");break;
              case 3:  fin = fopen_errchk ("./results/mcdisp.qem", "rb");break;
              case 4:  fprintf(stderr,"mcdisp: output of individual moment oscillation in eigenvector file mcdisp.qemi not yet implemented - thus exiting program spins\n");
                       fin = fopen_errchk ("./results/mcdisp.qemi", "rb");break;
             }
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_oscillation",gp.spins_show_oscillation);
              }

               j=fseek(fin,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               double delta,dd,ddT,ddHa,ddHb,ddHc,ddh,ddk,ddl,ddE;
               for (delta=1000.0;feof(fin)==0                      //end of file
                    &&(n=inputline(fin,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin); 
                 spincf ev_real(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,3);
                 spincf ev_imag(spinconf.na(),spinconf.nb(),spinconf.nc(),spinconf.nofatoms,3);
                 ev_real.load(fin);ev_imag.load(fin);
                 ddT=T-numbers[4];ddT*=ddT;
                 ddHa=Hext(1)-numbers[1];ddHa*=ddHa;
                 ddHb=Hext(2)-numbers[2];ddHb*=ddHb;
                 ddHc=Hext(3)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5+os],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6+os],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7+os],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8+os],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;checkdd=fabs(T-numbers[4])+fabs(Hext(1)-numbers[1])+fabs(Hext(2)-numbers[2])+fabs(Hext(3)-numbers[3]);
                  sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];//E=numbers[9]; 
                  spinconfev_real=ev_real;
                  spinconfev_imag=ev_imag;                  
                 }
                 pos=ftell(fin); 
                 fgets(instr,MAXNOFCHARINLINE,fin); 
                 while (instr[strspn(instr," \t")]=='#'&&feof(fin)==0) // pointer to 'ltrimstring' 
                  {pos=ftell(fin);fgets(instr,MAXNOFCHARINLINE,fin);}
                 j=fseek(fin,pos,SEEK_SET);
               }
              fclose (fin);
             if(checkdd>1e-7)
               {fprintf(stderr,"Error program spins - inconsistent output files mcphas.mf and mcdisp.*:  temperature/magnetic field in static configuration\n" 
                               "(from results/mcphas.mf) is different from mcdisp output files results/mcdisp.qes,qem,qeo:\n %s\n"
                               "... probably you need to rerun setup_mcdisp_mf and mcdisp with T=%g Ha=%g Hb=%g Hc=%g",outstr,T,Hext(1),Hext(2),Hext(3));exit(EXIT_FAILURE);}
               fprintf(stdout,"#%s - moment oscillation - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              spinconfev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              spinconfev_imag.print(stdout);
             }//arrow
     
//----------------------------------------------------------------------------------------------------------
            if(density){double checkdd=1e7;
             switch(argv[1][1]) // dimension definition from jjjpar.hpp
                {case 'c':fin = fopen_errchk ("./results/mcdisp.qee", "rb");break;
                 case 's':fin = fopen_errchk ("./results/mcdisp.qsd", "rb");break;
                 case 'o':fin = fopen_errchk ("./results/mcdisp.qod", "rb");break;
                 case 'm':fprintf(stderr,"Error spins: magnetic moment density oscillation not yet implemented\n");exit(1);break;
                             // would have to look into qsd and qod files !!
                 case 'j':fin = fopen_errchk ("./results/mcdisp.qod", "rb");break;
                }
             // input file header ------------------------------------------------------------------
             instr[0]='#';
              while (instr[strspn(instr," \t")]=='#') // pointer to 'ltrimstring' 
              { pos=ftell(fin); 
                if (pos==-1) 
                {fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
                fgets(instr,MAXNOFCHARINLINE,fin); 
                // inserted 4.4.08 in order to format output correctly (characterstring 13 spoiled output string)
                for(i=0;(unsigned int)i<=strlen(instr);++i){if(instr[i]==13)instr[i]=32;}
               // load evs and check which one is nearest -------------------------------   
               extract(instr,"spins_wave_amplitude",gp.spins_wave_amplitude);
               extract(instr,"spins_show_ellipses",gp.spins_show_ellipses);
               extract(instr,"spins_show_oscillation",gp.spins_show_oscillation);
               extract(instr,"extended_eigenvector_dimension",extended_eigenvector_dimension);
              }

               j=fseek(fin,pos,SEEK_SET); 
               if (j!=0){fprintf(stderr,"Error: wrong qev file format\n");exit (EXIT_FAILURE);}
   
               double delta,dd,ddT,ddHa,ddHb,ddHc,ddh,ddk,ddl,ddE;
               for (delta=1000.0;feof(fin)==0                      //end of file
                    &&(n=inputline(fin,numbers))>=8   //error in line reading (8 old format, 9 new format)
		    ;)

               { fgets(instr,MAXNOFCHARINLINE,fin); 
                 spincf ev_real(densitycf.na(),densitycf.nb(),densitycf.nc(),densitycf.nofatoms,extended_eigenvector_dimension);
                 spincf ev_imag(densitycf.na(),densitycf.nb(),densitycf.nc(),densitycf.nofatoms,extended_eigenvector_dimension);
                 ev_real.load(fin);ev_imag.load(fin);
                 ddT=T-numbers[4];ddT*=ddT;
                 ddHa=Hext(1)-numbers[1];ddHa*=ddHa;
                 ddHb=Hext(2)-numbers[2];ddHb*=ddHb;
                 ddHc=Hext(3)-numbers[3];ddHc*=ddHc;
                 ddh=strtod(argv[5+os],NULL)-numbers[5];ddh*=ddh;
                 ddk=strtod(argv[6+os],NULL)-numbers[6];ddk*=ddk;
                 ddl=strtod(argv[7+os],NULL)-numbers[7];ddl*=ddl;
                 ddE=strtod(argv[8+os],NULL)-numbers[9];ddE*=ddE;
                 
                 dd=sqrt(ddT+ddHa+ddHb+ddHc+ddh+ddk+ddl+ddE+0.000001);
                 if (dd<delta)
                 {delta=dd;checkdd=fabs(T-numbers[4])+fabs(Hext(1)-numbers[1])+fabs(Hext(2)-numbers[2])+fabs(Hext(3)-numbers[3]);
                  sprintf(outstr,"T=%g Ha=%g Hb=%g Hc=%g h=%g k=%g l=%g E=%g",numbers[4],numbers[1],numbers[2],numbers[3],numbers[5],numbers[6],numbers[7],numbers[9]);
                  hkl(1)=numbers[5];hkl(2)=numbers[6];hkl(3)=numbers[7];//E=numbers[9]; 

            switch(argv[1][1]) // dimension definition from jjjpar.hpp
            {case 's': 
             case 'o': 
                      if(doijk==3){// moments=xx*momentsx+yy*momentsy+zz*momentsz;
                   for(ii=1;ii<=inputpars.nofatoms;++ii)
                   for (i=1;i<=savmf.na();++i)for(j=1;j<=savmf.nb();++j)for(k=1;k<=savmf.nc();++k)
                  for(nt=1;nt<=dim;++nt){densityev_real.m(i,j,k)(nt+dim*(ii-1))=xx*ev_real.m(i,j,k)(nt+3*dim*(ii-1))+yy*ev_real.m(i,j,k)(nt+dim+3*dim*(ii-1))+zz*ev_real.m(i,j,k)(nt+2*dim+3*dim*(ii-1));
                                         densityev_imag.m(i,j,k)(nt+dim*(ii-1))=xx*ev_imag.m(i,j,k)(nt+3*dim*(ii-1))+yy*ev_imag.m(i,j,k)(nt+dim+3*dim*(ii-1))+zz*ev_imag.m(i,j,k)(nt+2*dim+3*dim*(ii-1));
                                        }
                      }
                      else{densityev_real=ev_real;
                           densityev_imag=ev_imag;
                     }
                       break;
            case 'c': 
            case 'j': densityev_real=ev_real;
                      densityev_imag=ev_imag;
                     break;
            default: help_and_exit();
                     }
                    }
               }
              fclose (fin);
             if(checkdd>1e-7)
               {fprintf(stderr,"Error program spins - inconsistent output files mcphas.mf and mcdisp.*:  temperature/magnetic field in static configuration\n" 
                               "(from results/mcphas.mf) is different from mcdisp output files results/mcdisp.qee,qsd,qod:\n %s\n"
                               "... probably you need to rerun setup_mcdisp_mf and mcdisp with T=%g Ha=%g Hb=%g Hc=%g",outstr,T,Hext(1),Hext(2),Hext(3));exit(EXIT_FAILURE);}
              fprintf(stdout,"#%s - density oscillation - eigenvector\n",outstr);
              fprintf(stdout,"#real\n");
              densityev_real.print(stdout);
              fprintf(stdout,"#imag\n");
              densityev_imag.print(stdout);
             }//gp.show_density
//----------------------------------------------------------------------------------------------------------           
              gp.read();// read graphic parameters which are set by user in file results/graphic_parameters.set
                        // in case he wants to overwrite some default settings
              if(arrow==0)gp.spins_scale_moment=0;
              if(density==0)gp.show_density=0;
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              double phase;
              complex <double> im(0,1);
              for(i=0;i<16;++i)
              {phase=2*3.1415*i/15;
               printf("\n********************************************\n");
               printf(" calculating movie sequence %i(16)\n",i+1);
               printf("********************************************\n");
               char filename[MAXNOFCHARINLINE];
               sprintf(filename,"./results/spins.%i.jvx",i+1);
               fin = fopen_errchk (filename, "w");gp.showprim=0;
                     densitycf.jvx_cd(fin,outstr,cs,gp,
                                  phase,densityev_real,densityev_imag,hkl,T,hh,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag);
               fclose (fin);
               sprintf(filename,"./results/spins_prim.%i.jvx",i+1);
               fin = fopen_errchk (filename, "w");gp.showprim=1;
                     densitycf.jvx_cd(fin,outstr,cs,gp,
                                  phase,densityev_real,densityev_imag,hkl,T,hh,Hext,cs4,spinconf,spinconfev_real,spinconfev_imag);
               fclose (fin);
              }
          printf("# %s\n",outstr);
          }
fprintf(stderr,"# ************************************************************************\n");
fprintf(stderr,"# *             end of program spins\n");
fprintf(stderr,"# * Reference: M. Rotter PRB 79 (2009) 140405R\n");
fprintf(stderr,"# * \n");
fprintf(stderr,"# * view jvx file by:\n");
fprintf(stderr,"# * javaview results/spins.jvx\n");
fprintf(stderr,"# * java javaview \"model=results/spins.*.jvx\" Animation.LastKey=16 background=\"255 255 255\" \n");
fprintf(stderr,"# * saved density mesh in results/spins.grid\n");
fprintf(stderr,"# ************************************************************************\n");

  for(i=1;i<=cs.nofatoms;++i){  delete cs.sipffilenames[i];}
  return 0;
}


