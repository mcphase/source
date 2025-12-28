

//output for javaview *******************************************************************
void spincf::jvx_cd(FILE * fout,char * text,cryststruct & cs,graphic_parameters & gp,
                    double phase,spincf & densityev_real,spincf & densityev_imag,Vector & hkl,double & T, Vector &  gjmbHxc,
                    Vector & Hext,cryststruct & cs4,
                    spincf & magmom,spincf & magmomev_real, spincf & magmomev_imag,
                    spincf & phonon,spincf & pev_real, spincf & pev_imag)
{ int l;
 // some checks
 if(nofatoms!=densityev_real.nofatoms||nofa!=densityev_real.na()||nofb!=densityev_real.nb()||nofc!=densityev_real.nc()||
    nofatoms!=densityev_imag.nofatoms||nofa!=densityev_imag.na()||nofb!=densityev_imag.nb()||nofc!=densityev_imag.nc()||
    nofcomponents<densityev_real.nofcomponents||densityev_real.nofcomponents!=densityev_imag.nofcomponents)
    {fprintf(stderr,"Error creating jvx movie files: eigenvector read from .qev file dim %i does not match  dimension %i of spins structure read from sps file\n",nofcomponents,densityev_real.nofcomponents);exit(1);}

  Vector maxv(1,3),minv(1,3),ijkmax(1,3),ijkmin(1,3),max_min(1,3),dd(1,3),dd0(1,3),c(1,3),xyz(1,3);
  Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
  Matrix abc_in_ijk_Inverse(1,3,1,3); abc_in_ijk_Inverse=abc_in_ijk.Inverse();
  
  calc_minmax_scale_relabc(minv,maxv,ijkmin,ijkmax,cs.r,cs.abc,gp.scale_view_1,gp.scale_view_2,gp.scale_view_3);
   if(gp.showprim==1){ijkmin(1)=1;ijkmin(2)=1;ijkmin(3)=1;ijkmax(1)=-1;ijkmax(2)=-1;ijkmax(3)=-1;} // show only primitive magnetic unit cell
  max_min=maxv-minv;
  int * pl;  pl=new int[magmom.nofatoms+1];
  for(l=1;l<=magmom.nofatoms;++l)pl[l]=l;

fprintf(fout,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n");
fprintf(fout,"<jvx-model>\n");
fprintf(fout,"  <title>%s</title>\n",text);
fprintf(fout,"  <geometries>\n");

if(gp.show_abc_unitcell>0)jvx_show_abc_unitcell(fout,gp,abc_in_ijk);

if(gp.show_primitive_crystal_unitcell>0)jvx_show_primitive_crystal_unitcell(fout,gp,cs);

if(gp.show_magnetic_unitcell>0)jvx_show_magnetic_unitcell(fout,gp,cs);

if(gp.show_atoms>0)jvx_show_atoms(fout,gp,cs,ijkmin,ijkmax,hkl,maxv,minv,abc_in_ijk_Inverse,phase,phonon,pev_real,pev_imag,pl);

if(gp.spins_scale_moment>0)jvx_show_magnetic_moments(fout,gp,cs,cs4,ijkmin,ijkmax,hkl,maxv,minv,abc_in_ijk_Inverse,
                            phase,magmom,magmomev_real,magmomev_imag,phonon,pev_real,pev_imag,pl);

if(gp.spins_show_static_moment_direction>0)jvx_show_static_magnetic_moments(fout,gp,cs,cs4,ijkmin,ijkmax,hkl,maxv,minv,abc_in_ijk_Inverse,
                             phase,magmom,phonon,pev_real,pev_imag,pl);

if(gp.spins_show_ellipses>0)jvx_spins_show_ellipses(fout,gp,cs,cs4,ijkmin,ijkmax,hkl,maxv,minv,abc_in_ijk_Inverse,
                             phase,magmom,magmomev_real,magmomev_imag,phonon,pev_real,pev_imag,pl);

if(gp.scale_density_vectors>0)jvx_density_vectors(fout,gp,cs,ijkmin,ijkmax,hkl,maxv,minv,abc_in_ijk_Inverse,
                             phase,magmom,densityev_real,densityev_imag,phonon,pev_real,pev_imag,pl,T,gjmbHxc,Hext);

if(gp.show_density>0)jvx_density(fout,gp,cs,ijkmin,ijkmax,hkl,maxv,minv,abc_in_ijk_Inverse,
                             phase,magmom,densityev_real,densityev_imag,phonon,pev_real,pev_imag,pl,T,gjmbHxc,Hext);

fprintf(fout,"  </geometries>\n");
fprintf(fout,"</jvx-model>\n");
delete []pl;
}


//output for javaview (old) *******************************************************************
void spincf::jvx_cd(FILE * fout,char * text,cryststruct & cs,graphic_parameters & gp,
                    double phase,spincf & densityev_real,spincf & densityev_imag,Vector & hkl,double & T, Vector &  gjmbHxc,
                    Vector & Hext,spincf & magmom,spincf & magmomev_real, spincf & magmomev_imag)
{spincf phonon(magmom.na(),magmom.nb(),magmom.nc(),magmom.nofatoms,3);phonon=phonon * 0.0;
 spincf pev_r(magmomev_real.na(),magmomev_real.nb(),magmomev_real.nc(),magmomev_real.nofatoms,3);pev_r=pev_r * 0.0;
 spincf pev_i(magmomev_imag.na(),magmomev_imag.nb(),magmomev_imag.nc(),magmomev_imag.nofatoms,3);pev_i=pev_i * 0.0;
 jvx_cd(fout,text,cs,gp,phase,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hext,cs,magmom,magmomev_real,magmomev_imag,phonon,pev_r,pev_i);}

//output for javaview *******************************************************************
void spincf::jvx_cd(FILE * fout,char * text,cryststruct & cs,graphic_parameters & gp,
                    double phase,spincf  densityev_real,spincf  densityev_imag,Vector & hkl,double & T, Vector &  gjmbHxc,
                    Vector & Hext,cryststruct & cs4,spincf  magmom,spincf  magmomev_real, spincf  magmomev_imag,spincf  phonon)
{//spincf phonon(magmom.na(),magmom.nb(),magmom.nc(),magmom.nofatoms,3);//phonon=phonon * 0.0;
 spincf pev_r(magmomev_real.na(),magmomev_real.nb(),magmomev_real.nc(),magmomev_real.nofatoms,3);//pev_r=pev_r * 0.0;
 spincf pev_i(magmomev_imag.na(),magmomev_imag.nb(),magmomev_imag.nc(),magmomev_imag.nofatoms,3);//pev_i=pev_i * 0.0;
// fprintf(stderr,"densityev_real.nofcomponents=%i\n",densityev_real.nofcomponents);
 jvx_cd(fout,text,cs,gp,phase,densityev_real,densityev_imag,hkl,T,gjmbHxc,Hext,cs4,magmom,magmomev_real,magmomev_imag,phonon,pev_r,pev_i);}


void spincf::calc_minmax_scale_relabc(Vector & minv,Vector & maxv,Vector & ijkmin,Vector & ijkmax,Matrix & r,Vector & abc,double scale_view_1,double scale_view_2,double scale_view_3)
{// determine maxv(1,2,3) minv(1,2,3) (vector in units of A direction of abc describing
 //a parallelepiped) for viewing magnetic unit cell
  Matrix p(1,3,1,3);
  Vector nofabc(1,3);nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  Matrix rijk(1,3,1,3);dadbdc2ijk(rijk,r,abc);
  for (int i=1;i<=3;++i)for(int j=1;j<=3;++j)p(i,j)=nofabc(j)*r(i,j)*abc(i);
   
  Vector ddd(1,8),dd0(1,3),dd(1,3);
  int i;
 double t;
  for (i=1;i<=3;++i)
  {ddd(1)=p.Column(1)(i);
   ddd(2)=p.Column(2)(i);
   ddd(3)=p.Column(3)(i);
   ddd(4)=p.Column(1)(i)+p.Column(2)(i);
   ddd(5)=p.Column(1)(i)+p.Column(3)(i);
   ddd(6)=p.Column(2)(i)+p.Column(3)(i);
   ddd(7)=0;
   ddd(8)=p.Column(1)(i)+p.Column(2)(i)+p.Column(3)(i);
   minv(i)=Min(ddd);maxv(i)=Max(ddd);
   t=minv(i)/abc(i);if(abs(t-int(t))>0.0001){minv(i)=(int(t)-1.0)*abc(i);}
   t=maxv(i)/abc(i);if(abs(t-int(t))>0.0001){maxv(i)=(int(t)+1.0)*abc(i);}
  }
  maxv(1)=minv(1)+(maxv(1)-minv(1))*scale_view_1;
  maxv(2)=minv(2)+(maxv(2)-minv(2))*scale_view_2;
  maxv(3)=minv(3)+(maxv(3)-minv(3))*scale_view_3;
  // determine ijkmin ijkmax by calculating the 8 corners of the  quader
  // in terms of primitive lattice
  // i*p.Column(1)+j*p.Column(2)+k*p.Column(3)=cornerpointvector ... i,j,k =?
  // ijk=p^-1*corerpointvector
  for (i=1;i<=3;++i)
  {dd0=minv;               dd=p.Inverse()*dd0;ddd(1)=dd(i);
   dd0=minv;dd0(1)=maxv(1);dd=p.Inverse()*dd0;ddd(2)=dd(i);
   dd0=minv;dd0(2)=maxv(2);dd=p.Inverse()*dd0;ddd(3)=dd(i);
   dd0=minv;dd0(3)=maxv(3);dd=p.Inverse()*dd0;ddd(4)=dd(i);
   dd0=maxv;               dd=p.Inverse()*dd0;ddd(5)=dd(i);
   dd0=maxv;dd0(1)=minv(1);dd=p.Inverse()*dd0;ddd(6)=dd(i);
   dd0=maxv;dd0(2)=minv(2);dd=p.Inverse()*dd0;ddd(7)=dd(i);
   dd0=maxv;dd0(3)=minv(3);dd=p.Inverse()*dd0;ddd(8)=dd(i);
   ijkmin(i)=Min(ddd);ijkmax(i)=Max(ddd);
  }
  for(i=1;i<=3;++i){minv(i)/=abc(i);maxv(i)/=abc(i);}
}

// creates a box corresponding to the abc unit cell  ********************************************************
void spincf::jvx_show_abc_unitcell(FILE * fout,graphic_parameters & gp,Matrix & abc_in_ijk)
 {  // plot frame around crystallographic unit cell
fprintf(fout,"    <geometry name=\"crystallographic unit cell\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
fprintf(fout,"          <p>  %g       %g       %g </p>\n",0.0,0.0,0.0);
fprintf(fout,"          <p name=\"a\">  %g       %g       %g </p>\n",abc_in_ijk(1,1),abc_in_ijk(2,1),abc_in_ijk(3,1));
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc_in_ijk(1,1)+abc_in_ijk(1,2),abc_in_ijk(2,1)+abc_in_ijk(2,2),abc_in_ijk(3,1)+abc_in_ijk(3,2));
fprintf(fout,"          <p name=\"b\"> %g       %g       %g </p>\n",abc_in_ijk(1,2),abc_in_ijk(2,2),abc_in_ijk(3,2));
fprintf(fout,"          <p name=\"c\"> %g       %g       %g </p>\n",abc_in_ijk(1,3),abc_in_ijk(2,3),abc_in_ijk(3,3));
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc_in_ijk(1,1)+abc_in_ijk(1,3),abc_in_ijk(2,1)+abc_in_ijk(2,3),abc_in_ijk(3,1)+abc_in_ijk(3,3));
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc_in_ijk(1,1)+abc_in_ijk(1,2)+abc_in_ijk(1,3),abc_in_ijk(2,1)+abc_in_ijk(2,2)+abc_in_ijk(2,3),abc_in_ijk(3,1)+abc_in_ijk(3,2)+abc_in_ijk(3,3));
fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",abc_in_ijk(1,2)+abc_in_ijk(1,3),abc_in_ijk(2,2)+abc_in_ijk(2,3),abc_in_ijk(3,2)+abc_in_ijk(3,3));
fprintf(fout,"          <thickness>0.0</thickness>\n");
fprintf(fout,"          <colorTag type=\"rgb\">255 0 0</colorTag>\n");
fprintf(fout,"			<labelAtt horAlign=\"head\" visible=\"show\" font=\"fixed\" verAlign=\"top\">\n");
fprintf(fout,"				<xOffset>0</xOffset>\n");
fprintf(fout,"				<yOffset>0</yOffset>\n");
fprintf(fout,"			</labelAtt>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\">\n");
fprintf(fout,"        <lines>\n");
fprintf(fout,"          <l>0 1</l>\n");
fprintf(fout,"          <l>1 2</l>\n");
fprintf(fout,"          <l>2 3</l>\n");
fprintf(fout,"          <l>3 0</l>\n");
fprintf(fout,"          <l>0 4</l>\n");
fprintf(fout,"          <l>1 5</l>\n");
fprintf(fout,"          <l>2 6</l>\n");
fprintf(fout,"          <l>3 7</l>\n");
fprintf(fout,"          <l>4 5</l>\n");
fprintf(fout,"          <l>5 6</l>\n");
fprintf(fout,"          <l>6 7</l>\n");
fprintf(fout,"          <l>7 4</l>\n");
fprintf(fout,"          <thickness>1.0</thickness>\n");
fprintf(fout,"        <color type=\"rgb\">%i 0 0</color>\n",(int)(255*gp.show_abc_unitcell));
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
 }

// creates a box corresponding to the primitive unit cell   ********************************************************
void spincf::jvx_show_primitive_crystal_unitcell(FILE * fout,graphic_parameters & gp,cryststruct & cs)
 { Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
   Vector dd(1,3);
 // plot frame around primitive crystallographic unit cell
fprintf(fout,"    <geometry name=\"primitive crystallographic unit cell\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
dd=0;           fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\"r1\">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(2)/(double)nofb;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd-=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\"r2\">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd =p.Column(3)/(double)nofc;fprintf(fout,"          <p name=\"r3\">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(2)/(double)nofb;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd-=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
fprintf(fout,"			<labelAtt horAlign=\"head\" visible=\"show\" font=\"fixed\" verAlign=\"bottom\">\n");
fprintf(fout,"				<xOffset>0</xOffset>\n");
fprintf(fout,"				<yOffset>0</yOffset>\n");
fprintf(fout,"                          <colorTag type=\"rgb\">0 255 0</colorTag>\n");
fprintf(fout,"			</labelAtt>\n");
fprintf(fout,"          <thickness>0.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n");
fprintf(fout,"        <lines>\n");
fprintf(fout,"          <l>0 1</l>\n");
fprintf(fout,"          <l>1 2</l>\n");
fprintf(fout,"          <l>2 3</l>\n");
fprintf(fout,"          <l>3 0</l>\n");
fprintf(fout,"          <l>0 4</l>\n");
fprintf(fout,"          <l>1 5</l>\n");
fprintf(fout,"          <l>2 6</l>\n");
fprintf(fout,"          <l>3 7</l>\n");
fprintf(fout,"          <l>4 5</l>\n");
fprintf(fout,"          <l>5 6</l>\n");
fprintf(fout,"          <l>6 7</l>\n");
fprintf(fout,"          <l>7 4</l>\n");
fprintf(fout,"          <thickness>1.0</thickness>\n");
fprintf(fout,"        <color type=\"rgb\">0 %i 0</color>\n",(int)(255*gp.show_primitive_crystal_unitcell));
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
}

  // creates a box corresponding to the magnetic unit cell********************************************************
void spincf::jvx_show_magnetic_unitcell(FILE * fout,graphic_parameters & gp,cryststruct & cs)
 {Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
Vector dd(1,3);
  // plot frame around primitive magnetic unit cell
fprintf(fout,"    <geometry name=\"magnetic unit cell\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n");
fprintf(fout,"        <points>\n");
dd=0;           fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(2);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd-=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd =p.Column(3);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd+=p.Column(2);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
dd-=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n");
fprintf(fout,"        <lines>\n");
fprintf(fout,"          <l>0 1</l>\n");
fprintf(fout,"          <l>1 2</l>\n");
fprintf(fout,"          <l>2 3</l>\n");
fprintf(fout,"          <l>3 0</l>\n");
fprintf(fout,"          <l>0 4</l>\n");
fprintf(fout,"          <l>1 5</l>\n");
fprintf(fout,"          <l>2 6</l>\n");
fprintf(fout,"          <l>3 7</l>\n");
fprintf(fout,"          <l>4 5</l>\n");
fprintf(fout,"          <l>5 6</l>\n");
fprintf(fout,"          <l>6 7</l>\n");
fprintf(fout,"          <l>7 4</l>\n");
fprintf(fout,"          <thickness>%g</thickness>\n",gp.show_magnetic_unitcell);
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
}

// creates points at atom positions ********************************************************************
void spincf::jvx_show_atoms(FILE * fout,graphic_parameters & gp,cryststruct & cs,Vector & ijkmin, Vector & ijkmax,
Vector & hkl,Vector &maxv,Vector &minv,Matrix & abc_in_ijk_Inverse,
 double & phase,spincf & phonon,spincf & pev_real,spincf & pev_imag,int * pl)
 { Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
int n=0,ctr=0; // limits the  number of atoms to be plotted   
 // plot atoms in region xmin to xmax (quader)
fprintf(fout,"    <geometry name=\"ions\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
 int i1,j1,k1,i,j,k,l; Vector dd(1,3),dd0(1,3);
  for (i1=int(ijkmin(1)-1.0);(i1<=int(ijkmax(1)+1))&(n<NOMORE);++i1){
   for (j1=int(ijkmin(2)-1.0);(j1<=int(ijkmax(2)+1))&(n<NOMORE);++j1){
    for (k1=int(ijkmin(3)-1.0);(k1<=int(ijkmax(3)+1))&(n<NOMORE);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;(i<=nofa)&(n<NOMORE);++i){for (j=1;(j<=nofb)&(n<NOMORE);++j){for (k=1;(k<=nofc)&(n<NOMORE);++k){
         for(l=1;(l<=nofatoms)&(n<NOMORE);++l)
	 {dd=pos(i,j,k,l,cs); int showdd=1;// the following is to remove atoms too close to
            // each other (because e.g. phonon and so1ion modules 'sit' at nearly the same position 
            // and refer to the same atom
         for(int ll=1;ll<=nofatoms;++ll){Vector ddd(1,3);ddd=pos(i,j,k,ll,cs);
            if(Norm(dd-ddd)<0.2&&ll!=l){
                                        // here we decide what to do with two ions at the same position
                                        // static: magmom(1..3)  moment (arrow)
                                        //         phonon (1..3) displacement (shift of position)
                                        // dynamic: pev_real,imag contains nuclear movement
                                        //          magmomev_real,imag   moment oscillation eigenvector
                                        //          densityev_real,imag  chargedensity oscillation eigenvector

                                        // if neighbour ll is displacement (radius!=0) and  l is moment - 
                                       // do not show point and set oscillation of l  equal to ll
                                         double radius=0;extract(cs.sipffilenames[ll],"radius",radius);
                                         if(radius!=0){radius=0;
                                         extract(cs.sipffilenames[l],"radius",radius);
                                         if(radius==0){//printf("%i %i %i taking for magnetic atom nr %i position from atom nr %i\n",i,j,k,l,ll);
                                                   pl[l]=ll;
                                                   showdd=0;}}
                                       }

            }

         dd+=dd0;
            if(showdd==1)if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
            (gp.showprim==1&&i<=1+(nofa-1)*gp.scale_view_1&&j<=1+(nofb-1)*gp.scale_view_2&&k<=1+(nofc-1)*gp.scale_view_3))
            {double QR;  QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI;
             Vector p(1,3),xyz(1,3);p=0;
             double radius=0;extract(cs.sipffilenames[l],"radius",radius);if(radius!=0){p=phonon.moment(i,j,k,l);}
             if(pl[l]!=l)p=phonon.moment(i,j,k,pl[l]); 
             xyz=gp.phonon_scale_static_displacements*p;
             xyz+=gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,l)+sin(phase-QR)*pev_imag.moment(i,j,k,l));
//printf("i1=%i j1=%i k1=%i i=%i j=%i k=%i l=%i pl=%i dd=%g %g %g xyz = %g %g %g \n ",i1,j1,k1,i,j,k,l,pl[l],dd(1),dd(2),dd(3),xyz(1),xyz(2),xyz(3));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));             
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)+xyz(1)),myround(dd(2)+xyz(2)),myround(dd(3)+xyz(3)));
     ++ctr;++n;
	     }
	  }
       }}}
  }}}
fprintf(fout,"          <thickness>3.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"        <colors type=\"rgb\">\n");n=0;
  for (i1=int(ijkmin(1)-1.0);(i1<=int(ijkmax(1)+1))&(n<NOMORE);++i1){
   for (j1=int(ijkmin(2)-1.0);(j1<=int(ijkmax(2)+1))&(n<NOMORE);++j1){
    for (k1=int(ijkmin(3)-1.0);(k1<=int(ijkmax(3)+1))&(n<NOMORE);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;(i<=nofa)&(n<NOMORE);++i){for (j=1;(j<=nofb)&(n<NOMORE);++j){for (k=1;(k<=nofc)&(n<NOMORE);++k){
         for(l=1;(l<=nofatoms)&(n<NOMORE);++l)
	   {dd=pos(i,j,k,l,cs); int showdd=1;// the following is to remove atoms too close to
            // each other (because e.g. phonon and so1ion modules 'sit' at nearly the same position 
            // and refer to the same atom
         for(int ll=1;ll<=nofatoms;++ll){Vector ddd(1,3);ddd=pos(i,j,k,ll,cs);
            if(Norm(dd-ddd)<0.15&&ll!=l){ double radius=0;extract(cs.sipffilenames[ll],"radius",radius);
                                         if(radius!=0){radius=0;
                                         extract(cs.sipffilenames[l],"radius",radius);
                                         if(radius==0){
                                          showdd=0;}}
        }
            }
         dd+=dd0;
            if(showdd==1)if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
             (gp.showprim==1&&i<=1+(nofa-1)*gp.scale_view_1&&j<=1+(nofb-1)*gp.scale_view_2&&k<=1+(nofc-1)*gp.scale_view_3))
            {
fprintf(fout,"          <c>  %i       %i       %i </c>\n",255,255,255);
int r=(int)(255*gp.show_atoms),g=(int)(gp.show_atoms*((l*97)%256)),b=0;
 extract(cs.sipffilenames[l],"r",r);
 extract(cs.sipffilenames[l],"g",g);
 extract(cs.sipffilenames[l],"b",b);
// here should be rgb color of the ion !!
fprintf(fout,"          <c>  %i       %i       %i </c>\n",r,g,b);
//printf("          <c> %s  %i       %i       %i </c>\n",cs.sipffilenames[l],r,g,b);
	     }
	  }
       }}}
  }}}
fprintf(fout,"        </colors>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n");
fprintf(fout,"        <lines>\n");
  for(i=0;i<ctr;++i)fprintf(fout,"          <l>%i %i</l>\n",2*i,2*i+1);
fprintf(fout,"          <thickness>1.0</thickness>\n");
fprintf(fout,"          <colorTag type=\"rgb\">255 255 0</colorTag>\n");
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
if(n>=NOMORE){fprintf(stderr,"# Warning - jvx output truncated because maximum number of atoms %i reached - continuing\n",NOMORE);}

}


// creates arrows corresponding to magnetic moments (or ...)  ***************************************************
void spincf::jvx_show_magnetic_moments(FILE * fout,graphic_parameters & gp,cryststruct & cs,cryststruct & cs4,Vector & ijkmin, Vector & ijkmax,
          Vector & hkl,Vector &maxv,Vector &minv,Matrix & abc_in_ijk_Inverse,
          double & phase,spincf & magmom,spincf & magmomev_real,spincf & magmomev_imag,
          spincf & phonon,spincf & pev_real,spincf & pev_imag, int * pl)
{Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
 if(gp.phonon_scale_static_displacements!=0&&gp.show_atoms==0)
 {fprintf(stderr,"# jvx Warning: moment arrows cannot be drawn, because phonon_scale_static_displacements not zero and show_atoms=0\n");
 }
 else 
 {
fprintf(fout,"    <geometry name=\"magnetic moments\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
 int ctr=0,n=0,i1,j1,k1,i,j,k,l;Vector dd(1,3),dd0(1,3),xyz(1,3);
 for (i1=int(ijkmin(1)-1.0);(i1<=int(ijkmax(1)+1))&(n<NOMORE);++i1){
 for (j1=int(ijkmin(2)-1.0);(j1<=int(ijkmax(2)+1))&(n<NOMORE);++j1){
 for (k1=int(ijkmin(3)-1.0);(k1<=int(ijkmax(3)+1))&(n<NOMORE);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;(i<=nofa)&(n<NOMORE);++i){for (j=1;(j<=nofb)&(n<NOMORE);++j){for (k=1;(k<=nofc)&(n<NOMORE);++k){
         for(l=1;(l<=magmom.nofatoms)&(n<NOMORE);++l)
	 {//printf("i=%i j=%i k=%i l=%i i!=%i j1=%i k1=%i\n",i,j,k,l,i1,j1,k1);
          dd=magmom.pos(i,j,k,l, cs4);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
                   (gp.showprim==1&&i<=1+(nofa-1)*gp.scale_view_1&&j<=1+(nofb-1)*gp.scale_view_2&&k<=1+(nofc-1)*gp.scale_view_3))
            {double QR; QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI;
             xyz=magmom.moment(i,j,k,l);
             if(gp.spins_show_oscillation){xyz+=gp.spins_wave_amplitude*(cos(-phase+QR)*magmomev_real.moment(i,j,k,l)+sin(phase-QR)*magmomev_imag.moment(i,j,k,l));}
              //if(pl[l]!=l)dd+=gp.phonon_scale_static_displacements * phonon.moment(i,j,k,pl[l]);
               dd+=gp.phonon_scale_static_displacements * phonon.moment(i,j,k,pl[l]);
               dd+=gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,l)+sin(phase-QR)*pev_imag.moment(i,j,k,l));
              //printf("gJ=%g magmom=%g %g %g %g %g %g %g %g %g\n",cs.gJ[l],mom[in(i,j,k)](1),mom[in(i,j,k)](2),mom[in(i,j,k)](3),mom[in(i,j,k)](4),mom[in(i,j,k)](5),mom[in(i,j,k)](6),xyz(1),xyz(2),xyz(3));
              //if(l==170||l==171){fprintf(stderr,"l=%i\n %4.4f + i %4.4f\n %4.4f + i %4.4f\n %4.4f + i %4.4f\n",
              //                    l,magmomev_real.moment(i,j,k,l)(1),magmomev_imag.moment(i,j,k,l)(1),
              //                      magmomev_real.moment(i,j,k,l)(2),magmomev_imag.moment(i,j,k,l)(2),
              //                      magmomev_real.moment(i,j,k,l)(3),magmomev_imag.moment(i,j,k,l)(3)
              //                          );
                               
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(densityev_real*cos(-phase) + densityev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
//fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)-xyz(1)*gp.spins_scale_moment),myround(dd(2)-xyz(2)*gp.spins_scale_moment),myround(dd(3)-xyz(3)*gp.spins_scale_moment));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)+xyz(1)*gp.spins_scale_moment),myround(dd(2)+xyz(2)*gp.spins_scale_moment),myround(dd(3)+xyz(3)*gp.spins_scale_moment));
	     ++ctr;++n;            
              //                  }
	     }
	  }
       }}}
  }}}
fprintf(fout,"          <thickness>6.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"show\" line=\"show\">\n");
fprintf(fout,"        <lines>\n");
  for(i=0;i<ctr;++i)fprintf(fout,"          <l>%i %i</l>\n",2*i,2*i+1);
switch((int)gp.spins_colour)
{case 4:fprintf(fout,"        <color type=\"rgb\">150  0 0</color>\n");break;  // pel
 case 3:fprintf(fout,"        <color type=\"rgb\">0  200 150</color>\n");break;  // S
 case 2:fprintf(fout,"        <color type=\"rgb\">200 153 0</color>\n");break;  // L
 case 1:
 default:fprintf(fout,"        <color type=\"rgb\">0 150 0</color>\n");break; // magmom
}
fprintf(fout,"          <thickness>4.0</thickness>\n");
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
if(n>=NOMORE){fprintf(stderr,"# Warning - jvx output truncated because maximum number of atoms %i reached - continuing\n",NOMORE);}
                     }
}

// creates lines corresponding to static magnetic moments (or ...) *******************************************
void spincf::jvx_show_static_magnetic_moments(FILE * fout,graphic_parameters & gp,cryststruct & cs,cryststruct & cs4,Vector & ijkmin, Vector & ijkmax,
          Vector & hkl,Vector &maxv,Vector &minv,Matrix & abc_in_ijk_Inverse,
          double & phase,spincf & magmom,
          spincf & phonon,spincf & pev_real,spincf & pev_imag, int * pl)
 {Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
 if(gp.phonon_scale_static_displacements!=0&&gp.show_atoms==0)
 {fprintf(stderr,"# jvx Warning: static moment lines cannot be drawn, because phonon_scale_static_displacements not zero and show_atoms=0\n");
 }
 else 
 {
// plot a line along static magnetic moments for comparison
fprintf(fout,"    <geometry name=\"static magnetic moments\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n");
fprintf(fout,"        <points>\n");
 int ctr=0,n=0,i1,j1,k1,i,j,k,l;Vector xyz(1,3),dd0(1,3),dd(1,3);
 for (i1=int(ijkmin(1)-1.0);(i1<=int(ijkmax(1)+1))&(n<NOMORE);++i1){
 for (j1=int(ijkmin(2)-1.0);(j1<=int(ijkmax(2)+1))&(n<NOMORE);++j1){
 for (k1=int(ijkmin(3)-1.0);(k1<=int(ijkmax(3)+1))&(n<NOMORE);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;(i<=nofa)&(n<NOMORE);++i){for (j=1;(j<=nofb)&(n<NOMORE);++j){for (k=1;(k<=nofc)&(n<NOMORE);++k){
         for(l=1;(l<=magmom.nofatoms)&(n<NOMORE);++l)
	 {dd=magmom.pos(i,j,k,l, cs4);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
                     (gp.showprim==1&&i<=1+(nofa-1)*gp.scale_view_1&&j<=1+(nofb-1)*gp.scale_view_2&&k<=1+(nofc-1)*gp.scale_view_3))
            {double QR;  QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI; 
                          xyz=magmom.moment(i,j,k,l);
            //if(pl[l]!=l)dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
            dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
            dd+=gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,l)+sin(phase-QR)*pev_imag.moment(i,j,k,l));
             
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)),myround(dd(2)),myround(dd(3)));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)+xyz(1)*gp.spins_scale_moment),myround(dd(2)+xyz(2)*gp.spins_scale_moment),myround(dd(3)+xyz(3)*gp.spins_scale_moment));
	     ++ctr;++n;

	     }
	  }
       }}}
  }}}
fprintf(fout,"          <thickness>6.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n");
fprintf(fout,"        <lines>\n");
  for(i=0;i<ctr;++i)fprintf(fout,"          <l>%i %i</l>\n",2*i,2*i+1);
fprintf(fout,"          <thickness>%g</thickness>\n",gp.spins_scale_moment<3?gp.spins_scale_moment:3);
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
if(n>=NOMORE){fprintf(stderr,"# Warning - jvx output truncated because maximum number of atoms %i reached - continuing\n",NOMORE);}
 }
}

// creates ellipses along the change of magnetic moments (or ...)  ***************************************************
void spincf::jvx_spins_show_ellipses(FILE * fout,graphic_parameters & gp,cryststruct & cs,cryststruct & cs4,Vector & ijkmin, Vector & ijkmax,
          Vector & hkl,Vector &maxv,Vector &minv,Matrix & abc_in_ijk_Inverse,
          double & phase,spincf & magmom,spincf & magmomev_real,spincf & magmomev_imag,
          spincf & phonon,spincf & pev_real,spincf & pev_imag, int * pl)
{Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
if(gp.phonon_scale_static_displacements!=0&&gp.show_atoms==0)
 {fprintf(stderr,"# jvx Warning: ellipses cannot be drawn, because phonon_scale_static_displacements not zero and show_atoms=0\n");
 }
 else 
 {// plot an ellipse along path of moment
fprintf(fout,"    <geometry name=\"ellipses\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n");
fprintf(fout,"        <points>\n");
 int ctr=0,n=0,i,j,k,l,i1,j1,k1; Vector dd0(1,3),dd(1,3),xyz(1,3);
 for (i1=int(ijkmin(1)-1.0);(i1<=int(ijkmax(1)+1))&(n<NOMORE);++i1){
 for (j1=int(ijkmin(2)-1.0);(j1<=int(ijkmax(2)+1))&(n<NOMORE);++j1){
 for (k1=int(ijkmin(3)-1.0);(k1<=int(ijkmax(3)+1))&(n<NOMORE);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;(i<=nofa)&(n<NOMORE);++i){for (j=1;(j<=nofb)&(n<NOMORE);++j){for (k=1;(k<=nofc)&(n<NOMORE);++k){
         for(l=1;(l<=magmom.nofatoms)&(n<NOMORE);++l)
	 {dd=magmom.pos(i,j,k,l, cs4);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
                    (gp.showprim==1&&i<=1+(nofa-1)*gp.scale_view_1&&j<=1+(nofb-1)*gp.scale_view_2&&k<=1+(nofc-1)*gp.scale_view_3))
            {double QR; QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI;
             //if(pl[l]!=l)dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
             dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
             dd+=gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,l)+sin(phase-QR)*pev_imag.moment(i,j,k,l));
                          int phi;
             for(phi=0;phi<=16;phi++)
             {
             xyz=magmom.moment(i,j,k,l)+gp.spins_wave_amplitude*(cos(-(double)phi*2*3.1415/16+QR)*magmomev_real.moment(i,j,k,l)+sin((double)phi*2*3.1415/16-QR)*magmomev_imag.moment(i,j,k,l));
              // <Jalpha>(i)=<Jalpha>0(i)+gp.spins_wave_amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(densityev_real*cos(-phase) + densityev_imag*sin(phase))*gp.spins_wave_amplitude; // Q ri not considered for test !!!
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dd(1)+xyz(1)*gp.spins_scale_moment),myround(dd(2)+xyz(2)*gp.spins_scale_moment),myround(dd(3)+xyz(3)*gp.spins_scale_moment));
	     }++ctr;++n;

	     }
	  }
       }}}
  }}}
fprintf(fout,"          <thickness>6.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"hide\" line=\"show\" color=\"show\">\n");
fprintf(fout,"        <lines>\n");
  for(i=0;i<ctr;++i)for(j=0;j<16;++j)fprintf(fout,"          <l>%i %i</l>\n",17*i+j,17*i+j+1);
fprintf(fout,"          <thickness>1.0</thickness>\n");
fprintf(fout,"        <color type=\"rgb\">0 %i %i</color>\n",(int)(100*gp.spins_show_ellipses),(int)(100*gp.spins_show_ellipses));
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
if(n>=NOMORE){fprintf(stderr,"# Warning - jvx output truncated because maximum number of atoms %i reached - continuing\n",NOMORE);}
 }
}

// creates density vectors in primitive magnetic unit cell ***************************************************
void spincf::jvx_density_vectors(FILE * fout,graphic_parameters & gp,cryststruct & cs,Vector & ijkmin, Vector & ijkmax,
          Vector & hkl,Vector &maxv,Vector &minv,Matrix & abc_in_ijk_Inverse,
          double & phase,spincf & magmom,spincf & densityev_real,spincf & densityev_imag,
          spincf & phonon,spincf & pev_real,spincf & pev_imag, int * pl,double & T, Vector &  gjmbHxc,Vector & Hext)
{int ii;
  double dtheta=0.2; //stepwidth to step surface
  double dfi=0.2;
if(gp.phonon_scale_static_displacements!=0&&gp.show_atoms==0)
 {fprintf(stderr,"# jvx Warning: density vector arrows cannot be drawn, because phonon_scale_static_displacements not zero and show_atoms=0\n");
 }
 else 
 {
for(int l=1;l<=nofatoms;++l)
 {fprintf(fout,"    <geometry name=\"density vectors in primitive magnetic unit cell - atom %i\">\n",l);
  fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"show\">\n");
  fprintf(fout,"<points >\n");
  double radius=0;double dx,dy,dz,R,fi,theta;
  int ctr=0,i,j,k;
  extract(cs.sipffilenames[l],"radius",radius);
  if(radius==0) // this is a trick: if radius is given as sipffilename then a sphere with this is radius is generated (pointcharge)
  {jjjpar ionpar(cs.x[l],cs.y[l],cs.z[l],cs.sipffilenames[l],1);
   density cd(gp.title,dtheta,dfi);int ndd;
   for (i=1;i<=1+(nofa-1)*gp.scale_view_1;++i){for(j=1;j<=1+(nofb-1)*gp.scale_view_2;++j){for(k=1;k<=1+(nofc-1)*gp.scale_view_2;++k){
   Vector dd(1,3);dd=pos(i,j,k,l, cs);
   Vector moments(1,nofcomponents);//printf("nofcomp=%i\n",nofcomponents);
   double QR;  QR=(hkl*abc_in_ijk_Inverse)*dd;
   QR*=2*PI;
      //if(pl[l]!=l)dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
      dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
      dd+=gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,l)+sin(phase-QR)*pev_imag.moment(i,j,k,l));
                             for(ndd=1;ndd<=densityev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i,j,k,l)(ndd)+gp.spins_wave_amplitude*(cos(-phase+QR)*densityev_real.moment(i,j,k,l)(ndd)+sin(phase-QR)*densityev_imag.moment(i,j,k,l)(ndd));}
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(densityev_real*cos(-phase) + densityev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
 // here we calculate the chargedensity of ion
   cd.calc_cd_surface(moments,ionpar,gp.threshhold,T,gjmbHxc,Hext);
   for(ii=1;ii<=cd.nofpoints();++ii)
     {R=cd.rtf(ii)(1);theta=cd.rtf(ii)(2);fi=cd.rtf(ii)(3);
     if(ionpar.orientation==abc_yzx){// mind abc||yzx in module cfield
     dx=R*sin(theta)*sin(fi)+dd(1);dy=R*cos(theta)+dd(2);dz=R*sin(theta)*cos(fi)+dd(3);
                              }
     else
                              {// mind abc||xyz in other cases ...
     dx=R*sin(theta)*cos(fi)+dd(1);dy=R*sin(theta)*sin(fi)+dd(2);dz=R*cos(theta)+dd(3);
                              }
     
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dx),myround(dy),myround(dz));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",myround(dx+cd.rtf(ii)(4)*gp.scale_density_vectors),myround(dy+cd.rtf(ii)(5)*gp.scale_density_vectors),myround(dz+cd.rtf(ii)(6)*gp.scale_density_vectors));
    ++ctr; }
   }}}
  }
fprintf(fout,"          <thickness>6.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"      <lineSet  arrow=\"show\" line=\"show\" color=\"show\">\n");
fprintf(fout,"        <lines>\n");
  for(i=0;i<ctr;++i)fprintf(fout,"          <l>%i %i</l>\n",2*i,2*i+1);
if (strncmp(gp.title+14,"currdensity",10)==0){
  fprintf(fout,"<color type=\"rgb\">%i %i %i </color>\n",220,153,0);
                                             }
 else
    {fprintf(fout,"        <color type=\"rgb\">0 255 0</color>\n");}
fprintf(fout,"          <thickness>%g</thickness>\n",gp.scale_density_vectors);
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");

 }
}
}

// creates density surface in primitive magnetic unit cell ************************************************
void spincf::jvx_density(FILE * fout,graphic_parameters & gp,cryststruct & cs,Vector & ijkmin, Vector & ijkmax,
          Vector & hkl,Vector &maxv,Vector &minv,Matrix & abc_in_ijk_Inverse,
          double & phase,spincf & magmom,spincf & densityev_real,spincf & densityev_imag,
          spincf & phonon,spincf & pev_real,spincf & pev_imag, int * pl,double & T, Vector &  gjmbHxc,Vector & Hext)
{int ii,tt,ff,i,j,k,l;
 if(gp.phonon_scale_static_displacements!=0&&gp.show_atoms==0)
 {fprintf(stderr,"# jvx Warning: density surfaces cannot be drawn, because phonon_scale_static_displacements not zero and show_atoms=0\n");
 }
 else 
 { double dtheta=gp.density_dtheta; //stepwidth to step surface
  double dfi=gp.density_dfi;Vector dd(1,3);
for(l=1;l<=nofatoms;++l)
 {double radius=0;double dx,dy,dz,R,fi,theta;
  extract(cs.sipffilenames[l],"radius",radius);
  if(radius!=0) // this is a trick: if radius is given as sipffilename then a sphere with this is radius is generated (pointcharge)
  {   if(gp.show_pointcharges>0)
      { fprintf(fout,"    <geometry name=\"densities in primitive magnetic unit cell - atom %i\">\n",l);
  fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
  fprintf(fout,"<points >\n");
   double rp=abs(radius);
        for (i=1;i<=(1+(nofa-1)*gp.scale_view_1);++i){for(j=1;j<=(1+(nofb-1)*gp.scale_view_2);++j){for(k=1;k<=(1+(nofc-1)*gp.scale_view_3);++k){
        dd=pos(i,j,k,l, cs); 
      double QR; QR=(hkl*abc_in_ijk_Inverse)*dd;
      QR*=2*PI;
        dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,l)+gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,l)+sin(phase-QR)*pev_imag.moment(i,j,k,l));             
        for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
             theta=(double)tt*dtheta;fi=(double)ff*dfi;
             dx=rp*sin(theta)*cos(fi)+dd(1);dy=rp*sin(theta)*sin(fi)+dd(2);dz=rp*cos(theta)+dd(3);
             fprintf(fout,"<p>%4g %4g %4g</p>\n",myround(dx),myround(dy),myround(dz));
             if(tt==0){ff=(int)(2*3.1415/dfi+1);}
             }}
        }}}
  } }
  else
  {fprintf(fout,"    <geometry name=\"densities in primitive magnetic unit cell - atom %i\">\n",l);
  fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
  fprintf(fout,"<points >\n");
  jjjpar ionpar(cs.x[l],cs.y[l],cs.z[l],cs.sipffilenames[l],1);
   density cd(gp.title,dtheta,dfi);int ndd;
   for (i=1;i<=1+(nofa-1)*gp.scale_view_1;++i){for(j=1;j<=1+(nofb-1)*gp.scale_view_2;++j){for(k=1;k<=1+(nofc-1)*gp.scale_view_2;++k){
   dd=pos(i,j,k,l, cs); 
   Vector moments(1,nofcomponents);
   double QR; QR=(hkl*abc_in_ijk_Inverse)*dd;
   QR*=2*PI;//printf("dd=%g",Norm(dd));
   //if(pl[l]!=l)dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
   dd+=gp.phonon_scale_static_displacements*phonon.moment(i,j,k,pl[l]);
   dd+=gp.phonon_wave_amplitude*(cos(-phase+QR)*pev_real.moment(i,j,k,pl[l])+sin(phase-QR)*pev_imag.moment(i,j,k,pl[l]));
    for(ndd=1;ndd<=densityev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i,j,k,l)(ndd)+gp.spins_wave_amplitude*(cos(-phase+QR)*densityev_real.moment(i,j,k,l)(ndd)+sin(phase-QR)*densityev_imag.moment(i,j,k,l)(ndd));
   }
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(densityev_real*cos(-phase) + densityev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
 // here we calculate the chargedensity of ion
   cd.calc_cd_surface(moments,ionpar,gp.threshhold,T,  gjmbHxc,Hext);
   for(ii=1;ii<=cd.nofpoints();++ii)
     {R=cd.rtf(ii)(1);theta=cd.rtf(ii)(2);fi=cd.rtf(ii)(3);
     if(ionpar.orientation==abc_yzx){// mind abc||yzx in module cfield
     dx=R*sin(theta)*sin(fi)+dd(1);dy=R*cos(theta)+dd(2);dz=R*sin(theta)*cos(fi)+dd(3);
                              }
     else
                              {// mind abc||xyz in other cases ...
     dx=R*sin(theta)*cos(fi)+dd(1);dy=R*sin(theta)*sin(fi)+dd(2);dz=R*cos(theta)+dd(3);
                              }
     fprintf(fout,"<p>%4g %4g %4g</p>\n",myround(dx),myround(dy),myround(dz));
     }
   }}}
  }
 if(radius==0||gp.show_pointcharges>0)
  {fprintf(fout,"<thickness>0.0</thickness><color type=\"rgb\">255 0 0</color><colorTag type=\"rgb\">255 0 255</colorTag>\n");
    fprintf(fout,"</points>			</pointSet>\n");
    fprintf(fout,"<faceSet face=\"show\" edge=\"show\">\n");
    fprintf(fout,"<faces >\n");
    int offset=0;
    for(i=1;i<=(1+(nofa-1)*gp.scale_view_1)*(1+(nofb-1)*gp.scale_view_2)*(1+(nofc-1)*gp.scale_view_3);++i)
    {int ntt,nff,pointnr,ffnr,p1,p2,p3,p4;
    ntt=(int)(3.1415/dtheta);
    nff=(int)(2*3.1415/dfi);
    pointnr=ntt*(nff+1);
    ffnr=nff+1;
    for(tt=1;tt<=ntt;++tt){for(ff=0;ff<=nff;++ff){
    p1 = ff + 1 + (tt - 2) * ffnr+offset;
    p2 = ff + 2 + (tt - 2) * ffnr+offset;
    p3 = ff + 2 + (tt - 1) * ffnr+offset;
    p4 = ff + 1 + (tt - 1) * ffnr+offset;
    if (ff==nff){p3 = p3 - ffnr; p2 = p2 - ffnr;}
    if (tt==1) {p1 = offset; p2 = offset;}
    fprintf(fout,"<f> %i %i %i %i </f>\n",p1,p2,p3,p4);
    }}
    offset+=pointnr+1;
 }
 //fprintf(fout,"<color type=\"rgb\">100 230 255</color>\n");
 if(radius>0||(strncmp(gp.title,"divergence",10)==0&&gp.threshhold>0))
 {fprintf(fout,"<color type=\"rgb\"> 255 0 0</color>\n");}
 else if (radius<0||(strncmp(gp.title,"divergence",10)==0&&gp.threshhold<0))
 {fprintf(fout,"<color type=\"rgb\">0  0 255</color>\n");}
 else
 {
 if(strncmp(gp.title,"chargedensity",10)==0){
  fprintf(fout,"<color type=\"rgb\">%i %i %i </color>\n",0,(int)(gp.show_density*((l*97)%256)),(int)(255*gp.show_density));
                                             }
 else if (strncmp(gp.title+14,"currdensity",10)==0){
  fprintf(fout,"<color type=\"rgb\">%i %i %i </color>\n",200+(int)(gp.show_density*((l*97)%56)),153,0);
                                             }
 else
    {
  fprintf(fout,"<color type=\"rgb\">%i %i %i </color>\n",0,200+(int)(gp.show_density*((l*97)%56)),0);
                                             }
 }
 fprintf(fout,"<colorTag type=\"rgb\">255 0 255</colorTag>\n");
 fprintf(fout,"</faces></faceSet>\n");
 fprintf(fout,"    </geometry>\n");
 }
 }
 }
}