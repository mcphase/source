//-----------------------------------------------------------------------
// create eps output of spinconfiguration
void spincf::eps(FILE * fout) //print spinconfiguration to stream
{eps(fout,"no title");}

void spincf::eps(FILE * fout,const char * text ) //print spinconfiguration to stream
{//viewport [-1,1,-1,1] ... distribute spins on that
  int i,j,k,l,m;
  float x0,y0,compoffset,atomoffset;
  Vector a(1,2);
  Vector b(1,2);
  float scale,d;

  scale=0;
  for (i=1;i<=nofa;++i)
    for (j=1;j<=nofb;++j)
     for (k=1;k<=nofc;++k)
      for (l=1;l<=nofcomponents*nofatoms;++l)
      {if ((d=fabs(mom[in(i,j,k)](l)))>scale)scale=d;
      }
  scale=0.2/(scale+0.01)/(double)nofc;

  fprintf(fout,"%s!PS-Adobe-2.0 EPSF-2.0\n","%");
  fprintf(fout,"%sBoundingBox:0 0 549 %i\n","%%",150*nofcomponents*nofatoms);
  fprintf(fout,"%sTitle: %s\n","%%",text);
  fprintf(fout,"%sEndComments\n","%%");
  fprintf(fout,"/mm {72 mul 25.4 div} bind def\n");
  fprintf(fout,"/mx {1.2 add 70 mul mm} bind def\n");
  fprintf(fout,"/my {1.2 add 70 mul mm} bind def\n");
   fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
   fprintf(fout,"-0.89 mx 0.38 my moveto \n (r1) show \n");
   fprintf(fout,"-0.93 mx 0.46 my moveto \n (r2) show \n");
   fprintf(fout,"-0.99 mx 0.52 my moveto \n (r3) show \n");
   fprintf(fout,"2 mm 2 mm moveto \n (%s) show \n",text);

   a(1)=-0.98;a(2)=0.4;
   b(1)=a(1);b(2)=0.5;epsarrow(fout,a,b);
   b(1)=-0.94;b(2)=0.45;epsarrow(fout,a,b);
   b(1)=-0.92;b(2)=0.4;epsarrow(fout,a,b);

//  char ss='a';
  for(l=1;l<=nofatoms;++l)
  {atomoffset=(l-0)*2*nofcomponents/3;
   for(m=1;m<=nofcomponents;++m)
    {compoffset=atomoffset-1.4-0.6*(m-1);
     fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
     fprintf(fout,"-0.98 mx %g my moveto \n (J%c%i) show \n",compoffset,'a'-1+m,l);

     for (i=1;i<=nofa;++i)
      for (j=1;j<=nofb;++j)
       for (k=1;k<=nofc;++k)
       { x0=(2.0*(i-1)/nofa-1+1.7*j/nofb/nofa)*1.1+0.3;
         y0=(2.0*(k-1)/nofc-1+1.2*j/nofb/nofc)*0.2+compoffset;
         a(1)=x0;a(2)=y0;
         b(1)=x0;b(2)=mom[in(i,j,k)](nofcomponents*(l-1)+m)*scale+y0;
         epsarrow(fout,a,b);
       }
    }
  }
fprintf(fout,"showpage\n");
}
void spincf::epsarrow(FILE * fout,Vector x,Vector y)
 {  double l=0.15*Norm(y-x);
    Vector y1(1,2),y2(1,2),unn(1,2),upn(1,2);
    if (y==x){y(2)=x(2)+0.0001;}
    upn=(y-x)/Norm(y-x);unn(1)=upn(2);unn(2)=-upn(1);
    y1=y-l*upn-0.3*l*unn;
    y2=y-l*upn+0.3*l*unn;

  fprintf(fout,"0.7 setlinewidth\n");
  fprintf(fout,"%g mx  %g my moveto\n",x(1),x(2));
  fprintf(fout,"%g mx  %g my lineto\n",y(1),y(2));
  fprintf(fout,"%g mx  %g my lineto\n",y1(1),y1(2));
  fprintf(fout,"%g mx  %g my lineto\n",y(1),y(2));
  fprintf(fout,"%g mx  %g my lineto\n",y2(1),y2(2));
  fprintf(fout,"stroke\n");

  }

void spincf::calc_prim_mag_unitcell(Matrix & p,Vector & abc, Matrix & r)
{ int i,j;
  Vector nofabc(1,3);nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  Matrix rijk(1,3,1,3);
  dadbdc2ijk(rijk,r,abc);
  for (i=1;i<=3;++i)for(j=1;j<=3;++j)p(i,j)=nofabc(j)*rijk(i,j);// old: dd(j)=nofabc(j)*r(i,j)*abc(i);
  // pa=p.Column(1);  //primitive magnetic unit cell
 // pb=p.Column(2);
 // pc=p.Column(3);
}



void spincf::calc_minmax(Vector & minv,Vector & maxv,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc)
{calc_minmax_scale(minv,maxv,ijkmin,ijkmax,p,abc,1.0,1.0,1.0);
}

void spincf::calc_minmax_scale(Vector & minv,Vector & maxv,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc,double scale_view_1,double scale_view_2,double scale_view_3)
{// determine max(1,2,3) min(1,2,3) (vector in units of A direction of abc describing
 //a parallelepiped) for viewing magnetic unit cell
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
}


Vector spincf::xy(Vector xyz,int orientation,Vector minv,Vector maxv,float bbwidth,float bbheight)
 {Vector p(1,2);
  switch(orientation)
  {case 1: p(1)=(xyz(1)-minv(1))/(maxv(1)-minv(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(2)-minv(2))/(maxv(2)-minv(2))*bbheight*0.8+bbheight*0.15;
   break;
   case 2: p(1)=(xyz(1)-minv(1))/(maxv(1)-minv(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)-minv(3))/(maxv(3)-minv(3))*bbheight*0.8+bbheight*0.15;
    break;
   case 3: p(1)=(xyz(2)-minv(2))/(maxv(2)-minv(2))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)-minv(3))/(maxv(3)-minv(3))*bbheight*0.8+bbheight*0.15;
    break;
   case 4: p(1)=(xyz(1)+(xyz(3)-minv(3))*0.1-minv(1))/(maxv(1)+(maxv(3)-minv(3))*0.1-minv(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(2)+(xyz(3)-minv(3))*0.15-minv(2))/(maxv(2)+(maxv(3)-minv(3))*0.15-minv(2))*bbheight*0.8+bbheight*0.15;
    break;
   case 5: p(1)=(xyz(1)+(xyz(2)-minv(2))*0.1-minv(1))/(maxv(1)+(maxv(2)-minv(2))*0.1-minv(1))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)+(xyz(2)-minv(2))*0.15-minv(3))/(maxv(3)+(maxv(2)-minv(2))*0.15-minv(3))*bbheight*0.8+bbheight*0.15;
    break;
   case 6: p(1)=(xyz(2)+(xyz(1)-minv(1))*0.1-minv(2))/(maxv(2)+(maxv(1)-minv(1))*0.1-minv(2))*bbwidth*0.8+bbwidth*0.15;
           p(2)=(xyz(3)+(xyz(1)-minv(1))*0.15-minv(3))/(maxv(3)+(maxv(1)-minv(1))*0.15-minv(3))*bbheight*0.8+bbheight*0.15;
    break;
   default: p=0;
  }
  return p;
 }


void spincf::eps3d(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,int orientation, spincf & magmom)
 {// function to plot spins in a 3d manner
  // orientation:1 ab 2 ac 3 bc projection
  //             4 ab 5 ac 6 bc side view
  int i,j,k,l;char r1,r2,r3;

  Vector a(1,2);
  Vector b(1,2),c(1,3);
  double scale,d,bbheight,bbwidth;

 // determine scale factor of moments
  scale=0;
  for (i=1;i<=nofa;++i)
    for (j=1;j<=nofb;++j)
     for (k=1;k<=nofc;++k)
      for(l=1;l<=nofatoms;++l)
      {c=magmom.moment(i,j,k,l);
       if ((d=Norm(c))>scale)scale=d;
      }
  scale=0.5/(scale+0.01);



  // determine max(1,2,3) min(1,2,3) (vector in Angstroem describing a quader) for viewing magnetic unit cell
  Vector maxv(1,3),minv(1,3),dd(1,3),max_min(1,3);
  Vector ddd(1,8),xyz(1,3),dd0(1,3),ijkmax(1,3),ijkmin(1,3);

  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  calc_minmax(minv,maxv,ijkmin,ijkmax,p,abc);
  max_min=maxv-minv;

 //determine bounding box for  specific view
  bbwidth=700;
  switch(orientation)
       { case 1 :
    		 {bbheight=max_min(2)/max_min(1)*bbwidth;r1='a';r3='b';}
	         break;
	 case 2 :
	         {bbheight=max_min(3)/max_min(1)*bbwidth;r1='a';r3='c';}
                 break;
	 case 3 :
	         {bbheight=max_min(3)/max_min(2)*bbwidth;r1='b';r3='c';}
                 break;
	 case 4 :
	         {bbheight=(max_min(2)+max_min(3)*0.1)/(max_min(3)*0.15+max_min(1))*bbwidth;r1='a';r3='b';r2='c';}
                  break;
	 case 5 :
	         {bbheight=(max_min(3)+max_min(2)*0.1)/(max_min(2)*0.15+max_min(1))*bbwidth;r1='a';r3='c';r2='b';}
                  break;
	 case 6 :
	         {bbheight=(max_min(3)+max_min(1)*0.1)/(max_min(1)*0.15+max_min(2))*bbwidth;r1='b';r3='c';r2='a';}
	          break;
	 default:  return;
	}

  fprintf(fout,"%s!PS-Adobe-2.0 EPSF-2.0\n","%");
  if (abc(4)!=90||abc(5)!=90||abc(6)!=90)
  {fprintf(fout,"%sBoundingBox:0 0 %i %i","%%",(int)bbwidth,(int)bbheight);
   fprintf(fout,"%sTitle: Nonorthogonal Lattice - Postscript output not supported\n","%%");
   fprintf(fout,"%sEndComments\n","%%");
   fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
   fprintf(fout,"/mm {72 mul 25.4 div} bind def\n");
   fprintf(fout,"/mx {1.2 add 50 mul mm} bind def\n");
   fprintf(fout,"/my {-0.3 add 50 mul mm} bind def\n");
   fprintf(fout,"-0.7 mx 0.38 my moveto \n (Nonorthogonal Lattice - Postscript output not supported) show \n");
  }
  else
  {
  fprintf(fout,"%sBoundingBox:0 0 %i %i\n","%%",(int)bbwidth,(int)bbheight);
  fprintf(fout,"%sTitle: %s\n","%%",text);
  fprintf(fout,"%sEndComments\n","%%");
  fprintf(fout,"/mm {72 mul 25.4 div} bind def\n");
  fprintf(fout,"/mx {1.2 add 50 mul mm} bind def\n");
  fprintf(fout,"/my {-0.3 add 50 mul mm} bind def\n");


  // draw abc coordinate label
   fprintf(fout,"/Helvetica findfont\n15 scalefont setfont\n");
   fprintf(fout,"-0.89 mx 0.38 my moveto \n (%c) show \n",r1);
   fprintf(fout,"-0.99 mx 0.52 my moveto \n (%c) show \n",r3);
   a(1)=-0.98;a(2)=0.4;
   b(1)=a(1);b(2)=0.5;epsarrow(fout,a,b);
   b(1)=-0.92;b(2)=0.4;epsarrow(fout,a,b);
   if (orientation>3){ fprintf(fout,"-0.93 mx 0.46 my moveto \n (%c) show \n",r2);
                      b(1)=-0.94;b(2)=0.45;epsarrow(fout,a,b);
		     }
  fprintf(fout,"-0.7 mx 0.38 my moveto \n (%s) show \n",text);

  fprintf(fout,"/mx {} bind def\n");
  fprintf(fout,"/my {} bind def\n");

  // draw frame around min vs max   (quader)
  fprintf(fout,"0.3 setlinewidth\n");
   a=xy(minv,orientation, minv, maxv,bbwidth,bbheight);
   dd=minv;dd(1)=maxv(1);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=minv;dd(2)=maxv(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=minv;dd(3)=maxv(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(maxv,orientation, minv, maxv,bbwidth,bbheight);
   dd=maxv;dd(1)=minv(1);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=maxv;dd(2)=minv(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=maxv;dd(3)=minv(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=minv;dd(2)=maxv(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=minv;dd(1)=maxv(1);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=maxv;dd(2)=minv(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=minv;dd(3)=maxv(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=maxv;dd(1)=minv(1);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=minv;dd(2)=maxv(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));


  // draw frame around primitive unit cell
  fprintf(fout,"1 setlinewidth\n");
   dd=0;
   a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   b=xy(p.Column(1),orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   b=xy(p.Column(2),orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   b=xy(p.Column(3),orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(2)+p.Column(3);a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=p.Column(1)+p.Column(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(2)+p.Column(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(2);a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=p.Column(1);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(1)+p.Column(3);a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=p.Column(1);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(2)+p.Column(3);a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
   dd=p.Column(2);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));
   dd=p.Column(3);b=xy(dd,orientation, minv, maxv,bbwidth,bbheight);fprintf(fout,"%g %g moveto\n",a(1),a(2));fprintf(fout,"%g %g lineto\n stroke \n",b(1),b(2));


  // plot atoms and moments in region xmin to xmax (quader)
int i1,j1,k1;
//i2,k2,j2,i1true,j1true,k1true;
   fprintf(fout,"/Helvetica findfont\n %i scalefont setfont\n",(int)(1000/nofa/nofb/nofc+1));

//these lines do not work if primitive lattice angles are > 90 deg ...
//i1true=1;for (i1=0;i1true==1;++i1){i1true=0;for(i2=-1;i2<=1;i2+=2){if (i1==0){i2=2;}
//j1true=1;for (j1=0;j1true==1;++j1){j1true=0;for(j2=-1;j2<=1;j2+=2){if (j1==0){j2=2;}
//k1true=1;for (k1=0;k1true==1;++k1){k1true=0;for(k2=-1;k2<=1;k2+=2){if (k1==0){k2=2;}
//   dd0=p.Column(1)*(double)(i2*i1)+p.Column(2)*(double)(j2*j1)+p.Column(3)*(double)(k2*k1);
for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
//printf("%i %i %i %i %i %i %i %i %i\n",i1,j1,k1,(int)ijkmin(1),(int)ijkmin(2),(int)ijkmin(3),(int)ijkmax(1),(int)ijkmax(2),(int)ijkmax(3));
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd+=dd0;
	    if(dd(1)<=maxv(1)+0.0001&&dd(1)>=minv(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=maxv(2)+0.0001&&dd(2)>=minv(2)-0.0001&&
            dd(3)<=maxv(3)+0.0001&&dd(3)>=minv(3)-0.0001)
            {c=magmom.moment(i,j,k,l);
              xyz(1)=dd(1)+scale*c(1);
              xyz(2)=dd(2)+scale*c(2);
              xyz(3)=dd(3)+scale*c(3);
	      a=xy(xyz,orientation, minv, maxv,bbwidth,bbheight);
              xyz(1)=dd(1)-scale*c(1);
              xyz(2)=dd(2)-scale*c(2);
              xyz(3)=dd(3)-scale*c(3);
              b=xy(xyz,orientation, minv, maxv,bbwidth,bbheight);
              epsarrow(fout,a,b);
	     }
	  }
       }}}
 }}}
 }

fprintf(fout,"showpage\n");


 }

int check_atom_in_big_unitcell(Vector & dd,Vector & maxv1,Vector & minv1,Matrix  &abc_in_ijk_Inverse){
            Vector dd1(1,3); dd1=abc_in_ijk_Inverse*dd;
//	    Vector minv1(1,3); minv1=minv*abc_in_ijk_Inverse;
//            Vector maxv1(1,3); maxv1=maxv*abc_in_ijk_Inverse;
        if((dd1(1)<=maxv1(1)+0.0001&&dd1(1)>=minv1(1)-0.0001&&   //if atom is in big unit cell
            dd1(2)<=maxv1(2)+0.0001&&dd1(2)>=minv1(2)-0.0001&&
            dd1(3)<=maxv1(3)+0.0001&&dd1(3)>=minv1(3)-0.0001))
           {return 1;}
           else
           {return 0;}
}

//***********************************************************************************************************************************
// output of chargedensity on grid as ascii file points are equally spaced as specified
// nofpoints*
void spincf::cd(FILE * fout,cryststruct & cs, graphic_parameters & gp,
                spincf & densityev_real,spincf & densityev_imag,double phase,Vector & hkl,double & T, Vector &  gjmbHxc,Vector & Hext)
{// some checks
 if(nofatoms!=densityev_real.nofatoms||nofa!=densityev_real.na()||nofb!=densityev_real.nb()||nofc!=densityev_real.nc()||
    nofatoms!=densityev_imag.nofatoms||nofa!=densityev_imag.na()||nofb!=densityev_imag.nb()||nofc!=densityev_imag.nc()||
    nofcomponents<densityev_real.nofcomponents||densityev_real.nofcomponents!=densityev_imag.nofcomponents)
    {fprintf(stderr,"Error creating density grid: eigenvector read from .qev file does not match dimension of spins structure read from sps file\n");exit(1);}
  int nofpointsi=gp.gridi;int nofpointsj=gp.gridj; int nofpointsk=gp.gridk;

  Vector maxv(1,3),minv(1,3),ijkmax(1,3),ijkmin(1,3),max_min(1,3),dd(1,3),dd0(1,3),c(1,3),xyz(1,3);
  Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
  Matrix abc_in_ijk_Inverse(1,3,1,3); abc_in_ijk_Inverse=abc_in_ijk.Inverse();
  Matrix p(1,3,1,3); calc_prim_mag_unitcell(p,cs.abc,cs.r);
  Matrix p_inverse (1,3,1,3); p_inverse=p.Inverse();
  calc_minmax_scale(minv,maxv,ijkmin,ijkmax,p,cs.abc,gp.scale_view_1,gp.scale_view_2,gp.scale_view_3);
  max_min=maxv-minv;

  int i,j,k,i1,j1,k1,imin,imax,jmin,jmax,kmin,kmax;Vector rijk(1,3);
  double *ro;ro=new double[nofpointsi*nofpointsj*nofpointsk];
  for(i=0;i<=nofpointsi*nofpointsj*nofpointsk-1;++i)ro[i]=0;
    // calculate density contribution of each ion around
  int l;
  for(l=1;l<=nofatoms;++l)
  {
  double radius=0;                            
  extract(cs.sipffilenames[l],"radius",radius);
  if(radius!=0) // this is a trick: if radius is given as sipffilename then a sphere with this is radius is generated (pointcharge)
  { if(gp.show_pointcharges>0)
    { double rp=abs(radius);
        // here we should introduce another loop to go around +-1 around the primitive
   // magnetic unit cell so that we see also atoms at the borders in the density map:
       int i0,j0,k0;
       for(i0=-1;i0<=1;++i0){for(j0=-1;j0<=1;++j0){for(k0=-1;k0<=1;++k0){
       for (i1=1;i1<=nofa;++i1){for(j1=1;j1<=nofb;++j1){for(k1=1;k1<=nofc;++k1){
        dd0=pos(i0*nofa+i1,j0*nofb+j1,k0*nofc+k1,l, cs);
        // here the ijk range is be more special according to the sphere radius rp
        imax=1+(int)((dd0(1)+rp-minv(1))*nofpointsi/max_min(1)+0.5);
        imin=-1+(int)((dd0(1)-rp-minv(1))*nofpointsi/max_min(1)+0.5);
        jmax=1+(int)((dd0(2)+rp-minv(2))*nofpointsj/max_min(2)+0.5);
        jmin=-1+(int)((dd0(2)-rp-minv(2))*nofpointsj/max_min(2)+0.5);
        kmax=1+(int)((dd0(3)+rp-minv(3))*nofpointsk/max_min(3)+0.5);
        kmin=-1+(int)((dd0(3)-rp-minv(3))*nofpointsk/max_min(3)+0.5);
        if (imin<1){imin=1;}if (jmin<1){jmin=1;}if (kmin<1){kmin=1;}if (imax>nofpointsi){imax=nofpointsi;}if (jmax>nofpointsj){jmax=nofpointsj;}if (kmax>nofpointsk){kmax=nofpointsk;}
        for (i=imin;i<=imax;++i){for (j=jmin;j<=jmax;++j){for (k=kmin;k<=kmax;++k){
        // set position vector
        rijk=minv; rijk(1)+=(2*i-1)*max_min(1)/nofpointsi/2;rijk(2)+=(2*j-1)*max_min(2)/nofpointsj/2;rijk(3)+=(2*k-1)*max_min(3)/nofpointsk/2;
        dd=p_inverse*rijk;
        if((dd(1)>0)&&(dd(1)<1)&&(dd(2)>0)&&(dd(2)<1)&&(dd(3)>0)&&(dd(3)<1))
        {dd=dd0-rijk;
        if(Norm(dd)<rp)ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+= copysign(1.6110481,radius);
        // this is the chargedensity of a homogeneous sphere with 1 electron/(4pi a0^3/3) with a0=0.529177 A
        }
        }}}
        }}}
        }}}
  } }
  else
  {jjjpar ionpar(cs.x[l],cs.y[l],cs.z[l],cs.sipffilenames[l],1);
   int ndd;
   // here we should introduce another loop to go around +-1 around the primitive
   // magnetic unit cell so that we see also atoms at the borders in the density map:
   int i0,j0,k0;
   for(i0=-1;i0<=1;++i0){for(j0=-1;j0<=1;++j0){for(k0=-1;k0<=1;++k0){
   for (i1=1;i1<=nofa;++i1){for(j1=1;j1<=nofb;++j1){for(k1=1;k1<=nofc;++k1){
   dd0=pos(i0*nofa+i1,j0*nofb+j1,k0*nofc+k1,l, cs);
 
   Vector moments(1,nofcomponents);
   density cd(gp.title,6,6);
//   Vector momSx(1,49),momLx(1,49),momSy(1,49),momLy(1,49),momSz(1,49),momLz(1,49);
   double QR; // old: QR=hkl(1)*dd0(1)/cs.abc(1)+hkl(2)*dd0(2)/cs.abc(2)+hkl(3)*dd0(3)/cs.abc(3);
   QR=(hkl*abc_in_ijk_Inverse)*dd0;
   QR*=2*PI;int i1r=i1,j1r=j1,k1r=k1;

                for(ndd=1;ndd<=densityev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i1r,j1r,k1r,l)(ndd)+gp.spins_wave_amplitude*(cos(-phase+QR)*densityev_real.moment(i1r,j1r,k1r,l)(ndd)+sin(phase-QR)*densityev_imag.moment(i1r,j1r,k1r,l)(ndd));}
    cd.moments_init(moments);
//   if(strncmp(gp.title+14,"momdensity",10)==0){if(nofcomponents>=3*49)
//                                            {int i1i;for(i1i=1;i1i<=49;++i1i){momSx(i1i)=moments(i1i);momSy(i1i)=moments(i1i+49);momSz(i1i)=moments(i1i+2*49);momLx(i1i)=moments(i1i+3*49);momLy(i1i)=moments(i1i+4*49);momLz(i1i)=moments(i1i+5*49);}}
//                                            else
//                                            {int i1i;for(i1i=1;i1i<=49;++i1i){momSx(i1i)=moments(i1i);momLx(i1i)=moments(i1i+49);}}
//                                           }
//   else if (strncmp(gp.title+14,"orbmomdensity",10)==0&&nofcomponents>=3*49)
//                                           {int i1i;for(i1i=1;i1i<=49;++i1i){momLx(i1i)=moments(i1i);momLy(i1i)=moments(i1i+49);momLz(i1i)=moments(i1i+2*49);}}
//   else if (strncmp(gp.title+14,"spindensity",10)==0&&nofcomponents>=3*49)
//                                           {int i1i;for(i1i=1;i1i<=49;++i1i){momSx(i1i)=moments(i1i);momSy(i1i)=moments(i1i+49);momSz(i1i)=moments(i1i+2*49);}}
//   else if(strncmp(gp.title+14,"currdensity",10)==0){
//                                            int i1i;for(i1i=1;i1i<=49;++i1i){momLx(i1i)=moments(i1i);momLy(i1i)=moments(i1i+49);momLz(i1i)=moments(i1i+2*49);}
//                                                 }

          // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(densityev_real*cos(-phase) + densityev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
        // here the ijk range should be more special according to the maximum sphere radius 3A - to get speed up!!!!
        // dd0 is the center of the atom ...
        radius=3.0;// only pixels nearer maxR (A) to the center of an atom will be considered
        imax=1+(int)((dd0(1)+radius-minv(1))*nofpointsi/max_min(1)+0.5);
        imin=-1+(int)((dd0(1)-radius-minv(1))*nofpointsi/max_min(1)+0.5);
        jmax=1+(int)((dd0(2)+radius-minv(2))*nofpointsj/max_min(2)+0.5);
        jmin=-1+(int)((dd0(2)-radius-minv(2))*nofpointsj/max_min(2)+0.5);
        kmax=1+(int)((dd0(3)+radius-minv(3))*nofpointsk/max_min(3)+0.5);
        kmin=-1+(int)((dd0(3)-radius-minv(3))*nofpointsk/max_min(3)+0.5);
        if (imin<1){imin=1;}if (jmin<1){jmin=1;}if (kmin<1){kmin=1;}if (imax>nofpointsi){imax=nofpointsi;}if (jmax>nofpointsj){jmax=nofpointsj;}if (kmax>nofpointsk){kmax=nofpointsk;}
        for (i=imin;i<=imax;++i){for (j=jmin;j<=jmax;++j){for (k=kmin;k<=kmax;++k){
        // set position vector
        rijk=minv; rijk(1)+=(2*i-1)*max_min(1)/nofpointsi/2;rijk(2)+=(2*j-1)*max_min(2)/nofpointsj/2;rijk(3)+=(2*k-1)*max_min(3)/nofpointsk/2;
        // we should check here if rijk is in primitive unitcell otherwise take next rijk
        dd=p_inverse*rijk;
        if((dd(1)>0)&&(dd(1)<1)&&(dd(2)>0)&&(dd(2)<1)&&(dd(3)>0)&&(dd(3)<1))
        {
        dd=dd0-rijk;
        // get theta phi R from dd
    double R,Rxy,theta,fi;
    R=Norm(dd);
    if(R<radius){ // do not consider any pixels further away than maxR
    if(ionpar.orientation==abc_yzx){// mind abc||yzx in module cfield
     //dx=R*sin(theta)*sin(fi);dy=R*cos(theta);dz=R*sin(theta)*cos(fi);
     theta=acos(dd(2)/R);Rxy=sqrt(dd(1)*dd(1)+dd(3)*dd(3));if(Rxy>SMALL){fi=acos(dd(3)/Rxy);}else{fi=0;}
                         if (dd(1)<0)fi=-fi;
                              }
     else
                              {// mind abc||xyz in other cases ...
     //dx=R*sin(theta)*cos(fi);dy=R*sin(theta)*sin(fi);dz=R*cos(theta);
     theta=acos(dd(3)/R);Rxy=sqrt(dd(1)*dd(1)+dd(2)*dd(2));if(Rxy>SMALL){fi=acos(dd(1)/Rxy);}else{fi=0;}
                         if (dd(2)<0)fi=-fi;
                              }

    ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=cd.denscalc(theta,fi,R,moments,ionpar,T,gjmbHxc,Hext);
//    if(strncmp(gp.title+14,"spindensity",10)==0)
    // here we calculate the spindensity of ion  (negative sign, because rocalc does give positive values)
//    {if(nofcomponents>=3*49)
//     {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=Norm(ionpar.spindensity_calc(theta,fi,R,momSx,momSy,momSz));}
//     else
//     {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=ionpar.spindensity_calc(theta,fi,R,moments);}
//    }

//    if(strncmp(gp.title+14,"orbmomdensity",10)==0)
    // here we calculate the spindensity of ion  (negative sign, because rocalc does give positive values)
//    {if(nofcomponents>=3*49)
//     {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=Norm(ionpar.orbmomdensity_calc(theta,fi,R,momLx,momLy,momLz));}
//     else
//     {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=ionpar.orbmomdensity_calc(theta,fi,R,moments);}
//    }
//    if(strncmp(gp.title+14,"momdensity",10)==0)
    // here we calculate the spindensity of ion  (negative sign, because rocalc does give positive values)
//    {if(nofcomponents>=3*49)
//     {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=Norm(ionpar.spindensity_calc(theta,fi,R,momSx,momSy,momSz)+ionpar.orbmomdensity_calc(theta,fi,R,momLx,momLy,momLz));}
//     else
//     {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=ionpar.spindensity_calc(theta,fi,R,momSx)
//                                                  +ionpar.orbmomdensity_calc(theta,fi,R,momLx);
//    }}
//    if(strncmp(gp.title,"abs value  of currdensity",25)==0)
    // here we calculate the currdensity of ion
//    {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=Norm(ionpar.currdensity_calc(theta,fi,R,momLx,momLy,momLz));
//    }
//    if(strncmp(gp.title,"projection of currdensity",25)==0)
    // here we calculate the currdensity of ion
//    {Vector pr(1,3); extract(gp.title,"i",pr(1));extract(gp.title,"j",pr(2));extract(gp.title,"k",pr(3));
//     ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=pr*ionpar.currdensity_calc(theta,fi,R,momLx,momLy,momLz);
//    }

//    if(strncmp(gp.title,"chargedensity",10)==0)
    // here we calculate the chargedensity of ion  (negative sign, because rocalc does give positive values)
//    {ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]-=ionpar.rocalc(theta,fi,R,moments);}
   // printf("%g %g %g %g\n",R,theta,fi,ro);
        }
   } // end if rijk is in primitive magnetic unit cell
   }}}
   }}}
   }}}
  }
  }

  // here starts printout density loop
  fprintf(fout,"#density map \n");
  fprintf(fout,"#ri[A] rj[A] rk[A] density[|e|/A^3] (or [mb/A^3]) (or  milliAmpere/A^2)\n");
  for (i=1;i<=nofpointsi;++i){for (j=1;j<=nofpointsj;++j){for (k=1;k<=nofpointsk;++k){
  // set position vector
  rijk=minv; rijk(1)+=(2*i-1)*max_min(1)/nofpointsi/2;rijk(2)+=(2*j-1)*max_min(2)/nofpointsj/2;rijk(3)+=(2*k-1)*max_min(3)/nofpointsk/2;
  // print out density
  fprintf(fout,"%10.7f %10.7f %10.7f %g\n",rijk(1),rijk(2),rijk(3),ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]);
                           }}}
 delete []ro;
}

//***********************************************************************************************************************************

// output for fullprof studio
void spincf::fst(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, spincf & magmom) //print std file to stream
{int i,j,k,l,ctr=1;


  Vector maxv(1,3),minv(1,3),dd(1,3),max_min(1,3);
  Vector xyz(1,3),dd0(1,3),ijkmax(1,3),ijkmin(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);
  calc_minmax(minv,maxv,ijkmin,ijkmax,p,abc);


  max_min=maxv-minv;


fprintf(fout,"!   FILE for FullProf Studio: generated automatically by McPhase\n");
fprintf(fout,"!Title: %s \n",text);
fprintf(fout,"SPACEG P 1           \n");
fprintf(fout,"CELL     %g    %g    %g  %g %g %g   DISPLAY MULTIPLE\n",myround(max_min(1)),myround(max_min(2)),myround(max_min(3)),myround(abc(4)),myround(abc(5)),myround(abc(6)));
fprintf(fout,"BOX   -0.15  1.15   -0.15  1.15    -0.15  1.15 \n");

  // plot atoms in region xmin to xmax (quader)
int i1,j1,k1;
for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){

   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);

      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if(dd(1)<=maxv(1)+0.0001&&dd(1)>=minv(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=maxv(2)+0.0001&&dd(2)>=minv(2)-0.0001&&
            dd(3)<=maxv(3)+0.0001&&dd(3)>=minv(3)-0.0001)
            {dd(1)/=max_min(1);dd(2)/=max_min(2);dd(3)/=max_min(3);

fprintf(fout,"ATOM DY%i    RE       %g       %g       %g        \n",ctr,myround(dd(1)),myround(dd(2)),myround(dd(3)));

	     ++ctr;

	     }
	  }
       }}}
 }}}

fprintf(fout," \n");
fprintf(fout,"{\n");
fprintf(fout,"LATTICE P\n");
fprintf(fout,"K     0.00000   0.00000   0.00000\n");
fprintf(fout,"SYMM  x,y,z\n");
fprintf(fout,"MSYM  u,v,w,0.0\n");

// plot moments in region xmin to xmax (quader)
for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){

   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);

      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
          dd+=dd0;
	    if(dd(1)<=maxv(1)+0.0001&&dd(1)>=minv(1)-0.0001&&   //if atom is in big unit cell
            dd(2)<=maxv(2)+0.0001&&dd(2)>=minv(2)-0.0001&&
            dd(3)<=maxv(3)+0.0001&&dd(3)>=minv(3)-0.0001)
            {dd(1)/=max_min(1);dd(2)/=max_min(2);dd(3)/=max_min(3);
//             i1true=1;j1true=1;k1true=1;

//	    a=xy(dd,orientation, minv, maxv,bbwidth,bbheight);
            xyz=magmom.moment(i,j,k,l);
fprintf(fout,"MATOM DY%i    DY      %g       %g       %g   GROUP\n",ctr,myround(dd(1)),myround(dd(2)),myround(dd(3)));
fprintf(fout,"SKP           1  1  %g       %g       %g       0.00000  0.00000  0.00000    0.00000\n",myround(xyz(1)),myround(xyz(2)),myround(xyz(3)));
	     ++ctr;

	     }
	  }
       }}}
 }}}
//}}}
fprintf(fout,"}\n");
}


void spincf::fstprim(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, spincf & magmom) //print std file to stream
{int i,j,k,l,ctr=1;

double alpha,beta,gamma;
  // determine max(1,2,3) min(1,2,3) (vector in Angstroem describing a quader) for viewing magnetic unit cell
  Vector ddd(1,8),xyz(1,3),xyz0(1,3),dd0(1,3),dd(1,3);

  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,abc,r);

gamma=180/3.1415926*acos(p.Column(1)*p.Column(2)/Norm(p.Column(1))/Norm(p.Column(2)));
beta=180/3.1415926*acos(p.Column(1)*p.Column(3)/Norm(p.Column(1))/Norm(p.Column(3)));
alpha=180/3.1415926*acos(p.Column(2)*p.Column(3)/Norm(p.Column(2))/Norm(p.Column(3)));


fprintf(fout,"!   FILE for FullProf Studio: generated automatically by McPhase\n");
fprintf(fout,"!Title: %s \n",text);
fprintf(fout,"SPACEG P 1           \n");
fprintf(fout,"CELL     %g    %g    %g  %g %g %g   DISPLAY MULTIPLE\n",Norm(p.Column(1)),Norm(p.Column(2)),Norm(p.Column(3)),alpha,beta,gamma);
fprintf(fout,"BOX   -0.15  1.15   -0.15  1.15    -0.15  1.15 \n");

  // plot atoms
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd0=p.Inverse()*dd;
fprintf(fout,"ATOM DY%i    RE       %g       %g       %g        \n",ctr,myround(dd0(1)),myround(dd0(2)),myround(dd0(3)));
	     ++ctr;

	     }
	  }
       }}

fprintf(fout," \n");
fprintf(fout,"{\n");
fprintf(fout,"LATTICE P\n");
fprintf(fout,"K     0.00000   0.00000   0.00000\n");
fprintf(fout,"SYMM  x,y,z\n");
fprintf(fout,"MSYM  u,v,w,0.0\n");

// plot moments
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, abc, r,x,y,z);
         dd0=p.Inverse()*dd;
          xyz=magmom.moment(i,j,k,l);
          xyz0=p.Inverse()*xyz; xyz0(1)*=Norm(p.Column(1));xyz0(2)*=Norm(p.Column(2));xyz0(3)*=Norm(p.Column(3));

fprintf(fout,"MATOM DY%i    DY      %g       %g       %g   GROUP\n",ctr,myround(dd0(1)),myround(dd0(2)),myround(dd0(3)));
fprintf(fout,"SKP           1  1  %g       %g       %g       0.00000  0.00000  0.00000    0.00000\n",myround(xyz0(1)),myround(xyz0(2)),myround(xyz0(3)));
	     ++ctr;

	     }
	  }
       }}
fprintf(fout,"}\n");
}


//-----------------------------------------------------------------------
//  numeric output of spinconfiguration to file
void spincf::print(FILE * fout) //print spinconfiguration to stream
{print(fout,nofcomponents);}

void spincf::print(FILE * fout,int nofcomp) //print spinconfiguration to stream
{int i,j,k,l;
 if(nofcomp>nofcomponents){fprintf(stderr,"Error spincf::print: nofcomp=%i > nofcomponents=%i\n",nofcomp,nofcomponents);exit(1);}
 if(nofcomp<1){fprintf(stderr,"Error spincf::print: nofcomp=%i <1 \n",nofcomp);exit(1);}

 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)if((l-1)%nofcomponents<nofcomp)
   {for (i=1;i<=nofa;++i)
      {fprintf(fout," %4.4f",myround(1e-5,mom[in(i,j,k)](l)));
       }
    fprintf(fout,"\n");
    }
   }
 fprintf(fout,"\n"); //new line to separate ab planes
 }
// fprintf(fout,"\n"); //new line to end spinconfiguration - removed aug 07
}
//-----------------------------------------------------------------------
//  numeric output of spinconfiguration to file with comments and only large spin-coponents
// components in interval [min,max] are shown

void spincf::print_commented(FILE * fout,const char * string,int min, int max, int maxnofpars, double & absvallimit)
{int i,j,k,l;div_t result;
 if (maxnofpars<1){fprintf(stderr,"Error spincf print_commented: maxnofpars <1\n");exit(1);}
 else
 {   if(min<1||max>nofcomponents){
  fprintf(stderr,"Error spincf::print_commmented: index interval is [%i,%i] -  must be in range [%i,%i]\n",min,max,1,nofcomponents);exit(1);}

  fprintf(fout,"# Output of  %s[%i-%i], at maximum %i numbers (in primitive unit cell) larger than %4.4f:\n",string,min,max,maxnofpars, absvallimit);
  for (i=1;i<=nofa;++i)
  for (j=1;j<=nofb;++j)
  for (k=1;k<=nofc;++k)
  {Vector ml(1,maxnofpars);ml=0;
   Vector mmom(1,maxnofpars);mmom=0;
   int pout=0;fprintf(fout,"#supercell-index %i %i %i:\n",i,j,k);
      for (l=1;l<=nofcomponents*nofatoms;++l)
      {result=div(l-1,nofcomponents);
       if(fabs(mom[in(i,j,k)](l))>absvallimit&&min<=(result.rem+1)&&(result.rem+1)<=max){if(pout<maxnofpars){ ++pout;
                                                                     mmom(pout)=mom[in(i,j,k)](l);
                                                                     ml(pout)=l;
                                                                  }
                                               else { int ii; Min(mmom%mmom,ii);
                                                      mmom(ii)=mom[in(i,j,k)](l);
                                                      ml(ii)=l;
                                                    }
                                              }
       }
       //myPrintVector(stdout,mmom);
       Sort(ml,mmom);
       //myPrintVector(stdout,mmom);
    for (int p=maxnofpars-pout+1;p<=maxnofpars;++p)
    {result=div((int)ml(p)-1,nofcomponents); 
     fprintf(fout,"#!    atom %i: %s%i=%s_%i_%i_%i_%i_%i=%4.4f\n",result.quot+1,string,result.rem+1,string,i,j,k,result.quot+1,result.rem+1,myround(1e-5,mmom(p)));
    }
      
// fprintf(fout,"\n"); //new line to separate ab planes
}
// fprintf(fout,"\n"); //new line to end spinconfiguration - removed aug 07
 }
}