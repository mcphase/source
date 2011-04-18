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
  Vector dd(1,3),nofabc(1,3);
  Matrix rijk(1,3,1,3);
  dadbdc2ijk(rijk,r,abc);
  nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  for (i=1;i<=3;++i)
  {for(j=1;j<=3;++j) {dd(j)=nofabc(j)*rijk(i,j);// old: dd(j)=nofabc(j)*r(i,j)*abc(i);
                      p(i,j)=dd(j);}
  }
 // pa=p.Column(1);  //primitive magnetic unit cell
 // pb=p.Column(2);
 // pc=p.Column(3);
}

void spincf::calc_prim_mag_unitcell_old(Matrix & p,Vector & abc, Matrix & r)
{ int i,j;
  Vector dd(1,3),nofabc(1,3);
  Matrix rijk(1,3,1,3);
  dadbdc2ijk(rijk,r,abc);
  nofabc(1)=nofa;nofabc(2)=nofb;nofabc(3)=nofc;
  for (i=1;i<=3;++i)
  {for(j=1;j<=3;++j) {dd(j)=nofabc(j)*r(i,j)*abc(i);
                      p(i,j)=dd(j);}
  }
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

void calc_minmax_scale_relabc(Vector & minv,Vector & maxv,Vector & ijkmin,Vector & ijkmax,Matrix & p,Vector & abc,double scale_view_1,double scale_view_2,double scale_view_3)
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
  for(i=1;i<=3;++i){minv(i)/=abc(i);maxv(i)/=abc(i);}
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


void spincf::eps3d(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z,int orientation, Vector & gJ)
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
      {c=magmom(i,j,k,l,gJ(l));
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
            {c=magmom(i,j,k,l,gJ(l));
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

//output for javaview
void spincf::jvx_cd(FILE * fout,char * text,cryststruct & cs,graphic_parameters & gp,
                    double phase,spincf & savev_real,spincf & savev_imag,Vector & hkl,double & T, Vector &  gjmbH)
{ int i,j,k,l,ctr=0;int i1,j1,k1;
 // some checks
 if(nofatoms!=savev_real.nofatoms||nofa!=savev_real.na()||nofb!=savev_real.nb()||nofc!=savev_real.nc()||
    nofatoms!=savev_imag.nofatoms||nofa!=savev_imag.na()||nofb!=savev_imag.nb()||nofc!=savev_imag.nc()||
    nofcomponents<savev_real.nofcomponents||savev_real.nofcomponents!=savev_imag.nofcomponents)
    {fprintf(stderr,"Error creating jvx movie files: eigenvector read from .qev file does not match dimension of spins structure read from sps file\n");exit(1);}

  Vector maxv(1,3),minv(1,3),ijkmax(1,3),ijkmin(1,3),max_min(1,3),dd(1,3),dd0(1,3),c(1,3),xyz(1,3);
  Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
  Matrix abc_in_ijk_Inverse(1,3,1,3); abc_in_ijk_Inverse=abc_in_ijk.Inverse();
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell_old(p,cs.abc,cs.r);
  calc_minmax_scale_relabc(minv,maxv,ijkmin,ijkmax,p,cs.abc,gp.scale_view_1,gp.scale_view_2,gp.scale_view_3);
   if(gp.showprim==1){ijkmin(1)=1;ijkmin(2)=1;ijkmin(3)=1;ijkmax(1)=-2+(int)(gp.scale_view_1);ijkmax(2)=-2+(int)(gp.scale_view_2);ijkmax(3)=-2+(int)(gp.scale_view_3);} // show only primitive magnetic unit cell
  max_min=maxv-minv;
  calc_prim_mag_unitcell(p,cs.abc,cs.r);
  

fprintf(fout,"<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"no\"?>\n");
fprintf(fout,"<jvx-model>\n");
fprintf(fout,"  <title>%s</title>\n",text);
fprintf(fout,"  <geometries>\n");

if(gp.show_abc_unitcell>0)
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
if(gp.show_primitive_crystal_unitcell>0)
 {
 // plot frame around primitive crystallographic unit cell
fprintf(fout,"    <geometry name=\"primitive crystallographic unit cell\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
dd=0;           fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\"r1\">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(2)/(double)nofb;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd-=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\"r2\">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd =p.Column(3)/(double)nofc;fprintf(fout,"          <p name=\"r3\">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(2)/(double)nofb;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd-=p.Column(1)/(double)nofa;fprintf(fout,"          <p name=\" \">  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
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
if(gp.show_magnetic_unitcell>0)
 {
  // plot frame around primitive magnetic unit cell
fprintf(fout,"    <geometry name=\"magnetic unit cell\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n");
fprintf(fout,"        <points>\n");
dd=0;           fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(2);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd-=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd =p.Column(3);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd+=p.Column(2);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
dd-=p.Column(1);fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
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
  // plot atoms in region xmin to xmax (quader)
fprintf(fout,"    <geometry name=\"ions\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"show\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
  for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l,cs);
         dd+=dd0;
            if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
            (gp.showprim==1&&gp.scale_view_1>(double)(i1*nofa+i)/nofa&&gp.scale_view_2>(double)(j1*nofb+j)/nofb&&gp.scale_view_3>(double)(k1*nofc+k)/nofc))
            {
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
     ++ctr;
	     }
	  }
       }}}
  }}}
fprintf(fout,"          <thickness>3.0</thickness>\n");
fprintf(fout,"        </points>\n");
fprintf(fout,"        <colors type=\"rgb\">\n");
  for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, cs);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
             (gp.showprim==1&&gp.scale_view_1>(double)(i1*nofa+i)/nofa&&gp.scale_view_2>(double)(j1*nofb+j)/nofb&&gp.scale_view_3>(double)(k1*nofc+k)/nofc))
            {
fprintf(fout,"          <c>  %i       %i       %i </c>\n",(int)(255*gp.show_atoms),(int)(gp.show_atoms*((l*97)%256)),0);
	     }
	  }
       }}}
  }}}
fprintf(fout,"        </colors>\n");
fprintf(fout,"      </pointSet>\n");
fprintf(fout,"    </geometry>\n");




if(gp.spins_scale_moment>0){
fprintf(fout,"    <geometry name=\"magnetic moments\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"show\">\n");
fprintf(fout,"        <points>\n");
 ctr=0;
 for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
 for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
 for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, cs);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
                   (gp.showprim==1&&gp.scale_view_1>(double)(i1*nofa+i)/nofa&&gp.scale_view_2>(double)(j1*nofb+j)/nofb&&gp.scale_view_3>(double)(k1*nofc+k)/nofc))
            {double QR; // old: QR=hkl(1)*dd(1)/cs.abc(1)+hkl(2)*dd(2)/cs.abc(2)+hkl(3)*dd(3)/cs.abc(3);
             QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI;
             xyz=magmom(i,j,k,l,cs.gJ[l])+gp.spins_wave_amplitude*(cos(-phase+QR)*savev_real.magmom(i,j,k,l,cs.gJ[l])+sin(phase-QR)*savev_imag.magmom(i,j,k,l,cs.gJ[l]));
              //printf("gJ=%g magmom=%g %g %g %g %g %g %g %g %g\n",cs.gJ[l],mom[in(i,j,k)](1),mom[in(i,j,k)](2),mom[in(i,j,k)](3),mom[in(i,j,k)](4),mom[in(i,j,k)](5),mom[in(i,j,k)](6),xyz(1),xyz(2),xyz(3));
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
//fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)-xyz(1)*gp.spins_scale_moment,dd(2)-xyz(2)*gp.spins_scale_moment,dd(3)-xyz(3)*gp.spins_scale_moment);
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)+xyz(1)*gp.spins_scale_moment,dd(2)+xyz(2)*gp.spins_scale_moment,dd(3)+xyz(3)*gp.spins_scale_moment);
	     ++ctr;

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
fprintf(fout,"        <color type=\"rgb\">0 150 0</color>\n");
fprintf(fout,"          <thickness>4.0</thickness>\n");
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
                     }

if(gp.spins_show_static_moment_direction>0)
 {
// plot a line along static magnetic moments for comparison
fprintf(fout,"    <geometry name=\"static magnetic moments\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n");
fprintf(fout,"        <points>\n");
 ctr=0;
 for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
 for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
 for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, cs);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
                     (gp.showprim==1&&gp.scale_view_1>(double)(i1*nofa+i)/nofa&&gp.scale_view_2>(double)(j1*nofb+j)/nofb&&gp.scale_view_3>(double)(k1*nofc+k)/nofc))
            {double QR; // old: QR=hkl(1)*dd(1)/cs.abc(1)+hkl(2)*dd(2)/cs.abc(2)+hkl(3)*dd(3)/cs.abc(3);
             QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI;
                          xyz=magmom(i,j,k,l,cs.gJ[l]);
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1),dd(2),dd(3));
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)+xyz(1)*gp.spins_scale_moment,dd(2)+xyz(2)*gp.spins_scale_moment,dd(3)+xyz(3)*gp.spins_scale_moment);
	     ++ctr;

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
fprintf(fout,"          <thickness>%g</thickness>\n",gp.spins_scale_moment);
fprintf(fout,"        </lines>\n");
fprintf(fout,"      </lineSet>\n");
fprintf(fout,"    </geometry>\n");
}
if(gp.spins_show_ellipses>0)
 {
// plot an ellipse along path of moment
fprintf(fout,"    <geometry name=\"ellipses\">\n");
fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"hide\">\n");
fprintf(fout,"        <points>\n");
 ctr=0;
 for (i1=int(ijkmin(1)-1.0);i1<=int(ijkmax(1)+1);++i1){
 for (j1=int(ijkmin(2)-1.0);j1<=int(ijkmax(2)+1);++j1){
 for (k1=int(ijkmin(3)-1.0);k1<=int(ijkmax(3)+1);++k1){
   dd0=p.Column(1)*(double)(i1)+p.Column(2)*(double)(j1)+p.Column(3)*(double)(k1);
      for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, cs);
          dd+=dd0;if(check_atom_in_big_unitcell(dd,maxv,minv,abc_in_ijk_Inverse)||
                    (gp.showprim==1&&gp.scale_view_1>(double)(i1*nofa+i)/nofa&&gp.scale_view_2>(double)(j1*nofb+j)/nofb&&gp.scale_view_3>(double)(k1*nofc+k)/nofc))
            {double QR; // old: QR=hkl(1)*dd(1)/cs.abc(1)+hkl(2)*dd(2)/cs.abc(2)+hkl(3)*dd(3)/cs.abc(3);
             QR=(hkl*abc_in_ijk_Inverse)*dd;
             QR*=2*PI;
             int phi;
             for(phi=0;phi<=16;phi++)
             {
             xyz=magmom(i,j,k,l,cs.gJ[l])+gp.spins_wave_amplitude*(cos(-(double)phi*2*3.1415/16+QR)*savev_real.magmom(i,j,k,l,cs.gJ[l])+sin((double)phi*2*3.1415/16-QR)*savev_imag.magmom(i,j,k,l,cs.gJ[l]));
              // <Jalpha>(i)=<Jalpha>0(i)+gp.spins_wave_amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*gp.spins_wave_amplitude; // Q ri not considered for test !!!
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dd(1)+xyz(1)*gp.spins_scale_moment,dd(2)+xyz(2)*gp.spins_scale_moment,dd(3)+xyz(3)*gp.spins_scale_moment);
	     }++ctr;

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
}

if(gp.scale_density_vectors>0)
{int ii;
  double dtheta=0.2; //stepwidth to step surface
  double dfi=0.2;

for(l=1;l<=nofatoms;++l)
 {fprintf(fout,"    <geometry name=\"density vectors in primitive magnetic unit cell - atom %i\">\n",l);
  fprintf(fout,"      <pointSet dim=\"3\" point=\"hide\" color=\"show\">\n");
  fprintf(fout,"<points >\n");
  double radius=0;double dx,dy,dz,R,fi,theta;
  int ctr=0;
  extract(cs.cffilenames[l],"radius",radius);
  if(radius==0) // this is a trick: if radius is given as cffilename then a sphere with this is radius is generated (pointcharge)
  {jjjpar ionpar(cs.x[l],cs.y[l],cs.z[l],cs.cffilenames[l]);
   density cd(gp.title,dtheta,dfi);int ndd;
   for (i=1;i<=1+(nofa-1)*gp.scale_view_1;++i){for(j=1;j<=1+(nofb-1)*gp.scale_view_2;++j){for(k=1;k<=1+(nofc-1)*gp.scale_view_2;++k){
   dd=pos(i,j,k,l, cs);
   Vector moments(1,nofcomponents);//printf("nofcomp=%i\n",nofcomponents);
   double QR; // old: QR=hkl(1)*dd(1)/cs.abc(1)+hkl(2)*dd(2)/cs.abc(2)+hkl(3)*dd(3)/cs.abc(3);
   QR=(hkl*abc_in_ijk_Inverse)*dd;
   QR*=2*PI;
                for(ndd=1;ndd<=savev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i,j,k,l)(ndd)+gp.spins_wave_amplitude*(cos(-phase+QR)*savev_real.moment(i,j,k,l)(ndd)+sin(phase-QR)*savev_imag.moment(i,j,k,l)(ndd));}
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
 // here we calculate the chargedensity of ion
   cd.calc_cd_surface(moments,ionpar,gp.threshhold,T,gjmbH);
   for(ii=1;ii<=cd.nofpoints();++ii)
     {R=cd.rtf(ii)(1);theta=cd.rtf(ii)(2);fi=cd.rtf(ii)(3);
     if(ionpar.module_type==2){// mind abc||yzx in module cfield
     dx=R*sin(theta)*sin(fi)+dd(1);dy=R*cos(theta)+dd(2);dz=R*sin(theta)*cos(fi)+dd(3);
                              }
     else
                              {// mind abc||xyz in other cases ...
     dx=R*sin(theta)*cos(fi)+dd(1);dy=R*sin(theta)*sin(fi)+dd(2);dz=R*cos(theta)+dd(3);
                              }
     
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dx,dy,dz);
fprintf(fout,"          <p>  %g       %g       %g </p>\n",dx+cd.rtf(ii)(4)*gp.scale_density_vectors,dy+cd.rtf(ii)(5)*gp.scale_density_vectors,dz+cd.rtf(ii)(6)*gp.scale_density_vectors);
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

if(gp.show_density>0)
{int ii,tt,ff;
  double dtheta=gp.density_dtheta; //stepwidth to step surface
  double dfi=gp.density_dfi;

for(l=1;l<=nofatoms;++l)
 {fprintf(fout,"    <geometry name=\"densities in primitive magnetic unit cell - atom %i\">\n",l);
  fprintf(fout,"<pointSet color=\"hide\" point=\"show\" dim=\"1\">\n");
  fprintf(fout,"<points >\n");
  double radius=0;double dx,dy,dz,R,fi,theta;

  extract(cs.cffilenames[l],"radius",radius);
  if(radius!=0) // this is a trick: if radius is given as cffilename then a sphere with this is radius is generated (pointcharge)
  {   if(gp.show_pointcharges>0)
      {  double rp=abs(radius);
        for (i=1;i<=(1+(nofa-1)*gp.scale_view_1);++i){for(j=1;j<=(1+(nofb-1)*gp.scale_view_2);++j){for(k=1;k<=(1+(nofc-1)*gp.scale_view_3);++k){
        dd=pos(i,j,k,l, cs);
        for(tt=0;tt<=3.1415/dtheta;++tt){for(ff=0;ff<=2*3.1415/dfi;++ff){
             theta=(double)tt*dtheta;fi=(double)ff*dfi;
             dx=rp*sin(theta)*cos(fi)+dd(1);dy=rp*sin(theta)*sin(fi)+dd(2);dz=rp*cos(theta)+dd(3);
             fprintf(fout,"<p>%4g %4g %4g</p>\n",dx,dy,dz);
             if(tt==0){ff=(int)(2*3.1415/dfi+1);}
             }}
        }}}
  } }
  else
  {jjjpar ionpar(cs.x[l],cs.y[l],cs.z[l],cs.cffilenames[l]);
   density cd(gp.title,dtheta,dfi);int ndd;
   for (i=1;i<=1+(nofa-1)*gp.scale_view_1;++i){for(j=1;j<=1+(nofb-1)*gp.scale_view_2;++j){for(k=1;k<=1+(nofc-1)*gp.scale_view_2;++k){
   dd=pos(i,j,k,l, cs);
   Vector moments(1,nofcomponents);
   double QR; // old: QR=hkl(1)*dd(1)/cs.abc(1)+hkl(2)*dd(2)/cs.abc(2)+hkl(3)*dd(3)/cs.abc(3);
   QR=(hkl*abc_in_ijk_Inverse)*dd;
   QR*=2*PI;
                for(ndd=1;ndd<=savev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i,j,k,l)(ndd)+gp.spins_wave_amplitude*(cos(-phase+QR)*savev_real.moment(i,j,k,l)(ndd)+sin(phase-QR)*savev_imag.moment(i,j,k,l)(ndd));}
              // <Jalpha>(i)=<Jalpha>0(i)+amplitude * real( exp(-i omega t+ Q ri) <ev_alpha>(i) )
              // omega t= phase
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
 // here we calculate the chargedensity of ion
   cd.calc_cd_surface(moments,ionpar,gp.threshhold,T,  gjmbH);
   for(ii=1;ii<=cd.nofpoints();++ii)
     {R=cd.rtf(ii)(1);theta=cd.rtf(ii)(2);fi=cd.rtf(ii)(3);
     if(ionpar.module_type==2){// mind abc||yzx in module cfield
     dx=R*sin(theta)*sin(fi)+dd(1);dy=R*cos(theta)+dd(2);dz=R*sin(theta)*cos(fi)+dd(3);
                              }
     else
                              {// mind abc||xyz in other cases ...
     dx=R*sin(theta)*cos(fi)+dd(1);dy=R*sin(theta)*sin(fi)+dd(2);dz=R*cos(theta)+dd(3);
                              }
     fprintf(fout,"<p>%4g %4g %4g</p>\n",dx,dy,dz);
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



fprintf(fout,"  </geometries>\n");
fprintf(fout,"</jvx-model>\n");
}
//***********************************************************************************************************************************
// output of chargedensity on grid as ascii file points are equally spaced as specified
// nofpoints*
void spincf::cd(FILE * fout,cryststruct & cs, graphic_parameters & gp,
                spincf & savev_real,spincf & savev_imag,double phase,Vector & hkl,double & T, Vector &  gjmbH)
{// some checks
 if(nofatoms!=savev_real.nofatoms||nofa!=savev_real.na()||nofb!=savev_real.nb()||nofc!=savev_real.nc()||
    nofatoms!=savev_imag.nofatoms||nofa!=savev_imag.na()||nofb!=savev_imag.nb()||nofc!=savev_imag.nc()||
    nofcomponents<savev_real.nofcomponents||savev_real.nofcomponents!=savev_imag.nofcomponents)
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
  extract(cs.cffilenames[l],"radius",radius);
  if(radius!=0) // this is a trick: if radius is given as cffilename then a sphere with this is radius is generated (pointcharge)
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
        if (imin<1)imin=1;if (jmin<1)jmin=1;if (kmin<1)kmin=1;if (imax>nofpointsi)imax=nofpointsi;if (jmax>nofpointsj)jmax=nofpointsj;if (kmax>nofpointsk)kmax=nofpointsk;
        for (i=imin;i<=imax;++i){for (j=jmin;j<=jmax;++j){for (k=kmin;k<=kmax;++k){
        // set position vector
        rijk=minv; rijk(1)+=(2*i-1)*max_min(1)/nofpointsi/2;rijk(2)+=(2*j-1)*max_min(2)/nofpointsj/2;rijk(3)+=(2*k-1)*max_min(3)/nofpointsk/2;
        dd=p_inverse*rijk;
        if((dd(1)>0)&(dd(1)<1)&(dd(2)>0)&(dd(2)<1)&(dd(3)>0)&(dd(3)<1))
        {dd=dd0-rijk;
        if(Norm(dd)<rp)ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+= copysign(1.6110481,radius);
        // this is the chargedensity of a homogeneous sphere with 1 electron/(4pi a0^3/3) with a0=0.529177 A
        }
        }}}
        }}}
        }}}
  } }
  else
  {jjjpar ionpar(cs.x[l],cs.y[l],cs.z[l],cs.cffilenames[l]);
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

                for(ndd=1;ndd<=savev_real.nofcomponents;++ndd)
   {moments(ndd)=moment(i1r,j1r,k1r,l)(ndd)+gp.spins_wave_amplitude*(cos(-phase+QR)*savev_real.moment(i1r,j1r,k1r,l)(ndd)+sin(phase-QR)*savev_imag.moment(i1r,j1r,k1r,l)(ndd));}
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
              //spins=savspins+(savev_real*cos(-phase) + savev_imag*sin(phase))*amplitude; // Q ri not considered for test !!!
        // here the ijk range should be more special according to the maximum sphere radius 3A - to get speed up!!!!
        // dd0 is the center of the atom ...
        radius=3.0;// only pixels nearer maxR (A) to the center of an atom will be considered
        imax=1+(int)((dd0(1)+radius-minv(1))*nofpointsi/max_min(1)+0.5);
        imin=-1+(int)((dd0(1)-radius-minv(1))*nofpointsi/max_min(1)+0.5);
        jmax=1+(int)((dd0(2)+radius-minv(2))*nofpointsj/max_min(2)+0.5);
        jmin=-1+(int)((dd0(2)-radius-minv(2))*nofpointsj/max_min(2)+0.5);
        kmax=1+(int)((dd0(3)+radius-minv(3))*nofpointsk/max_min(3)+0.5);
        kmin=-1+(int)((dd0(3)-radius-minv(3))*nofpointsk/max_min(3)+0.5);
        if (imin<1)imin=1;if (jmin<1)jmin=1;if (kmin<1)kmin=1;if (imax>nofpointsi)imax=nofpointsi;if (jmax>nofpointsj)jmax=nofpointsj;if (kmax>nofpointsk)kmax=nofpointsk;
        for (i=imin;i<=imax;++i){for (j=jmin;j<=jmax;++j){for (k=kmin;k<=kmax;++k){
        // set position vector
        rijk=minv; rijk(1)+=(2*i-1)*max_min(1)/nofpointsi/2;rijk(2)+=(2*j-1)*max_min(2)/nofpointsj/2;rijk(3)+=(2*k-1)*max_min(3)/nofpointsk/2;
        // we should check here if rijk is in primitive unitcell otherwise take next rijk
        dd=p_inverse*rijk;
        if((dd(1)>0)&(dd(1)<1)&(dd(2)>0)&(dd(2)<1)&(dd(3)>0)&(dd(3)<1))
        {
        dd=dd0-rijk;
        // get theta phi R from dd
    double R,Rxy,theta,fi;
    R=Norm(dd);
    if(R<radius){ // do not consider any pixels further away than maxR
    if(ionpar.module_type==2){// mind abc||yzx in module cfield
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

    ro[((i-1)*nofpointsj+(j-1))*nofpointsk+k-1]+=cd.denscalc(theta,fi,R,moments,ionpar,T,gjmbH);
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
void spincf::fst(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ) //print std file to stream
{int i,j,k,l,ctr=1;


  Vector maxv(1,3),minv(1,3),dd(1,3),max_min(1,3);
  Vector xyz(1,3),dd0(1,3),ijkmax(1,3),ijkmin(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell_old(p,abc,r);
  calc_minmax(minv,maxv,ijkmin,ijkmax,p,abc);


  max_min=maxv-minv;


fprintf(fout,"!   FILE for FullProf Studio: generated automatically by McPhase\n");
fprintf(fout,"!Title: %s \n",text);
fprintf(fout,"SPACEG P 1           \n");
fprintf(fout,"CELL     %g    %g    %g  %g %g %g   DISPLAY MULTIPLE\n",max_min(1),max_min(2),max_min(3),abc(4),abc(5),abc(6));
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

fprintf(fout,"ATOM DY%i    RE       %g       %g       %g        \n",ctr,dd(1),dd(2),dd(3));

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
            xyz=magmom(i,j,k,l,gJ(l));
fprintf(fout,"MATOM DY%i    DY      %g       %g       %g   GROUP\n",ctr,dd(1),dd(2),dd(3));
fprintf(fout,"SKP           1  1  %g       %g       %g       0.00000  0.00000  0.00000    0.00000\n",xyz(1),xyz(2),xyz(3));
	     ++ctr;

	     }
	  }
       }}}
 }}}
//}}}
fprintf(fout,"}\n");
}


void spincf::fstprim(FILE * fout,char * text,Vector & abc,Matrix & r,float * x,float *y,float*z, Vector & gJ) //print std file to stream
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
fprintf(fout,"ATOM DY%i    RE       %g       %g       %g        \n",ctr,dd0(1),dd0(2),dd0(3));
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
          xyz=magmom(i,j,k,l,gJ(l));
          xyz0=p.Inverse()*xyz; xyz0(1)*=Norm(p.Column(1));xyz0(2)*=Norm(p.Column(2));xyz0(3)*=Norm(p.Column(3));

fprintf(fout,"MATOM DY%i    DY      %g       %g       %g   GROUP\n",ctr,dd0(1),dd0(2),dd0(3));
fprintf(fout,"SKP           1  1  %g       %g       %g       0.00000  0.00000  0.00000    0.00000\n",xyz0(1),xyz0(2),xyz0(3));
	     ++ctr;

	     }
	  }
       }}
fprintf(fout,"}\n");
}


//-----------------------------------------------------------------------
//  numeric output of spinconfiguration to file
void spincf::print(FILE * fout) //print spinconfiguration to stream
{int i,j,k,l;
 for (k=1;k<=nofc;++k)
 {for (j=1;j<=nofb;++j)
  {for (l=1;l<=nofcomponents*nofatoms;++l)
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

void spincf::printall(FILE * fout,cryststruct & cs) //print spinconfiguration to stream
{ int i,j,k,l,lc,m,maxm;

 // determine primitive magnetic unit cell
  Vector dd(1,3),ddp(1,3);
  Vector xyz(1,3),dd0(1,3),mmm(1,3);
  Matrix p(1,3,1,3);
  calc_prim_mag_unitcell(p,cs.abc,cs.r);
  Matrix abc_in_ijk(1,3,1,3); get_abc_in_ijk(abc_in_ijk,cs.abc);
  Matrix abc_in_ijk_Inverse(1,3,1,3); abc_in_ijk_Inverse=abc_in_ijk.Inverse();


 fprintf(fout,"#!nr1=%i nr2=%i nr3=%i nat=%i atoms in primitive magnetic unit cell:\n",nofa,nofb,nofc,nofatoms*nofa*nofb*nofc);
 fprintf(fout,"#{atom file} da[a] db[b] dc[c] dr1[r1] dr2[r2] dr3[r3]  <Ma> <Mb> <Mc> [mb] <Ja> <Jb> <Jc> <Jd> <Je> ...\n");

   // output atoms and moments in primitive unit cell
  for (i=1;i<=nofa;++i){for (j=1;j<=nofb;++j){for (k=1;k<=nofc;++k){
         for(l=1;l<=nofatoms;++l)
	 {dd=pos(i,j,k,l, cs);
         dd0=p.Inverse()*dd;dd0(1)*=nofa;dd0(2)*=nofb;dd0(3)*=nofc;
         ddp=abc_in_ijk_Inverse*dd;
              fprintf(fout,"{%s} %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f ",
	              cs.cffilenames[l],ddp(1),ddp(2),ddp(3),dd0(1),dd0(2),dd0(3));
             if(cs.gJ[l]!=0)
              {fprintf(fout," %4.4f",myround(1e-5,cs.gJ[l]*mom[in(i,j,k)](1+nofcomponents*(l-1))));
               if(nofcomponents>=2){fprintf(fout," %4.4f",myround(1e-5,cs.gJ[l]*mom[in(i,j,k)](2+nofcomponents*(l-1))));}else{fprintf(fout," %4.4f",0.0);}
               if(nofcomponents>=2){fprintf(fout," %4.4f",myround(1e-5,cs.gJ[l]*mom[in(i,j,k)](3+nofcomponents*(l-1))));}else{fprintf(fout," %4.4f",0.0);}
              }
             else   // if gJ=0 it means we have so print out total moment
              { //load magnetic moment into vector mmm
               if(nofcomponents>6){maxm=6;}else{maxm=nofcomponents;}
                mmm=0;
                for(m=1;m<=maxm;++m){if(m==2||m==4||m==6){mmm((m+1)/2)+=mom[in(i,j,k)](nofcomponents*(l-1)+m);}
                                     else                {mmm((m+1)/2)+=2*mom[in(i,j,k)](nofcomponents*(l-1)+m);}
                                     }

               fprintf(fout," %4.4f %4.4f %4.4f",myround(1e-5,mmm(1)),myround(1e-5,mmm(2)),myround(1e-5,mmm(3)));
              }
             {for (lc=1;lc<=nofcomponents;++lc)
              {fprintf(fout," %4.4f",myround(1e-5,mom[in(i,j,k)](lc+nofcomponents*(l-1))));}
              fprintf(fout,"\n");
	     }

	 }
  }}}

}

