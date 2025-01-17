    // routines for calculation of intensities for program mcdiff



// different output data for user defined columns ...
double setcoloutput(int i,inimcdiff & ini,float & scale, double & ovallt,float & lorentzf,complex <double> & nsf,float & msf2,float & msf2dip,
                   complex <double> & msfx, complex <double> & msfy, complex <double> & msfz,
                   complex <double> & msfdipx, complex <double> &msfdipy, complex <double> &msfdipz,Vector & Qvec,
                   float & d,float & theta,double & Q,float &inuc,float	&Imag,float&Imagdip,Vector & hkl)
{double cosw,crossx,crossy,crossz,Ip,Im,sinw2,R;

 double hklprim;
         // here do some precalculations with formulas common to several options
         switch (i) {case 13: // beyond cases
                     case 14:
                     case 15:
                     case 16:
                     case 23:
                     case 24:crossx=imag(msfy*conj(msfz)-msfz*conj(msfy));
                             crossy=imag(-msfx*conj(msfz)+msfz*conj(msfx));
                             crossz=imag(msfx*conj(msfy)-msfy*conj(msfx));
                            Ip=abs(nsf) * abs(nsf)+msf2 * 3.65 / 4 / PI;
                            Ip-=(crossx*ini.Pxyz(1)+crossy*ini.Pxyz(2)+crossz*ini.Pxyz(3))* 3.65 / 4 / PI;
                            Ip+=sqrt(3.65/4/PI)*real(nsf*conj(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3))+conj(nsf)*(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3)));
                            Im=abs(nsf) * abs(nsf)+msf2 * 3.65 / 4 / PI;
                            Im+=(crossx*ini.Pxyz(1)+crossy*ini.Pxyz(2)+crossz*ini.Pxyz(3))* 3.65 / 4 / PI;
                            Im-=sqrt(3.65/4/PI)*real(nsf*conj(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3))+conj(nsf)*(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3)));
                             break;
                     case 17:  // dipolar cases
                     case 18:
                     case 19:
                     case 20:
                     case 25:
                     case 26: crossx=imag(msfdipy*conj(msfdipz)-msfdipz*conj(msfdipy));
                              crossy=imag(-msfdipx*conj(msfdipz)+msfdipz*conj(msfdipx));
                              crossz=imag(msfdipx*conj(msfdipy)-msfdipy*conj(msfdipx));
                             Ip=abs(nsf) * abs(nsf)+msf2dip * 3.65 / 4 / PI;
                             Ip-=(crossx*ini.Pxyz(1)+crossy*ini.Pxyz(2)+crossz*ini.Pxyz(3))* 3.65 / 4 / PI;
                             Ip+=sqrt(3.65/4/PI)*real(nsf*conj(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3))+conj(nsf)*(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3)));
                             Im=abs(nsf) * abs(nsf)+msf2dip * 3.65 / 4 / PI;
                             Im+=(crossx*ini.Pxyz(1)+crossy*ini.Pxyz(2)+crossz*ini.Pxyz(3))* 3.65 / 4 / PI;
                             Im-=sqrt(3.65/4/PI)*real(nsf*conj(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3))+conj(nsf)*(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3)));
                              break;
                             
                     default: break;
                    }

         switch (i) {
case 0:  return lorentzf;break;//   "LF          ",
case 1:  return abs(nsf);break;//    "|NSF|[b]    ",
case 2:  return real(nsf);break;//    "Re(NSF)[b]  ",
case 3:  return imag(nsf);break;//    "Im(NSF)[b]  ",
case 4:  return sqrt(msf2+1e-100);break;//    "|MSF|       ",
case 5:  return abs(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3));break;//    "|MSF.P|     ",
case 6:  return real(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3));break;//    "Re(MSF.P)   ",
case 7:  return imag(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3));break;//   "Im(MSF.P)   ",
case 8:  return sqrt(msf2dip+1e-100);break;//    "|MSFdip|    ",
case 9:  return abs(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3));break;//    "|MSFdip.P|  ",
case 10: return real(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3));break;//    "Re(MSFdip.P)",
case 11: return imag(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3));break;//    "Im(MSFdip.P)"
case 12: cosw=(ini.Pxyz/Norm(ini.Pxyz))*Qvec/Norm(Qvec);return 180.0 / PI * atan(sqrt(1 - cosw * cosw)/cosw); break; // "angl(Q.P)[�]"
case 13: return (crossx*ini.Pxyz(1)+crossy*ini.Pxyz(2)+crossz*ini.Pxyz(3));break; //"i(MSFxMSF*).P",
case 14: return Ip* lorentzf * scale * ovallt;
                     //              "I+          ",
case 15: return Im* lorentzf * scale * ovallt;
                     //              "I-          ",
case 16: return Ip/Im;     //              "I+/I-       "
case 17: return (crossx*ini.Pxyz(1)+crossy*ini.Pxyz(2)+crossz*ini.Pxyz(3));break; //i(MSFdip x MSFdip*).P
case 18: return Ip* lorentzf * scale * ovallt;
          //Idip+
case 19: return Im* lorentzf * scale * ovallt;
         // Idip-
case 20: return Ip/Im;         // Idip+/Idip-
case 21: cosw=(ini.Pxyz/Norm(ini.Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2.0*abs(msfx*ini.Pxyz(1)+msfy*ini.Pxyz(2)+msfz*ini.Pxyz(3))/sinw2;break;//|MSF.P|/sin^2(angl(Q,P)
case 22: cosw=(ini.Pxyz/Norm(ini.Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2.0*abs(msfdipx*ini.Pxyz(1)+msfdipy*ini.Pxyz(2)+msfdipz*ini.Pxyz(3))/sinw2;break;//|MSFdip.P|/sin^2(angl(Q,P)
case 23: R= Ip/Im;     //   I+/I-
         cosw=(ini.Pxyz/Norm(ini.Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))-sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
case 24: R= Ip/Im;     //   I+/I-
         cosw=(ini.Pxyz/Norm(ini.Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))+sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
case 25:  R= Ip/Im;         // Idip+/Idip-
          cosw=(ini.Pxyz/Norm(ini.Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
          return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))-sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
case 26:  R= Ip/Im;         // Idip+/Idip-
          cosw=(ini.Pxyz/Norm(ini.Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
          return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))+sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
case 27: return inuc;  //"Inuc(2t) ",    //  27
case 28: return Imag;  //"Imag(2t) ",    //  28
case 29: return inuc+Imag;  //"Itot(2t) ",    //  29
case 30: return Imagdip;  //"Imag_dip(2t) ",    //  30
case 31: return hkl(1);  //"h ",    //  31
case 32: return hkl(2);  //"k ",    //  32
case 33: return hkl(3);  //"l ",    //  33
case 34: return d;  //"d[A]  ",    //  34             
case 35: return Q;  //"|Q|[1/A] ",    //  35
case 36: return 2*theta;  //"2theta ",    //  36
case 37: return Qvec(1);//"Qi[1/A] ",    //  37    euclidean components of scattering vector 
case 38: return Qvec(2);//"Qj[1/A] ",    //  38    with j||b, k||(a x b) and i normal to k and j
case 39: return Qvec(3); //"Qk[1/A] ",    //  39
case 40: return ini.T;  //"T[K] ",    //  40
case 41: return ini.H(1);  //"Ha[T] ",    //  41
case 42: return ini.H(2);  //"Hb[T] ",    //  42
case 43: return ini.H(3);  //"Hc[T] ",    //  43
case 44:  // transform hkl to primitive lattice
          //q=ini.r1sr2sr3s.Transpose()*hkl1;
          hklprim=0;for(int j=1;j<=3;++j){hklprim+=ini.r1s(j)*hkl(j);}
          return hklprim; //"hprim ",    //  44
case 45:  // transform hkl to primitive lattice
          //q=ini.r1sr2sr3s.Transpose()*hkl1;
          hklprim=0;for(int j=1;j<=3;++j){hklprim+=ini.r2s(j)*hkl(j);}
          return hklprim; //"kprim ",    //  45
case 46:  // transform hkl to primitive lattice
          //q=ini.r1sr2sr3s.Transpose()*hkl1;
          hklprim=0;for(int j=1;j<=3;++j){hklprim+=ini.r3s(j)*hkl(j);}
          return hklprim;//"lprim ",    //  46
case 47: return inuc+Imagdip; //"Itotdip ",    //  47
   default: fprintf(stderr,"Error mcdiff: unknown column code\n");exit(EXIT_FAILURE);   
         }
                        
return 0;
}

int getint(inimcdiff & ini,int hi,int ki,int li,Vector & hkl,
 float scale,float & d,
 float & Imag,float & Imagdip,float & inuc,float * outn,complex <double> & mqx,
 complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,
 complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2)
{
// this routine calculates the intensity of elastic neutrons for a reflection (hi ki li)
//
// input:
// (*jjjpar[1...n]).xyz(1..3)         atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// (*ini.jjjpars[1...n]).DWF              debye waller factors [A^2]
// (*ini.jjjpars[1...n]).SLR,SLI          nuclear scattering length[10^-12cm]
// (*ini.jjjpars[1...n]).mom(1..3)(45)(67)(89)        atomic magnetic moment Ma Mb Mc [mb] and (if input) Sa La Sb Lb Sc Lc
//                                       ' (with respect to coordinates 1,2,3=yzx)
// (*ini.jjjpars[1...n]).gj		      Lande factor
//(*ini.jjjpars[1...n]).FF_type      1  -2    +2   +3   -3
            // (*ini.jjjpars[1...n]).FF_type // code for indicating if ion is FF_type=1: nonmagnetic,
            //=-2 GO BEYOND: using MQ function, for dipole intensities use mom(1-3)
            //                   rare earth expression (if gJ>0), spin formfactor only (if gJ=0)
            //=+2 DIPOLE ONLY: rare earth (if gJ>0), spin formfactor only (if gJ=0)
            //=+3 DIPOLE ONLY: gJ=0,general L and S moments given, use dipole approximation 
            //         and separate formfactor for spin and orbital moment
            //=-3 GO BEYOND: using MQ function, go beyond dipole approximation. 
            //         for dipole use L and S stored in mom(4-9)
// (*ini.jjjpars[1...n]).magFFj0(1..7)         formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// (*ini.jjjpars[1...n]).magFFj2(1..7)         formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// (*ini.jjjpars[1...n]).magFFj4(1..7)         formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// (*ini.jjjpars[1...n]).magFFj6(1..7)         formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// (*ini.jjjpars[1...n]).Zc		         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// (*ini.jjjpars[1...n]).eigenstates(1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)
// thetamax  			      maximum theta value, if theta larger, routine returns false
// ini.rez1,ini.rez2,ini.rez3         vectors of reciprocal lattice
// scale                              scaling factor for intensity
// T				         temperature [K] (needed to go beyond dipole approx)
// lambda                             wavelength[A]
// ovalltemp                          overall temperature factor [A^2]
// ini.lorenz                             code for lorentzfactor to be used
// n                                  number of atoms per unit cell
// Pxyz                               Projection Vector

// output:
// d                                  d spacing in A
// Imag, Imagdip, inuc                scattering intensities
// outn[4,5,6,10,11]                  output column 4,5,6,10,11
// mx,my,mz,mxmy,mxmz,mymz,mx2my2mz2[].. fouriertransform of momentunitvectors (for mag xray scattering)


            double s,Q,FQ,FQL,sintheta,qr,sin2theta,ovallt,mux,muy,muz;
            float Theta;
            int i;
            Vector Qvec(1,3);
            //calculate d spacing and intensity for h,k,l triple (d,intmag,ikern)************
            Qvec=ini.rez1*(double)(hi) + ini.rez2*(double)(ki)  + ini.rez3*(double)(li) ;
            Q = Norm(Qvec); //dspacing
            d = 2.0 * PI / Q;
            s=0.5 / d;
	    sintheta = ini.lambda * s;
            if (sintheta >= sin(ini.thetamax / 180 * PI)) return false;
               Theta = 180 / PI * atan(sintheta / sqrt(1 - sintheta * sintheta));
               //nuclear(|nsfr+i nsfc|^2) and magnetic structure factor(msf) calculation
               complex <double> nsf=0;
               complex <double> msfx=0,msfdipx=0;
               complex <double> msfy=0,msfdipy=0;
               complex <double> msfz=0,msfdipz=0;
               complex <double> im(0,1);
               for(i=1;i<=ini.n;++i){
                                 complex <double> scl((*ini.jjjpars[i]).SLR,(*ini.jjjpars[i]).SLI);
                                 qr=hi*(*ini.jjjpars[i]).xyz(1)+ki*(*ini.jjjpars[i]).xyz(2)+li*(*ini.jjjpars[i]).xyz(3);

                                 //nuclear structure factor nsfr,nsfc
                                 nsf+=scl*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);

                                //magnetic structure factors
                             //J[i]         1   0    -1   -2   -3
                             //FF_type      1  -2    +2   +3   -3
                               // if(J[i]<=0)   // i.e. atom is magnetic
                             if((*ini.jjjpars[i]).FF_type!=1){

                                             // formfactor F(Q)
                                             //if(J[i]==-1)
                                               if((*ini.jjjpars[i]).FF_type==+2)
                                               {if((*ini.jjjpars[i]).gJ==0)(*ini.jjjpars[i]).gJ=2.0;} // set gJ to 2 in case it is zero (non rare earth)
                                                                                                        // so that we get spin only formfactor
                                             FQ = (*ini.jjjpars[i]).F(Q); //rare earth

                                             //if(J[i]==0) // go beyond dipole approximation for rare earth
                                               if((*ini.jjjpars[i]).FF_type==-2){
                                                         ComplexVector MQ(1,3);(*ini.jjjpars[i]).MQ(MQ,Qvec);
					               msfx+=0.5*MQ(1)*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);//MQ(123)=MQ(xyz)
					               msfy+=0.5*MQ(2)*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					               msfz+=0.5*MQ(3)*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
                                                       msfdipx+=(*ini.jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);// mom(123)=mom(abc)=mom(yzx)
					               msfdipy+=(*ini.jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					               msfdipz+=(*ini.jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                }
 					      //if(J[i]==-1)// dipole approximation - use magnetic moments and rare earth formfactor
                                               if((*ini.jjjpars[i]).FF_type==+2){
                                                           //                        for transition metals always set gJ=2 (spin only moment)
                                                        msfx+=(*ini.jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfy+=(*ini.jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfz+=(*ini.jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
                                                        msfdipx+=(*ini.jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipy+=(*ini.jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*ini.jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
//printf (" hi=%i ki=%i li=%i Mx=%g My=%g Mz=%g qr=%g msfx=%g msfy=%g msfz=%g\n",hi,ki,li,(*ini.jjjpars[i]).mom(1),(*ini.jjjpars[i]).mom(2),(*ini.jjjpars[i]).mom(3),qr,msfx,msfy,msfz);
					               }
					      //if(J[i]==-2)// dipole approximation - use S and L moments (only if gJ=0)
                                               if((*ini.jjjpars[i]).FF_type==+3){
                                                        FQL = (*ini.jjjpars[i]).F(-Q); // orbital formfactor
                                                        msfx+=(*ini.jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // spin FF
					                msfy+=(*ini.jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfz+=(*ini.jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfx+=(*ini.jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                msfy+=(*ini.jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfz+=(*ini.jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
                                                        msfdipx+=(*ini.jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // spin FF
					                msfdipy+=(*ini.jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*ini.jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipx+=(*ini.jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                msfdipy+=(*ini.jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*ini.jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					               }
                                     //if(J[i]==-3) // go beyond dipole approximation for gJ=0 (intermediate coupling)
                                     if((*ini.jjjpars[i]).FF_type==-3){
                                                       ComplexVector MQ(1,3);(*ini.jjjpars[i]).MQ(MQ,Qvec);
                                             FQL = (*ini.jjjpars[i]).F(-Q); // orbital formfactor
                                                        msfdipx+=(*ini.jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // spin FF
					                msfdipy+=(*ini.jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*ini.jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipx+=(*ini.jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                msfdipy+=(*ini.jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*ini.jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
                                                       if (Q<SMALLPOSITIONDEVIATION){// for Q=0 put dipole results, because M(Q) givs NaN
                                                                 msfx+=(*ini.jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // spin FF
					                         msfy+=(*ini.jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                         msfz+=(*ini.jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                         msfx+=(*ini.jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                         msfy+=(*ini.jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                         msfz+=(*ini.jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
                                                           }
					               else{
                                                            msfx+=0.5*MQ(1)*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);//MQ(123)=MQ(xyz)
					                    msfy+=0.5*MQ(2)*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
					                    msfz+=0.5*MQ(3)*exp(-2*PI*qr*im)*(*ini.jjjpars[i]).debyewallerfactor(Q);
                                                           }
					               }
// myPrintVector(stdout,(*ini.jjjpars[i]).mom);//equivalent to moment ...

                                                         mux=(*ini.jjjpars[i]).mom(1); // this is still here because correlation functions are calculated
                                                         muy=(*ini.jjjpars[i]).mom(2); // only for orhtogonal lattices (see printeln sub) and so we take
                                                         muz=(*ini.jjjpars[i]).mom(3); // the convention of the mcdiff program (a||x,b||y,c||z)
                                                         mqx+=mux*exp(-2*PI*qr*im); // moreover think of using jjjpar RXMS function to calculate moments ???
                                                         mqy+=muy*exp(-2*PI*qr*im);
                                                         mqz+=muz*exp(-2*PI*qr*im);
                                                         mqx2+=mux*mux*exp(-2*PI*qr*im);
                                                         mqy2+=muy*muy*exp(-2*PI*qr*im);
                                                         mqz2+=muz*muz*exp(-2*PI*qr*im);
                                                         mqxy+=mux*muy*exp(-2*PI*qr*im);
                                                         mqxz+=mux*muz*exp(-2*PI*qr*im);
                                                         mqyz+=muy*muz*exp(-2*PI*qr*im);
                             }
                                }

             //magnetic structure factors + polarisation factor===>msf
             double msf2,msf2dip;
             msf2 = norm(msfx)+norm(msfy)+norm(msfz);
             msf2dip = norm(msfdipx)+norm(msfdipy)+norm(msfdipz);

            msf2 -=  2 * Qvec(1) * Qvec(2) / Q / Q * (real(msfx) * real(msfy) + imag(msfx) * imag(msfy));
            msf2 -=  2 * Qvec(1) * Qvec(3) / Q / Q * (real(msfx) * real(msfz) + imag(msfx) * imag(msfz));
            msf2 -=  2 * Qvec(2) * Qvec(3) / Q / Q * (real(msfy) * real(msfz) + imag(msfy) * imag(msfz));

            msf2 -=  Qvec(1) * Qvec(1) / Q / Q * norm(msfx);
            msf2 -=  Qvec(2) * Qvec(2) / Q / Q * norm(msfy);
            msf2 -=  Qvec(3) * Qvec(3) / Q / Q * norm(msfz);

            msf2dip -=  2 * Qvec(1) * Qvec(2) / Q / Q * (real(msfdipx) * real(msfdipy) + imag(msfdipx) * imag(msfdipy));
            msf2dip -=  2 * Qvec(1) * Qvec(3) / Q / Q * (real(msfdipx) * real(msfdipz) + imag(msfdipx) * imag(msfdipz));
            msf2dip -=  2 * Qvec(2) * Qvec(3) / Q / Q * (real(msfdipy) * real(msfdipz) + imag(msfdipy) * imag(msfdipz));

            msf2dip -=  Qvec(1) * Qvec(1) / Q / Q * norm(msfdipx);
            msf2dip -=  Qvec(2) * Qvec(2) / Q / Q * norm(msfdipy);
            msf2dip -=  Qvec(3) * Qvec(3) / Q / Q * norm(msfdipz);
             // alternative procedure: project msf normal to Q and then take norm:
             // msfperp= msf - Q (Q.msf)/Q^2
            complex <double> Qmsf;
            complex <double> Qmsfdip;
            Qmsf=Qvec(1)*msfx+Qvec(2)*msfy+Qvec(3)*msfz; Qmsf/=Q;
            msfx=msfx-Qvec(1)*Qmsf/Q;
            msfy=msfy-Qvec(2)*Qmsf/Q;
            msfz=msfz-Qvec(3)*Qmsf/Q;

            Qmsfdip=Qvec(1)*msfdipx+Qvec(2)*msfdipy+Qvec(3)*msfdipz; Qmsfdip/=Q;
            msfdipx=msfdipx-Qvec(1)*Qmsfdip/Q;
            msfdipy=msfdipy-Qvec(2)*Qmsfdip/Q;
            msfdipz=msfdipz-Qvec(3)*Qmsfdip/Q;

            if (fabs((norm(msfx)+norm(msfy)+norm(msfz)-fabs(msf2))/(fabs(msf2)+0.01))>0.01){fprintf(stderr,"Q=(%g %g %g) msf^2=%g |msfperp|^2=%g\n",Qvec(1),Qvec(2),Qvec(3),msf2,norm(msfx)+norm(msfy)+norm(msfz));
                                                                   fprintf(stderr,"ERROR mcdiff 1(%i %i %i): internal calculation of MSF wrong, contact Martin Rotter\n",hi,ki,li);exit(EXIT_FAILURE);}
            msf2=fabs(norm(msfx)+norm(msfy)+norm(msfz));
            if (fabs((norm(msfdipx)+norm(msfdipy)+norm(msfdipz)-fabs(msf2dip))/(fabs(msf2dip)+0.01))>0.01){fprintf(stderr,"Q=(%g %g %g) msfdip^2=%g |msfdipperp|^2=%g\n",Qvec(1),Qvec(2),Qvec(3),msf2dip,norm(msfdipx)+norm(msfdipy)+norm(msfdipz));
                                                                   fprintf(stderr,"ERROR mcdiff (%i %i %i)dipint: internal calculation of MSF wrong, contact Martin Rotter\n",hi,ki,li);exit(EXIT_FAILURE);}
            msf2dip=fabs(norm(msfdipx)+norm(msfdipy)+norm(msfdipz));
            
            //lorentzfactor*************************************************************
            float lorentzf=1;
            sin2theta = 2.0 * sintheta * sqrt(1.0 - sintheta * sintheta);
            if(ini.lorenz == 0){lorentzf = 1.0;} // no lorentzfactor
            if(ini.lorenz == 1){lorentzf = 1.0 / sin2theta / sin2theta;} // powder flat sample
            if(ini.lorenz == 2){lorentzf = 1.0 / sin2theta / sintheta;}  // powder cyl. sample
            if(ini.lorenz == 3){lorentzf = 1.0 / sin2theta;}             //single crystal
            if(ini.lorenz == 4){lorentzf = d * d * d;}      //TOF powder cyl sample... log scaled d-pattern
            if(ini.lorenz == 5){lorentzf = d * d * d * d;}  //TOF powder cyl sample... d-pattern

             //overall temperature factor*************************************************
             ovallt = exp(-2 * ini.ovalltemp * (sintheta * sintheta / ini.lambda / ini.lambda));
             //***************************************************************************

             //A)nuclear intensity
            inuc = abs(nsf) * abs(nsf) * lorentzf * scale * ovallt;

             //B)magnetic intensity
            Imag = msf2 * 3.65 / 4 / PI * lorentzf * scale * ovallt;
            Imagdip = msf2dip * 3.65 / 4 / PI * lorentzf * scale * ovallt;

             // output user defined columns 
            float msf2fl=msf2,msf2dipfl=msf2dip;
            
                              //  millerindizes refering to kristallographische einheitszelle
                              hkl(1)=Qvec*ini.rtoijk.Column(1);
                              hkl(2)=Qvec*ini.rtoijk.Column(2);
                              hkl(3)=Qvec*ini.rtoijk.Column(3); 
                              hkl/=2.0*PI;
            for(int i=1;i<= ini.nofoutputcolumns;++i)outn[i] =setcoloutput(ini.colcod[i],ini,scale,ovallt,lorentzf,nsf,msf2fl,msf2dipfl,msfx,msfy,msfz,msfdipx,msfdipy,msfdipz,Qvec,d,Theta,Q,inuc,Imag,Imagdip,hkl);

return true;
}

void neutint(inimcdiff & ini,int code,int & m,Vector *  hkl,float * D,
             float * totint,float ** out,complex <double>*mx,
             complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,
             complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2)
{//****************************************************************************
// this routine calculates the intensity of elastic neutrons
// for a given magnetic unit cell (crystal axis orthogonal)
// the magnetic scattering is treated in the dipole approximation
// input:
// code                               governs if (>0)a list of hkl given in hkl[] or (0)all hkls should be generated
// ini.lambda                             wavelength[A]
// ini.ovalltemp                          overall temperature factor [A^2]
// ini.r1(1..3),r2(),r3()                 vectors of primitive (magnetic) unit cell[A]
// ini.n                                  number of atoms per unit cell
// (*ini.jjjpar[1...n]).xyz(1..3)         atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// (*ini.jjjpars[1...n]).DWF              debye waller factors [A^2]
// (*ini.jjjpars[1...n]).SLR,SLI          nuclear scattering length[10^-12cm]
// (*ini.jjjpars[1...n]).mom(1..3)(45)(67)(89)        atomic magnetic moment Ma Mb Mc [mb] and (if input) Sa La Sb Lb Sc Lc
//                                       ' (with respect to coordinates 1,2,3=yzx)
// (*ini.jjjpars[1...n]).gj		      Lande factor
//(*ini.jjjpars[1...n]).FF_type      1  -2    +2   +3   -3
            // (*ini.jjjpars[1...n]).FF_type // code for indicating if ion is FF_type=1: nonmagnetic,
            //=-2 GO BEYOND: using MQ function, for dipole intensities use mom(1-3)
            //                   rare earth expression (if gJ>0), spin formfactor only (if gJ=0)
            //=+2 DIPOLE ONLY: rare earth (if gJ>0), spin formfactor only (if gJ=0)
            //=+3 DIPOLE ONLY: gJ=0,general L and S moments given, use dipole approximation 
            //         and separate formfactor for spin and orbital moment
            //=-3 GO BEYOND: using MQ function, go beyond dipole approximation. 
            //         for dipole use L and S stored in mom(4-9)
// (*ini.jjjpars[1...n]).magFFj0(1..7)         formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// (*ini.jjjpars[1...n]).magFFj2(1..7)         formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// (*ini.jjjpars[1...n]).magFFj4(1..7)         formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// (*ini.jjjpars[1...n]).magFFj6(1..7)         formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// T				         temperature [K] (needed to go beyond dipole approx)
// (*ini.jjjpars[1...n]).Zc		         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// (*ini.jjjpars[1...n]).eigenstates(1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)
// Pxyz                                  Projection Vector

// output
// m                                  number of calculated reflections
// hkl[1...m](1..3)                   hkl values with respect to crystallographic unit cell
// D[1...m]                           d spacing
// theta[1...m]                       scattering angle theta
// totint[1...m]                      total intenstiy (nuclear+  magnetic )
// out[1,..,12][1...m]            output column 1-12
// mx,my,mz,mxmy,mxmz,mymz,mx2my2mz2[].. fouriertransform of momentunitvectors (for mag xray scattering)

//****experimental parameters*************************************************
float scale,inuc;
float Imag,Imagdip,d;
float outn[ ini.nofoutputcolumns+1];
scale = 1 /(double)(ini.n) /(double)(ini.n); // scalingfactor of intensities
//***************************************************************************


D[0]=100000;
int i;
//calculate reciprocal lattice vectors from r1,r2,r3
  Vector nmin(1,3),nmax(1,3);
  
if(code==0){ m = 0;// reset m
 double qmax;//,rr;
 int msort,hi,ki,li;//htrue,ktrue,ltrue,
 qmax = 4.0 * PI * sin(ini.thetamax / 180 * PI) / ini.lambda;
// rr=r1*r1; hmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
// rr=r2*r2; kmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
// rr=r3*r3; lmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);

     Matrix pstar(1,3,1,3);  // inserted 10.5.10 MR to make sure all hkls are probed
     //pstar=(rez1,rez2,rez3);
     for(i=1;i<=3;++i){pstar(i,1)=ini.rez1(i);pstar(i,2)=ini.rez2(i);pstar(i,3)=ini.rez3(i);}
     nlimits_calc(nmin, nmax, qmax, pstar);
     // problem: we want to find all lattice vectors Rn=ni*ai which are within a
     // sphere of radius r from the origin (ai = column vectors of matrix a)
     // this routine returns the maximum and minimum values of ni i=1,2,3
     // by probing the corners of a cube
      for (hi=(int)nmin(1);hi<=nmax(1);++hi){
       for (ki=(int)nmin(2);ki<=nmax(2);++ki){
        for (li=(int)nmin(3);li<=nmax(3);++li){


        if(hi!=0||li!=0||ki!=0)//{htrue=1;ktrue=1;ltrue=1;} //goto 30
         {  complex <double> mqx=0,mqx2=0,mqxy=0;
            complex <double> mqy=0,mqy2=0,mqxz=0;
            complex <double> mqz=0,mqz2=0,mqyz=0;
            Vector hkl1(1,3);
          if(getint(ini,hi,ki,li,hkl1,scale,d,Imag,Imagdip,inuc,outn,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2))
          {// reflection was found below thetamax....
           //printf("%g %g %g inuc=%g\n",hkl1(1),hkl1(2),hkl1(3),inuc);
            //sort according to descending d spacing
             if((Imag + inuc) > SMALLINTENSITY||Imagdip > SMALLINTENSITY||abs(mqx)*sqrt(scale)>SMALLINTENSITY
              ||abs(mqy)*sqrt(scale)>SMALLINTENSITY||abs(mqz)*sqrt(scale)>SMALLINTENSITY
              ||abs(mqx2)*sqrt(scale)>SMALLINTENSITY||abs(mqy2)*sqrt(scale)>SMALLINTENSITY||abs(mqz2)*sqrt(scale)>SMALLINTENSITY
              ||abs(mqxy)*sqrt(scale)>SMALLINTENSITY||abs(mqxz)*sqrt(scale)>SMALLINTENSITY||abs(mqyz)*sqrt(scale)>SMALLINTENSITY){

               ++m; if(m > MAXNOFREFLECTIONS){fprintf(stderr,"ERROR mcdiff: out of memory - too many reflections - chose smaller thetamax or recompile program with larger MAXNOFREFLECTIONS\n");exit(EXIT_FAILURE);}
               msort = m;
               while(D[msort-1]<=d){
                totint[msort] = totint[msort - 1];
                for(int i=1;i<= ini.nofoutputcolumns;++i)out[i][msort] = out[i][msort - 1];
                D[msort] = D[msort - 1];
                 
                if(ini.colcod[0]<0){
                 hkl[msort] = hkl[msort - 1];
                 mx[msort] = mx[msort - 1];
                 my[msort] = my[msort - 1];
                 mz[msort] = mz[msort - 1];
                 mxmy[msort] = mxmy[msort - 1];
                 mxmz[msort] = mxmz[msort - 1];
                 mymz[msort] = mymz[msort - 1];
                 mx2[msort] = mx2[msort - 1];
                 my2[msort] = my2[msort - 1];
                 mz2[msort] = mz2[msort - 1];
                                    }  
                  --msort;
               }
               totint[msort] = Imag +inuc;
               for(int i=1;i<= ini.nofoutputcolumns;++i)out[i][msort] = outn[i];                

                D[msort]=d;
               hkl[msort]=hkl1; 
               if(ini.colcod[0]<0){
               mx[msort]=(double)sqrt(scale)*mqx;
               my[msort]=(double)sqrt(scale)*mqy;
               mz[msort]=(double)sqrt(scale)*mqz;
               mxmy[msort]=(double)sqrt(scale)*mqxy;
               mxmz[msort]=(double)sqrt(scale)*mqxz;
               mymz[msort]=(double)sqrt(scale)*mqyz;
               mx2[msort]=(double)sqrt(scale)*mqx2;
               my2[msort]=(double)sqrt(scale)*mqy2;
               mz2[msort]=(double)sqrt(scale)*mqz2;
                                  }  
              
              }
            }
          }
   }}// NEXT li NEXT ki
   printf("%i %s",100* (hi-(int)nmin(1))/((int)nmax(1)-(int)nmin(1)),"%");
   print_time_estimate_until_end((nmax(1)-hi)/(hi-(int)nmin(1)+1));
   fflush(stdout);
  }
 printf("\n");
 }
else
 {for(i=1;i<=m;++i){
                complex <double> mqx=0,mqx2=0,mqxy=0;
                complex <double> mqy=0,mqy2=0,mqxz=0;
                complex <double> mqz=0,mqz2=0,mqyz=0;
          if(fabs(rint(hkl[i](1))-hkl[i](1))>SMALLPOSITIONDEVIATION||fabs(rint(hkl[i](2))-hkl[i](2))>SMALLPOSITIONDEVIATION||fabs(rint(hkl[i](3))-hkl[i](3))>SMALLPOSITIONDEVIATION)
          {Imag=0;inuc=0;Imagdip=0; // i.e. here we treat the case where input hkl does not correspond to a reciprocal lattice point
          Vector Qvec(1,3);
            //calculate d spacing  ************
            Qvec=ini.rez1*hkl[i](1) + ini.rez2*hkl[i](2)  + ini.rez3*hkl[i](3) ;
             //   printf("%g %g %g d=%g\n",hkl[i](1),hkl[i](2),hkl[i](3),d); 
             double Q,s,sintheta,sin2theta;
             float msf2=0,msf2dip=0;
             complex<double> nsf=0,msfx=0,msfy=0,msfz=0,msfdipx=0,msfdipy=0,msfdipz=0;
            Q= Norm(Qvec);d = 2.0 * PI / Q; //dspacing
            s=0.5 / d;
	    sintheta = ini.lambda * s;
            if (sintheta >= sin(ini.thetamax / 180 * PI)) {fprintf(stderr,"ERROR mcdiff: theta for reflection number %i above thetamax=%g\n",i,ini.thetamax);exit(1);}
             float  Theta ; Theta = 180 / PI * atan(sintheta / sqrt(1 - sintheta * sintheta));
                          double ovallt;ovallt = exp(-2 * ini.ovalltemp * (sintheta * sintheta / ini.lambda / ini.lambda));
             float lorentzf=1;
            sin2theta = 2.0 * sintheta * sqrt(1.0 - sintheta * sintheta);
            if(ini.lorenz == 0){lorentzf = 1.0;} // no lorentzfactor
            if(ini.lorenz == 1){lorentzf = 1.0 / sin2theta / sin2theta;} // powder flat sample
            if(ini.lorenz == 2){lorentzf = 1.0 / sin2theta / sintheta;}  // powder cyl. sample
            if(ini.lorenz == 3){lorentzf = 1.0 / sin2theta;}             //single crystal
            if(ini.lorenz == 4){lorentzf = d * d * d;}      //TOF powder cyl sample... log scaled d-pattern
            if(ini.lorenz == 5){lorentzf = d * d * d * d;}  //TOF powder cyl sample... d-pattern
                              //  millerindizes refering to kristallographische einheitszelle
                              hkl[i](1)=Qvec*ini.rtoijk.Column(1);
                              hkl[i](2)=Qvec*ini.rtoijk.Column(2);
                              hkl[i](3)=Qvec*ini.rtoijk.Column(3); 
                              hkl[i]/=2.0*PI;
          for(int j=1;j<= ini.nofoutputcolumns;++j)outn[j]=
           setcoloutput(ini.colcod[j],ini,scale,ovallt,lorentzf,nsf,msf2,msf2dip,msfx,msfy,
                        msfz,msfdipx,msfdipy,msfdipz,Qvec,d,Theta,Q,inuc,Imag,Imagdip,hkl[i]);
          }
          else
          { if(!getint(ini,(int)hkl[i](1),(int)hkl[i](2),(int)hkl[i](3),hkl[i],scale,d,Imag,Imagdip,inuc,outn,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2))
                {fprintf(stderr,"ERROR mcdiff: theta for reflection number %i above thetamax=%g\n",i,ini.thetamax);exit(1);}
          }    totint[i] = Imag+ inuc;
               for(int j=1;j<= ini.nofoutputcolumns;++j)out[j][i] = outn[j];
               D[i]=d;
               if(code==1&&ini.colcod[0]<0){
               mx[i]=(double)sqrt(scale)*mqx;
               my[i]=(double)sqrt(scale)*mqy;
               mz[i]=(double)sqrt(scale)*mqz;
               mxmy[i]=(double)sqrt(scale)*mqxy;
               mxmz[i]=(double)sqrt(scale)*mqxz;
               mymz[i]=(double)sqrt(scale)*mqyz;
               mx2[i]=(double)sqrt(scale)*mqx2;
               my2[i]=(double)sqrt(scale)*mqy2;
               mz2[i]=(double)sqrt(scale)*mqz2;
                          }
              printf("%i %s",(int)(100* (double)i/(double)m),"%");
                  }
 }

return;}

