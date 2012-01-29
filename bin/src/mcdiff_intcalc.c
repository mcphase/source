// routines for calculation of intensities for program mcdiff

// different output data for columns 10 and 11
double setcoloutput(int i,float & scale, double & ovallt,float & lorentzf,complex <double> & nsf,float & msf2,float & msf2dip, Vector & Pxyz,
                   complex <double> & msfx, complex <double> & msfy, complex <double> & msfz,
                   complex <double> & msfdipx, complex <double> &msfdipy, complex <double> &msfdipz,Vector & Qvec)
{double cosw,crossx,crossy,crossz,Ip,Im,sinw2,R;
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
                            Ip-=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                            Ip+=sqrt(3.65/4/PI)*real(nsf*conj(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3))+conj(nsf)*(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3)));
                            Im=abs(nsf) * abs(nsf)+msf2 * 3.65 / 4 / PI;
                            Im+=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                            Im-=sqrt(3.65/4/PI)*real(nsf*conj(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3))+conj(nsf)*(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3)));
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
                             Ip-=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                             Ip+=sqrt(3.65/4/PI)*real(nsf*conj(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3))+conj(nsf)*(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3)));
                             Im=abs(nsf) * abs(nsf)+msf2dip * 3.65 / 4 / PI;
                             Im+=(crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3))* 3.65 / 4 / PI;
                             Im-=sqrt(3.65/4/PI)*real(nsf*conj(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3))+conj(nsf)*(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3)));
                              break;
default: break;
                    }

         switch (i) {
case 0:  return lorentzf;break;//   "LF          ",
case 1:  return abs(nsf);break;//    "|NSF|[b]    ",
case 2:  return real(nsf);break;//    "Re(NSF)[b]  ",
case 3:  return imag(nsf);break;//    "Im(NSF)[b]  ",
case 4:  return sqrt(msf2+1e-100);break;//    "|MSF|       ",
case 5:  return abs(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3));break;//    "|MSF.P|     ",
case 6:  return real(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3));break;//    "Re(MSF.P)   ",
case 7:  return imag(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3));break;//   "Im(MSF.P)   ",
case 8:  return sqrt(msf2dip+1e-100);break;//    "|MSFdip|    ",
case 9:  return abs(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3));break;//    "|MSFdip.P|  ",
case 10: return real(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3));break;//    "Re(MSFdip.P)",
case 11: return imag(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3));break;//    "Im(MSFdip.P)"
case 12: cosw=(Pxyz/Norm(Pxyz))*Qvec/Norm(Qvec);return 180.0 / PI * atan(sqrt(1 - cosw * cosw)/cosw); break; // "angl(Q.P)[°]"
case 13: return (crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3));break; //"i(MSFxMSF*).P",
case 14: return Ip* lorentzf * scale * ovallt;
                     //              "I+          ",
case 15: return Im* lorentzf * scale * ovallt;
                     //              "I-          ",
case 16: return Ip/Im;     //              "I+/I-       "
case 17: return (crossx*Pxyz(1)+crossy*Pxyz(2)+crossz*Pxyz(3));break; //i(MSFdip x MSFdip*).P
case 18: return Ip* lorentzf * scale * ovallt;
          //Idip+
case 19: return Im* lorentzf * scale * ovallt;
         // Idip-
case 20: return Ip/Im;         // Idip+/Idip-
case 21: cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2.0*abs(msfx*Pxyz(1)+msfy*Pxyz(2)+msfz*Pxyz(3))/sinw2;break;//|MSF.P|/sin^2(angl(Q,P)
case 22: cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2.0*abs(msfdipx*Pxyz(1)+msfdipy*Pxyz(2)+msfdipz*Pxyz(3))/sinw2;break;//|MSFdip.P|/sin^2(angl(Q,P)
case 23: R= Ip/Im;     //   I+/I-
         cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))-sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
case 24: R= Ip/Im;     //   I+/I-
         cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
         return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))+sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+I+/I-)/(1-I+/I-)
case 25:  R= Ip/Im;         // Idip+/Idip-
          cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
          return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))-sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|-sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
case 26:  R= Ip/Im;         // Idip+/Idip-
          cosw=(Pxyz/Norm(Pxyz))*(Qvec/Norm(Qvec));sinw2=1.00000001-cosw*cosw;
          return 2*abs(nsf)*sqrt(4*PI/3.65)*(fabs((1+R)/(1-R))+sqrt((1+2*R+R*R)/(1-2*R+R*R)-1.0/sinw2));
         //2|NSF|sqrt(4PI/3.65)(|g|+sqrt(g^2-1/sin(angl(Q,P))))_with_g=(1+Idip+/Idip-)/(1-Idip+/Idip-)
         }

return 0;
}

int getint(jjjpar ** jjjpars,int hi,int ki,int li,float thetamax,Vector rez1,Vector rez2, Vector rez3,
 float scale,double T,float lambda,float ovalltemp,int lorenz,int & n,int * J,float & d,float & Theta,
 float & Imag,float & Imagdip,float & inuc,float & OUT10,float & OUT11,complex <double> & mqx,
 complex <double> & mqy,complex <double> & mqz,complex <double> & mqxy,complex <double> & mqxz,
 complex <double> & mqyz,complex <double> & mqx2,complex <double> & mqy2,complex <double> & mqz2,
 int * colcode, Vector & Pxyz)
{
// this routine calculates the intensity of elastic neutrons for a reflection (hi ki li)
//
// input:
// (*jjjpar[1...n]).xyz(1..3)         atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// (*jjjpars[1...n]).DWF              debye waller factors [A^2]
// (*jjjpars[1...n]).SLR,SLI          nuclear scattering length[10^-12cm]
// (*jjjpars[1...n]).mom(1..3)(45)(67)(89)        atomic magnetic moment Ma Mb Mc [mb] and (if input) Sa La Sb Lb Sc Lc
//                                       ' (with respect to coordinates 1,2,3=yzx)
// (*jjjpars[1...n]).gj		      Lande factor
// J[1..n] // code for indicating if ion is J=1: nonmagnetic,
            //J[i]= 0  rare earth beyond dipole approx, but with given nonzero gJ (stevens-balcar formalism),
            //J[i]=-1 rare earth with dipole approx (if gJ>0), spin formfactor only (if gJ=0)
            //J[i]=-2 and gJ=0,general L and S moments given, use dipole approximation and separate formfactor for spin and orbital moment
            //J[i]=-3 intermediate coupling (gJ=0), go beyond dipole approximation
// (*jjjpars[1...n]).magFFj0(1..7)         formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// (*jjjpars[1...n]).magFFj2(1..7)         formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// (*jjjpars[1...n]).magFFj4(1..7)         formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// (*jjjpars[1...n]).magFFj6(1..7)         formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// (*jjjpars[1...n]).Zc		         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// (*jjjpars[1...n]).eigenstates(1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)
// thetamax  			      maximum theta value, if theta larger, routine returns false
// rez1,rez2,rez3                     vectors of reciprocal lattice
// scale                              scaling factor for intensity
// T				         temperature [K] (needed to go beyond dipole approx)
// lambda                             wavelength[A]
// ovalltemp                          overall temperature factor [A^2]
// lorenz                             code for lorentzfactor to be used
// n                                  number of atoms per unit cell
// Pxyz                               Projection Vector
// colcode                            code for output column type of data

// output:
// d                                  d spacing in A
// theta                              scattering angle
// Imag, Imagdip, inuc                scattering intensities
// OUT10                                 output column 10, default structure factor
// OUT11                           outpu column 11, default Lorentz Factor
// mx,my,mz,mxmy,mxmz,mymz,mx2my2mz2[].. fouriertransform of momentunitvectors (for mag xray scattering)


            double s,Q,FQ,FQL,sintheta,qr,sin2theta,ovallt,mux,muy,muz;
            int i;
            Vector Qvec(1,3);
            //calculate d spacing and intensity for h,k,l triple (d,intmag,ikern)************
            Qvec=rez1*(double)(hi) + rez2*(double)(ki)  + rez3*(double)(li) ;
            Q = Norm(Qvec); //dspacing
            d = 2.0 * PI / Q;
            s=0.5 / d;
	    sintheta = lambda * s;
            if (sintheta >= sin(thetamax / 180 * PI)) return false;
               Theta = 180 / PI * atan(sintheta / sqrt(1 - sintheta * sintheta));
               //nuclear(|nsfr+i nsfc|^2) and magnetic structure factor(msf) calculation
               complex <double> nsf=0;
               complex <double> msfx=0,msfdipx=0;
               complex <double> msfy=0,msfdipy=0;
               complex <double> msfz=0,msfdipz=0;
               complex <double> im(0,1);
               for(i=1;i<=n;++i){
                                 complex <double> scl((*jjjpars[i]).SLR,(*jjjpars[i]).SLI);
                                 qr=hi*(*jjjpars[i]).xyz(1)+ki*(*jjjpars[i]).xyz(2)+li*(*jjjpars[i]).xyz(3);

                                 //nuclear structure factor nsfr,nsfc
                                 nsf+=scl*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);

                                //magnetic structure factors
                                if(J[i]<=0){   // i.e. atom is magnetic

                                             // formfactor F(Q)
                                             if(J[i]==-1){if((*jjjpars[i]).gJ==0)(*jjjpars[i]).gJ=2.0;} // set gJ to 2 in case it is zero (non rare earth)
                                                                                                        // so that we get spin only formfactor
                                             FQ = (*jjjpars[i]).F(Q); //rare earth

                                             if(J[i]==0){ // go beyond dipole approximation for rare earth
                                                         ComplexVector MQ(1,3);MQ=(*jjjpars[i]).MQ(Qvec);
					               msfx+=0.5*MQ(1)*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);//MQ(123)=MQ(xyz)
					               msfy+=0.5*MQ(2)*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					               msfz+=0.5*MQ(3)*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
                                                       msfdipx+=(*jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);// mom(123)=mom(abc)=mom(yzx)
					               msfdipy+=(*jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					               msfdipz+=(*jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                }
 					      if(J[i]==-1){// dipole approximation - use magnetic moments and rare earth formfactor
                                                           //                        for transition metals always set gJ=2 (spin only moment)
                                                        msfx+=(*jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfy+=(*jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfz+=(*jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
                                                        msfdipx+=(*jjjpars[i]).mom(1)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipy+=(*jjjpars[i]).mom(2)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*jjjpars[i]).mom(3)*FQ/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					               }
					      if(J[i]==-2){// dipole approximation - use S and L moments (only if gJ=0)
                                                        FQL = (*jjjpars[i]).F(-Q); // orbital formfactor
                                                        msfx+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // spin FF
					                msfy+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfz+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfx+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                msfy+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfz+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
                                                        msfdipx+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // spin FF
					                msfdipy+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipx+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                msfdipy+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					               }
                                     if(J[i]==-3){ // go beyond dipole approximation for gJ=0 (intermediate coupling)
                                                       ComplexVector MQ(1,3);MQ=(*jjjpars[i]).MQ(Qvec);
                                             FQL = (*jjjpars[i]).F(-Q); // orbital formfactor
                                                        msfdipx+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // spin FF
					                msfdipy+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipx+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                msfdipy+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                msfdipz+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
                                                       if (Q<SMALL){// for Q=0 put dipole results, because M(Q) givs NaN
                                                                 msfx+=(*jjjpars[i]).mom(4)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // spin FF
					                         msfy+=(*jjjpars[i]).mom(6)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                         msfz+=(*jjjpars[i]).mom(8)*FQ*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                         msfx+=(*jjjpars[i]).mom(5)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q); // orbital FF
					                         msfy+=(*jjjpars[i]).mom(7)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                         msfz+=(*jjjpars[i]).mom(9)*FQL/2*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
                                                           }
					               else{
                                                            msfx+=0.5*MQ(1)*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);//MQ(123)=MQ(xyz)
					                    msfy+=0.5*MQ(2)*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
					                    msfz+=0.5*MQ(3)*exp(-2*PI*qr*im)*(*jjjpars[i]).debyewallerfactor(Q);
                                                           }
					               }
//printf("added ion %i MOUT10=(%8.6f %+8.6f i,%8.6f %+8.6f i,%8.6f %+8.6f i)\n         MSFdip=(%8.6f %+8.6f i,%8.6f %+8.6f i,%8.6f %+8.6f i)\n",i,real(msfx),imag(msfx),real(msfy),imag(msfy),real(msfz),imag(msfz),real(msfdipx),imag(msfdipx),real(msfdipy),imag(msfdipy),real(msfdipz),imag(msfdipz));
// myPrintVector(stdout,(*jjjpars[i]).mom);//equivalent to moment ...

                                                         mux=(*jjjpars[i]).mom(1); // this is still here because correlation functions are calculated
                                                         muy=(*jjjpars[i]).mom(2); // only for orhtogonal lattices (see printeln sub) and so we take
                                                         muz=(*jjjpars[i]).mom(3); // the convention of the mcdiff program (a||x,b||y,c||z)
                                                         mqx+=mux*exp(-2*PI*qr*im);
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
             float msf2,msf2dip;
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

            if (fabs((norm(msfx)+norm(msfy)+norm(msfz)-fabs(msf2))/(fabs(msf2)+0.0001))>0.01){fprintf(stderr,"Q=(%g %g %g) msf^2=%g |msfperp|^2=%g\n",Qvec(1),Qvec(2),Qvec(3),msf2,norm(msfx)+norm(msfy)+norm(msfz));
                                                                   fprintf(stderr,"ERROR mcdiff 1(%i %i %i): internal calculation of MSF wrong, contact Martin Rotter\n",hi,ki,li);exit(EXIT_FAILURE);}
            msf2=fabs(norm(msfx)+norm(msfy)+norm(msfz));
            if (fabs((norm(msfdipx)+norm(msfdipy)+norm(msfdipz)-fabs(msf2dip))/(fabs(msf2dip)+0.0001))>0.01){fprintf(stderr,"Q=(%g %g %g) msfdip^2=%g |msfdipperp|^2=%g\n",Qvec(1),Qvec(2),Qvec(3),msf2dip,norm(msfdipx)+norm(msfdipy)+norm(msfdipz));
                                                                   fprintf(stderr,"ERROR mcdiff (%i %i %i)dipint: internal calculation of MSF wrong, contact Martin Rotter\n",hi,ki,li);exit(EXIT_FAILURE);}
            msf2dip=fabs(norm(msfdipx)+norm(msfdipy)+norm(msfdipz));
            
            //lorentzfactor*************************************************************
            float lorentzf=1;
            sin2theta = 2.0 * sintheta * sqrt(1.0 - sintheta * sintheta);
            if(lorenz == 0){lorentzf = 100;} // no lorentzfactor
            if(lorenz == 1){lorentzf = 1.0 / sin2theta / sin2theta;} // powder flat sample
            if(lorenz == 2){lorentzf = 1.0 / sin2theta / sintheta;}  // powder cyl. sample
            if(lorenz == 3){lorentzf = 1.0 / sin2theta;}             //single crystal
            if(lorenz == 4){lorentzf = d * d * d;}      //TOF powder cyl sample... log scaled d-pattern
            if(lorenz == 5){lorentzf = d * d * d * d;}  //TOF powder cyl sample... d-pattern

             //overall temperature factor*************************************************
             ovallt = exp(-2 * ovalltemp * (sintheta * sintheta / lambda / lambda));
             //***************************************************************************

             //A)nuclear intensity
            inuc = abs(nsf) * abs(nsf) * lorentzf * scale * ovallt;

             //B)magnetic intensity
            Imag = msf2 * 3.65 / 4 / PI * lorentzf * scale * ovallt;
            Imagdip = msf2dip * 3.65 / 4 / PI * lorentzf * scale * ovallt;

             // output column 10
            OUT10 = setcoloutput(colcode[10],scale,ovallt,lorentzf,nsf,msf2,msf2dip,Pxyz,msfx,msfy,msfz,msfdipx,msfdipy,msfdipz,Qvec);
             // output column 11
            OUT11 = setcoloutput(colcode[11],scale,ovallt,lorentzf,nsf,msf2,msf2dip,Pxyz,msfx,msfy,msfz,msfdipx,msfdipy,msfdipz,Qvec);


return true;
}

void neutint(jjjpar ** jjjpars,int code,double T,float lambda, float thetamax, float ovalltemp,int lorenz,
             Vector r1,Vector r2,Vector r3,int & n,int * J,int & m,Vector *  hkl,float * D,float * theta,
             float * intmag,float * intmagdip,float * ikern,float * out10,float * out11,complex <double>*mx,
             complex <double>*my,complex <double>*mz,complex <double>*mxmy,complex <double>*mxmz,
             complex <double>*mymz,complex <double>*mx2,complex <double>*my2,complex <double>*mz2,
             int * colcode,Vector & Pxyz)
{//****************************************************************************
// this routine calculates the intensity of elastic neutrons
// for a given magnetic unit cell (crystal axis orthogonal)
// the magnetic scattering is treated in the dipole approximation
// input:
// code                               governs if a list of hkl given in hkl[] or all hkls should be generated
// lambda                             wavelength[A]
// ovalltemp                          overall temperature factor [A^2]
// r1(1..3),r2(),r3()                 vectors of primitive unit cell[A]
// n                                  number of atoms per unit cell
// (*jjjpar[1...n]).xyz(1..3)         atomic positional parameters dr1 dr2 dr3
//                                        '(with respect to primitive lattice)
// (*jjjpars[1...n]).DWF              debye waller factors [A^2]
// (*jjjpars[1...n]).SLR,SLI          nuclear scattering length[10^-12cm]
// (*jjjpars[1...n]).mom(1..3)(45)(67)(89)        atomic magnetic moment Ma Mb Mc [mb] and (if input) Sa La Sb Lb Sc Lc
//                                       ' (with respect to coordinates 1,2,3=yzx)
// (*jjjpars[1...n]).gj		      Lande factor
// J[1..n] // code for indicating if ion is nonmagnetic (J=1),
            //rare earth with dipole approx (J=-1),
            //rare earth beyond dipole approx, but with given nonzero gJ (stevens-balcar formalism) (J=0),
            //gJ=0,general L and S moments given, use dipole approximation and separate formfactor for spin and orbital moment (J=-2)
            //intermediate coupling (gJ=0), go beyond dipole approximation (J=-3)
// (*jjjpars[1...n]).magFFj0(1..7)         formfactor j0 for atom 1...n <j0(kr)>-terms A,a,B,b,C,c,D
// (*jjjpars[1...n]).magFFj2(1..7)         formfactor j2 for atom 1...n <j2(kr)>-terms A,a,B,b,C,c,D
//     <jl(kr)> is defined as = integral[0,inf] U^2(r) jl(kr) 4 pi r^2 dr
//     where U(r) is the Radial wave function for the unpaired electrons in the atom
// (*jjjpars[1...n]).magFFj4(1..7)         formfactor j4 for atom 1...n  (needed to go beyond dipole approx)
// (*jjjpars[1...n]).magFFj6(1..7)         formfactor j6 for atom 1...n  (needed to go beyond dipole approx)
// T				         temperature [K] (needed to go beyond dipole approx)
// (*jjjpars[1...n]).Zc		         Z-factors from Lovesey table 11.1 for Z(K) calc (needed to go beyond dipole approx)
// (*jjjpars[1...n]).eigenstates(1..2J+1,1..2J+1)   CF+MF eigenstates (needed to go beyond dipole approx)
// Pxyz                                  Projection Vector
// colcode                            code for output column type of data

// output
// m                                  number of calculated reflections
// hkl[1...m](1..3)                   hkl values
// D[1...m]                           d spacing
// theta[1...m]                       scattering angle theta
// intmag[1...m]                        magnetic intensity
// intmagdip[1...m]                     magnetic intensity in dipole approx
// ikern[1...m]                       nuclear intensity
// out10[1...m]                       output column 10, default   nuclear structurfactor |sf|
// out11[1...m]                       output column 11, default  lorentzfactor
// mx,my,mz,mxmy,mxmz,mymz,mx2my2mz2[].. fouriertransform of momentunitvectors (for mag xray scattering)

//****experimental parameters*************************************************
float scale,inuc;
float OUT10,Imag,Imagdip,Theta,OUT11,d;

scale = 1 /(double)(n) /(double)(n); // scalingfactor of intensities
//***************************************************************************


D[0] = 10000;
int i;
//calculate reciprocal lattice vectors from r1,r2,r3
  Vector rez1(1,3),rez2(1,3),rez3(1,3),nmin(1,3),nmax(1,3);
  rezcalc(r1, r2, r3, rez1, rez2, rez3);

if(code==0){ m = 0;// reset m
 double qmax;//,rr;
 int htrue,ktrue,ltrue,msort,hi,ki,li;
 qmax = 4.0 * PI * sin(thetamax / 180 * PI) / lambda;
// rr=r1*r1; hmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
// rr=r2*r2; kmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);
// rr=r3*r3; lmax =(int)( qmax / 2 / PI * sqrt(rr) + 1);

     Matrix pstar(1,3,1,3);  // inserted 10.5.10 MR to make sure all hkls are probed
     //pstar=(rez1,rez2,rez3);
     for(i=1;i<=3;++i){pstar(i,1)=rez1(i);pstar(i,2)=rez2(i);pstar(i,3)=rez3(i);}
     nlimits_calc(nmin, nmax, qmax, pstar);
     // problem: we want to find all lattice vectors Rn=ni*ai which are within a
     // sphere of radius r from the origin (ai = column vectors of matrix a)
     // this routine returns the maximum and minimum values of ni i=1,2,3
     // by probing the corners of a cube
      for (hi=(int)nmin(1);hi<=nmax(1);++hi){
       for (ki=(int)nmin(2);ki<=nmax(2);++ki){
        for (li=(int)nmin(3);li<=nmax(3);++li){


        if(hi==0&&li==0&&ki==0){htrue=1;ktrue=1;ltrue=1;} //goto 30
        else {  complex <double> mqx=0,mqx2=0,mqxy=0;
                complex <double> mqy=0,mqy2=0,mqxz=0;
                complex <double> mqz=0,mqz2=0,mqyz=0;

          if(getint(jjjpars,hi,ki,li,thetamax,rez1,rez2,rez3,scale,T,lambda,ovalltemp,lorenz,n,J,d,Theta,Imag,Imagdip,inuc,OUT10,OUT11,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2,colcode,Pxyz))
          {// reflection was found below thetamax....

            //sort according to descending d spacing
             if((Imag + inuc) > SMALLINT||Imagdip > SMALLINT||abs(mqx)*sqrt(scale)>SMALLINT
              ||abs(mqy)*sqrt(scale)>SMALLINT||abs(mqz)*sqrt(scale)>SMALLINT
              ||abs(mqx2)*sqrt(scale)>SMALLINT||abs(mqy2)*sqrt(scale)>SMALLINT||abs(mqz2)*sqrt(scale)>SMALLINT
              ||abs(mqxy)*sqrt(scale)>SMALLINT||abs(mqxz)*sqrt(scale)>SMALLINT||abs(mqyz)*sqrt(scale)>SMALLINT){

               ++m; if(m > MAXNOFREFLECTIONS){fprintf(stderr,"ERROR mcdiff: out of memory - too many reflections - chose smaller thetamax or recompile program with larger MAXNOFREFLECTIONS\n");exit(EXIT_FAILURE);}
               msort = m;
               while(D[msort-1]<=d){
                D[msort] = D[msort - 1];
                theta[msort]= theta[msort - 1];
                hkl[msort] = hkl[msort - 1];
                intmag[msort] = intmag[msort - 1];
                intmagdip[msort] = intmagdip[msort - 1];
                ikern[msort] = ikern[msort - 1];
                out10[msort] = out10[msort - 1];
                out11[msort] = out11[msort - 1];

                mx[msort] = mx[msort - 1];
                my[msort] = my[msort - 1];
                mz[msort] = mz[msort - 1];
                mxmy[msort] = mxmy[msort - 1];
                mxmz[msort] = mxmz[msort - 1];
                mymz[msort] = mymz[msort - 1];
                mx2[msort] = mx2[msort - 1];
                my2[msort] = my2[msort - 1];
                mz2[msort] = mz2[msort - 1];
                  --msort;
               }
               hkl[msort](1) = hi;hkl[msort](2) = ki; hkl[msort](3) = li;
               D[msort] = d; theta[msort] = Theta;
               intmag[msort] = Imag;intmagdip[msort] = Imagdip;  ikern[msort] = inuc;
               out10[msort] = OUT10; out11[msort] = OUT11;
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
   }}// NEXT li NEXT ki
   printf("%i %s",100* (hi-(int)nmin(1))/((int)nmax(1)-(int)nmin(1)),"%"); fflush(stdout);
  }
 printf("\n");
 }
else
 {for(i=1;i<=m;++i){
                complex <double> mqx=0,mqx2=0,mqxy=0;
                complex <double> mqy=0,mqy2=0,mqxz=0;
                complex <double> mqz=0,mqz2=0,mqyz=0;
//               printf("%g %g %g\n",hkl[i](1),hkl[i](2),hkl[i](3));
               hkl[i](1)=rint(hkl[i](1));
               hkl[i](2)=rint(hkl[i](2));
               hkl[i](3)=rint(hkl[i](3));
          if(!getint(jjjpars,(int)hkl[i](1),(int)hkl[i](2),(int)hkl[i](3),thetamax,rez1,rez2,rez3,scale,T,lambda,ovalltemp,lorenz,n,J,d,Theta,Imag,Imagdip,inuc,OUT10,OUT11,mqx,mqy,mqz,mqxy,mqxz,mqyz,mqx2,mqy2,mqz2,colcode,Pxyz))
                {fprintf(stderr,"ERROR mcdiff: theta for reflection number %i above thetamax=%g\n",i,thetamax);exit(1);}
               D[i] = d; theta[i] = Theta;
               intmag[i] = Imag;intmagdip[i] = Imagdip;  ikern[i] = inuc;
               out10[i] = OUT10; out11[i] = OUT11;
               if(code==1){
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

