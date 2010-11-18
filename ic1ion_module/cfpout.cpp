/* cfpout.cpp
 *
 * Outputs the coefficient of fractional parentage calculated by ic1ion, for comparison with Cowan's value using the 
 * perl script testcfp.pl
 *
 * Please use the testcfp.pl script, rather than compile and invoke this program separately.
 *
 * (c) 2009 Duc Le - duc.le@ucl.ac.uk
 * This program is licensed under the GNU General Purpose License, version 2. Please see the COPYING file
 */

#include "ic1ion.hpp"
#include <iomanip>

double sign(double val) { return val<0?-1.:1.; }

int main(int argc, char *argv[])
{
   int n,l,i,j,k,sz,szp,sc;
   bool df=false;
   if(argc>1) { n = atoi(argv[1]); } else { n = 2; }
   if(argc>2) { l = atoi(argv[2]); } else { l = 3; }
   if(l==3) df=true;
   fconf conf(n,(orbital)l);
   fconf confp(n-1,(orbital)l);
   double cfp; int count, nlines;
   std::vector<double> cfps; std::vector<int> jj;

   sz = (int)conf.states.size(); szp = (int)confp.states.size();
   for(j=0; j<szp; j++) std::cout << confp.states[j].id << " "; std::cout << "\n";
   for(i=0; i<sz; i++)
   {
      cfps.clear(); jj.clear();
      for(j=0; j<szp; j++)
      {
         if(df)
         cfp = racah_cfp(n,conf.states[i].U,conf.states[i].v,conf.states[i].S2,conf.states[i].L,confp.states[j].U,confp.states[j].v,confp.states[j].S2,confp.states[j].L);
	 else cfp = racah_cfp(n,conf.states[i].v,conf.states[i].S2,conf.states[i].L,confp.states[j].v,confp.states[j].S2,confp.states[j].L);
         if(fabs(cfp)>DBL_EPSILON) {
            cfps.push_back(cfp); jj.push_back(j); }
      }

      nlines = ((int)cfps.size()-1)/4 + 1;
      for(j=0; j<nlines; j++)
      {
         std::cout << conf.states[i].id; if(conf.states[i].id.size()==2) std::cout << "   "; else std::cout << "  ";  
         if(j==(nlines-1)) std::cout << -(j+1) << "  "; else std::cout << " " << j+1 << "  ";
         for(k=0; k<4; k++)
         {
            if((j*4+k+1)>cfps.size()) std::cout << "0   0.0000000000  ";
            else 
            {
               std::cout.precision(10);
               if((jj[j*4+k]+1)>9)
               {
                  if(cfps[j*4+k]<0) std::cout  << jj[j*4+k]+1 << " "  << std::fixed << cfps[j*4+k] << "  ";
                  else              std::cout  << jj[j*4+k]+1 << "  " << std::fixed << cfps[j*4+k] << "  ";
               }
               else
               {
                  if(cfps[j*4+k]<0) std::cout  << jj[j*4+k]+1 << "  "  << std::fixed << cfps[j*4+k] << "  ";
                  else              std::cout  << jj[j*4+k]+1 << "   " << std::fixed << cfps[j*4+k] << "  ";
               }
            }
         }
         if(df) std::cout << "F  "; else std::cout << "D  "; std::cout << n << "  " << i+1 << "\n";
      }
   }

   return 0;
}
