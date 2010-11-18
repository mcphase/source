/* checkgcc.cpp
 *
 * Routine to check gcc version, and returns the version number.
 *
 * If given no arguments the program returns the GCC version number as a string MAJOR.MINOR.PATCHLEVEL
 *
 * Otherwise it assumes all arguments are version numbers and parses it to determine if the compiler 
 * used is newer or older than the argument. Argument must be given as a string MAJOR.MINOR.PATCHLEVEL
 * The MINOR and PATCHLEVEL part may be omitted (in which case they are assumed to be zero).
 *
 * If the major version is positive, the program checks if the current compiler is newer than the 
 * argument. If the major version is negative, the program checks if the current compiler is older
 * than the argument. If the comparisons are true, it returns nothing. If false, it returns "false"
 *
 */

#include<iostream>
#include<cstdlib>

int main(int argc, char *argv[])
{
   int retval = 0;
#ifdef __GNUC__
   if(argc>1)
   {
      int wantver,curver,imajor,iminor,ipatch,ich=0,imj=0,imn=0,ipc=0;
      char ch='a',cmajor[200],cminor[200],cpatch[200];
      bool df=false, ddf=false;
      while(ch!='\0')
      {
         ch = argv[1][ich++];
         if(ch=='.') 
         {
            if(!df) df=true; else ddf=true;
         }
         else
         {
            if(ddf)     cpatch[ipc++] = ch;
            else if(df) cminor[imn++] = ch;
            else        cmajor[imj++] = ch;
         }
      }
      cpatch[ipc]='\0'; cminor[imn]='\0'; cmajor[imj]='\0';
      imajor = atoi(cmajor); iminor = atoi(cminor); ipatch = atoi(cpatch);
      wantver = abs(imajor)*10000 + iminor*100 + ipatch;
      curver = __GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__;
      if((imajor<0 && curver>wantver) || (imajor>0 && curver<wantver))
      {
         std::cout << "false\n"; retval = 1;
      }
   }
   else
      std::cout << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
#else
   std::cout << "Gnu Compiler not used.\n";
#endif
   return retval;
}
