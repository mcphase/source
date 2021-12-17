/*-----------------------------------------------------------------------------*\
| Include file isolating compatability information for                 compat.h |
| various machines and compilers. Additional code required for                  |
| compatibility is isolated in the file compat.cc                               |
|                                                                               |
| Last change: Mar 11, 2000							|
|                                                                               |
| Matpack Library Release 1.5.0                                                 |
| Copyright (C) 1991-2000 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without pri
or consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

#ifndef _compat_h
#define _compat_h 1


//-----------------------------------------------------------------------------//
// for all machines
//-----------------------------------------------------------------------------//


//-----------------------------------------------------------------------------//
// gcc version egcs-2.91.66 19990314 (egcs-1.1.2 release)
// on Sun Ultra 1, SunOS Release 5.5.1
//-----------------------------------------------------------------------------//

#if defined( __GNUC__ ) && defined( __sun__ )
// absolute value for long
inline long abs (long x){return(x<0)?-x:x;}
#include <ctime>
#include <cfloat>
// undefine because they are already in <limits.h>
//#undef DBL_MAX
//#undef DBL_MIN
//#undef FLT_MAX
//#undef FLT_MIN
// for calloc() , malloc(), ...
#include <cstdlib>
// sometimes defined !
#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC 1e6
#endif
extern "C" { 
    char *econvert(double value,int ndigit,int*decpt,int*sign,char*buf);
    char *fconvert(double value,int ndigit,int*decpt,int*sign,char*buf);
    char *gconvert(double value,int ndigit,int trailing,char*buf); 
    int usleep (unsigned);
}
inline void delay (long t)
   { usleep( (unsigned) t ); }

//-----------------------------------------------------------------------------//
// GNU g++ on DEC 3100/5000 MIPS, Ultrix 4.2 (UNIX)
//-----------------------------------------------------------------------------//

#elif defined( __GNUC__ ) && defined( __MIPSEL__ )
// absolute value for long
inline long abs (long x){return(x<0)?-x:x;}
extern "C" {  
    char *ecvt (double,int,int*,int*);
    char *gcvt (double,int,char*);
}
inline char *econvert (double value,int ndigit,int *decpt,int *sign,char *buf)
   { return strcpy(buf,ecvt(value,ndigit,decpt,sign)); }
inline char *gconvert (double value,int ndigit,int trailing,char *buf)
   { return gcvt(value,ndigit,buf); }

//-----------------------------------------------------------------------------//
// GNU g++ 2.7.2 compiler on the DEC Alpha, OSF 3.2 (UNIX)
//-----------------------------------------------------------------------------//

#elif defined( __GNUC__ ) && defined( __alpha__ )
// absolute value for long
inline long abs (long x){return(x<0)?-x:x;}
#include <cfloat>
extern "C" {  
    char *ecvt (double,int,int*,int*);
    char *gcvt (double,int,char*);
    unsigned usleep (unsigned);
    char* mktemp (char*);
}

inline char *econvert (double value,int ndigit,int *decpt,int *sign,char *buf)
   { return strcpy(buf,ecvt(value,ndigit,decpt,sign)); }
inline char *gconvert (double value,int ndigit,int trailing,char *buf)
   { return gcvt(value,ndigit,buf); }
inline void delay (long t)
   { usleep( (unsigned) t ); }

//-----------------------------------------------------------------------------//
// gcc version egcs-2.91.66 19990314/Linux (egcs-1.1.2 release)
// on Intel Pentium PC
//-----------------------------------------------------------------------------//

#elif defined( __GNUC__ ) && defined( __linux__ ) 

#if __GNUC__ > 2
#include <cstdlib>
#include <cstring>
#else
// absolute value for long
inline long abs (long x){return(x<0)?-x:x;}
//#ifndef _SVID_SOURCE
//#define _SVID_SOURCE // terrible hack to get the common ecvt/gcvt/usleep
//#endif
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#endif
/*  no longer necessary with GNU g++ 2.7.2.3
  extern "C" {  
  extern char* ecvt (double __value, size_t __ndigit, int *__decpt,
		     int *__sign);
  extern char* fcvt (double __value, size_t __ndigit, int *__decpt,
		     int *__sign);
  extern char* gcvt (double __value, size_t __ndigit, char *__buf);
  void usleep(unsigned long usec);
  }
*/
inline char *econvert (double value,int ndigit,int *decpt,int *sign,char *buf)
   { return strcpy(buf, ecvt(value,(size_t)ndigit,decpt,sign)); }
inline char *gconvert (double value,int ndigit,int,char *buf)
   { return gcvt(value,(size_t)ndigit,buf); }
//inline void delay (long t)
//   { usleep( (unsigned long) t ); }


//-----------------------------------------------------------------------------//
// add for other unix systems (This covers the HP - Ferdinand Evers)
//-----------------------------------------------------------------------------//

// absolute value for long
#inline long abs (long x){return(x<0)?-x:x;}
#else 
#include <cfloat>
#include <climits>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
extern "C" { char* mktemp (char*); }
extern "C" { char* ecvt (double,int,int*,int*); }
extern "C" { char* gcvt (double,int,char*); }
inline char *econvert (double value,int ndigit,int *decpt,int *sign,char *buf)
   { return strcpy(buf, ecvt(value,ndigit,decpt,sign)); }
inline char *gconvert (double value,int ndigit,int,char *buf)
   { return gcvt(value,ndigit,buf); }
//inline void delay (long t) { usleep( (unsigned) t ); }
#endif

//-----------------------------------------------------------------------------//
// for all machines except ...
//-----------------------------------------------------------------------------//

#if ! defined( __GNUC__ ) && ! defined( __alpha__ )
		 
inline float	abs	(float x){return(x<0)?-x:x;}
inline double	abs	(double x){return(x<0)?-x:x;}

#endif

//-----------------------------------------------------------------------------//

#endif
