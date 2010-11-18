/*-----------------------------------------------------------------------------*\
| Matpack - some utilities, error handling, ...                      matpack.cc |
|                                                                               |
| Last change: Jan 5, 1999							|
|                                                                               |
| Matpack Library Release 1.4                                                   |
| Copyright (C) 1991-1998 by Berndt M. Gammel. All rights reserved.             |
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
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "common.h"

#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <iostream> 
using namespace std;

//-----------------------------------------------------------------------------//
// define the instance of the mat class
//-----------------------------------------------------------------------------//

Mat Matpack;

//-----------------------------------------------------------------------------//
// char* Message(void) { return message; }
//
// INLINE FUNCTION
// returns the message string generated.
//
//-----------------------------------------------------------------------------//

//-----------------------------------------------------------------------------//
// Mat::Mat: CONSTRUCTOR
//-----------------------------------------------------------------------------//

Mat::Mat (void)  
{
  errnum = NoError;
  behave = AbortOnError;
  version = MP_VERSION_STRING ", " MP_COPYRIGHT_STRING;
}

//-----------------------------------------------------------------------------//
// int Mat::Behaviour (int how);
//
// INLINE FUNCTION
// Determines the behaviour of the function Error().
// For a description of possible arguments see there.
// Returns the previous behaviour setting.
//
//-----------------------------------------------------------------------------//

//-----------------------------------------------------------------------------//
// int Mat::Result (void)
//
// INLINE FUNCTION
// Returns the error status of the last Matpack function call
// and resets the error flag to NoError.
// 
//-----------------------------------------------------------------------------//

//-----------------------------------------------------------------------------//
// int Mat::Ok (void)
//
// INLINE FUNCTION
// Returns True if in the last function call there was no error
// False, otherwise. In this case you should call Matpack.Result()
// to get out what occured ! The error flag is not reset !
// 
//-----------------------------------------------------------------------------//

//-----------------------------------------------------------------------------//
// void Mat::Error (int err, char* msg, ...);
// void Mat::Error (char* msg, ...);
//
// Default error handler which is called whenever
// an error occurs in a Matpack function.
// The behaviour of this function can be controlled by
// the Matpack.Behaviour() function.
//
// The first argument 'err' sets the error flag, which
// can be tested later by Matpack.Result(). For this
// purpose it must be non-negative. Use the predefined values
// in the enum Mat::Errors for this purpose.
//
// The variable argument list defines the error message
// to be printed. It can be formated like in the printf()
// standard function. If not suppressed the error  message is sent to
// the standard file stream "cerr".
//
// If an error occurs this function does what you
// have specified with the Matpack::Behaviour(how) function.
// The following values for 'how' are possible:
//
//	Mat::AbortOnError	exit with core dump and return value 3
//	Mat::ExitOnError	exit with a return value of 1
//	Mat::SuppressWarning	dont't print warnings, but regard errors
//	Mat::DisregardError	disregard error, continue, but print warning
//	Mat::IgnoreError 	ignore error completely and don't print warning
// 
// Note: 
//	1) the last error message can always be retrieved by
//	   char* Matpack.Message(void), also if no message is printed to cerr.
//
//	2) If you specify one of the values listed above for the first argument 
// 	   'err' you will override the global behaviour for this time only.
//
//-----------------------------------------------------------------------------//

void Mat::Error (int err, const char* msg, ...)
{
    int behaviour = behave;

    // set the error number for calls of Mat::Result();
    errnum = err;

    // override the default behaviour
    if (err < 0) behaviour = err;

    // build error message
    va_list argptr;

    if (err > 0)
	sprintf(message,"Error %d: ", err);
    else
	sprintf(message,"Error: ");

    va_start(argptr,msg);
    int l = message ? strlen(message) : 0;
    vsprintf(message+l,msg,argptr);
    va_end(argptr);

    // don't bother about error any more
    if (behaviour == IgnoreError || behaviour == SuppressWarning) return;
    
    // else print error message
    cerr << message << endl;

    if (behaviour == DisregardError) return;

    // fatal error with exit
    cerr << "...exiting" << endl;

    if (behaviour == ExitOnError)
	exit(1);

    // fatal error with dump
    else if (behaviour == AbortOnError)
	abort();

    else // unknown behaviour specified
	exit(1);
}

//-----------------------------------------------------------------------------//

void Mat::Error (const char* msg, ...)
//
// Has same functionality as Error() above, but the first argument 
// defaults to 'UnspecifiedError'.
//
{
    int err = UnspecifiedError;

    int behaviour = behave;

    // set the error number for calls of Mat::Result();
    errnum = err;

    // override the default behaviour
    if (err < 0) behaviour = err;

    // build error message
    va_list argptr;

    if (err > 0)
	sprintf(message,"Error %d: ", err);
    else
	sprintf(message,"Error: ");

    va_start(argptr,msg);
    int l = message ? strlen(message) : 0;
    vsprintf(message+l,msg,argptr);
    va_end(argptr);

    // don't bother about error any more
    if (behaviour == IgnoreError || behaviour == SuppressWarning) return;
    
    // else print error message
    cerr << message << endl;

    if (behaviour == DisregardError) return;

    // fatal error with exit
    cerr << "...exiting" << endl;

    if (behaviour == ExitOnError)
	exit(1);

    // fatal error with dump
    else if (behaviour == AbortOnError)
	abort();

    else // unknown behaviour specified
	exit(1);




  // DOES NOT WORK:
  //    va_list argptr;
  //    va_start(argptr,msg);
  //    Error(UnspecifiedError,msg,argptr);
  //    va_end(argptr);

}

//-----------------------------------------------------------------------------//


void Mat::Warning (int err, const char* msg, ...)
{
    int behaviour = behave;

    // set the error number for calls of Mat::Result();
    errnum = err;

    // override the default behaviour
    if (err < 0) behaviour = err;

    // build error message
    va_list argptr;

    if (err > 0)
	sprintf(message,"Warning %d: ", err);
    else
	sprintf(message,"Warning: ");

    va_start(argptr,msg);
    int l = message ? strlen(message) : 0;
    vsprintf(message+l,msg,argptr);
    va_end(argptr);

    // don't bother about warning any more
    if (behaviour == IgnoreError || behaviour == SuppressWarning) return;
    
    // else print error message
    cerr << message << endl;
}

//-----------------------------------------------------------------------------//

void Mat::Warning (const char* msg, ...)
//
// Has same functionality as Warning() above, but the first argument 
// defaults to 'UnspecifiedError'.
//
{
    va_list argptr;
    va_start(argptr,msg);
    Warning(UnspecifiedWarning,msg,argptr);
    va_end(argptr);
}

//-----------------------------------------------------------------------------//
