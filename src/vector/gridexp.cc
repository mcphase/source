/*-----------------------------------------------------------------------------*\
| gridding with the Matpack exression parser                         gridexp.cc |
|                                                                               |
| Last change: Nov 1, 1997							|
|                                                                               |
| Matpack Library Release 1.0                                                   |
| Copyright (C) 1991-1997 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.                                                |
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.                                                    |
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.                       |
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration |
| and installation of Matpack.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "matpack.h"

//-----------------------------------------------------------------------------//
// routines for initializing the expression parser
//-----------------------------------------------------------------------------//

void    xve_init (void);
int     xve_eval_float (long id, char* string, double* value, char* error);
int     xve_delete_variables (long id);

/*
extern "C" {  
    void xve_init(void);
    int  xve_eval_float(long,char*,double*,char*);
    int  xve_delete_variables(long id);
}
*/

const int MAX_TEXT = 256;

// The id for the name space: here it is set only once to the value of id = 0.
// That means only one name space is supported upto now.
static long id;

// a static string that holds the error messages
static char error_string[MAX_TEXT];

//-----------------------------------------------------------------------------//


void SwitchExpressionParser (int on_off)
//
// Initialize/Finish the expression parser.
// Must be called with argument true before any expression is evaluated.
// To remove all dynamic structures it should be called with  argument false
// after all evaluations.
//
{
    id = 0;
    static int init = true; 
    if (init && on_off) {             // initialize
	xve_init();
	init = false;
    } else if (!init && !on_off) {    // finish
	xve_delete_variables(id);
	init = true;
    }
}

//-----------------------------------------------------------------------------//
// 3-dimensional gridding
//-----------------------------------------------------------------------------//

char* GridExpression (MpMatrix3d<double> &M, 
		      char* expression, 
		      double xstart, double xend, 
		      double ystart, double yend,
		      double zstart, double zend,
		      ProgressMeterType user_prog_meter,
		      void *user_prog_meter_args)
//
// Evaluate an expression w = f(x,y,z,w) given in the string "expression"
// on a grid that is supplied by the 3-dimensional matrix "M". That means the
// matrix row range corresponds to the grid in x-direction, the column range
// corresponds to the y-direction and the slice range cooresponds to the
// z-direction.
// The x,y,z-values are sampled within [xstart..xend], [ystart..yend] and
// [zstart..zend] respectively. 
// The string can contain an arbitrary expression as a function of (x,y,z,w).
// For each grid point the variables "x", "y", "z" are set to their 
// value at the grid point and "w" is set to the previous value M(x,y,z) at
// the grid point. After evaluating the expression the matrix element
// is set to the resulting value.
//
// Examples for expressions:
// 	"sin(x)" 
//	"w += cos(x) * sin(y) * j0(z) * sqrt(w)"
//
// Return Values:
//	The NULL pointer is returned on success, otherwise
//      a pointer to the (static) error message is returned.
//
{
    char  tmp_string[MAX_TEXT];

    int   xlo = M.Rlo(),
	  xhi = M.Rhi(),
	  ylo = M.Clo(),
	  yhi = M.Chi(),
	  zlo = M.Slo(),
	  zhi = M.Shi(),
	  x,y,z;

    double xincr,yincr,zincr,xval,yval,zval;
    double value;

    if (xend < xstart || yend < ystart || zend < zstart) {
	strcpy(error_string,"Starting point is greater than the end point");
	return error_string;
    }
   
    if (xhi-xlo > 0)
      xincr = (xend - xstart)/((double) (xhi-xlo));
    else
      xincr = (xend - xstart);
     
    if (yhi-ylo > 0)
      yincr = (yend - ystart)/((double) (yhi-ylo));
    else
      yincr = (yend - ystart);
    
    if (zhi-zlo > 0)
      zincr = (zend - zstart)/((double) (zhi-zlo));
    else
      zincr = (zend - zstart);
    
    // do parser initialization
    SwitchExpressionParser(true);

    // evaluate on grid
    for (z = zlo; z <= zhi; z++) {
	zval = (z-zlo)*zincr + zstart;

	// call the progress meter hook
	user_prog_meter(zlo,zhi,z,user_prog_meter_args);

	for (y = ylo; y <= yhi; y++) {
	    yval = (y-ylo)*yincr + ystart;
	    for (x = xlo; x <= xhi; x++) {
		xval = (x-xlo)*xincr + xstart;
		
		// set x,y,z and current w-value, evaluate expression
		sprintf(tmp_string,"x=%.15g; y=%.15g; z=%.15g; w=%.15g; %s", 
			xval, yval, zval, M(x,y,z), expression );

		// evaluate expression
		if (xve_eval_float(id,tmp_string,&value,error_string) == false) 
		    return error_string;
		
		// set matrix element
		M(x,y,z) = value;
	    }
	}
    }

    // don't switch off parser - defines may be used in later calls
    // the user has to switch it off explicitely

    // return NULL on success
    return (char*) NULL;
}


//-----------------------------------------------------------------------------//
// 2-dimensional gridding
//-----------------------------------------------------------------------------//

char* GridExpression (Matrix &M, 
		      char* expression, 
		      double xstart, double xend, 
		      double ystart, double yend,
		      ProgressMeterType user_prog_meter,
		      void *user_prog_meter_args)
//
// Evaluate an expression z = f(x,y,z) given in the string "expression"
// on a grid that is supplied by the 2-dimensional matrix "M". That means the
// matrix row range corresponds to the grid in x-direction and the column range
// corresponds to the y-direction.
// The x,y-values are sampled within [xstart..xend] and [ystart..yend]
// respectively. 
// The string can contain an arbitrary expression as a function of (x,y,z).
// For each grid point the variables "x", "y" are set to their 
// value at the grid point and "z" is set to the previous value M(x,y) at
// the grid point. After evaluating the expression the matrix element
// is set to the resulting value.
//
// Examples for expressions:
// 	"sin(x)" 
//	"z += cos(x) * sin(y) * sqrt(z)"
//
// Return Values:
//	The NULL pointer is returned on success, otherwise
//      a pointer to the (static) error message is returned.
//
{
    char  tmp_string[MAX_TEXT];

    int   xlo = M.Rlo(),
	  xhi = M.Rhi(),
	  ylo = M.Clo(),
	  yhi = M.Chi(),
	  x,y;

    double xincr,yincr,xval,yval;
    double value;

    if (xend < xstart || yend < ystart) {
	strcpy(error_string,"Starting point is greater than the end point");
	return error_string;
    }
   
    if (xhi-xlo > 0)
      xincr = (xend - xstart)/((double) (xhi-xlo));
    else
      xincr = (xend - xstart);
     
    if (yhi-ylo > 0)
      yincr = (yend - ystart)/((double) (yhi-ylo));
    else
      yincr = (yend - ystart);
    
    // do parser initialization
    SwitchExpressionParser(true);

    // evaluate on grid
    for (y = ylo; y <= yhi; y++) {
	yval = (y-ylo)*yincr + ystart;

	// call the progress meter hook
	user_prog_meter(ylo,yhi,y,user_prog_meter_args);

	for (x = xlo; x <= xhi; x++) {
	    xval = (x-xlo)*xincr + xstart;

            // set x,y and current z-value, evaluate expression
	    sprintf(tmp_string,"x=%.15g; y=%.15g; z=%.15g; %s", 
		    xval, yval, M(x,y), expression );

	    // evaluate expression
	    if (xve_eval_float(id,tmp_string,&value,error_string) == false)
		return error_string;
	    
	    // set matrix element
	    M(x,y) = value;
	}
    }

    // don't switch off parser - defines may be used in later calls
    // the user has to switch it off explicitely

    // return NULL on success
    return (char*) NULL;
}

//-----------------------------------------------------------------------------//

char* EvaluateExpression (double& value, const char* expression)
//
// Evaluate an arbitrary expression given in the string "expression".
// The numeric result is returned in "value".
//
// Return Values:
//	The NULL pointer is returned on success, otherwise
//      a pointer to the (static) error message is returned.
//
{
    double val;

    // do parser initialization - multiple switching on doesn't matter
    SwitchExpressionParser(true);

    // evaluate expression
    if (xve_eval_float(id,(char*)expression,&val,error_string) == false)
	return error_string;

    value = val;

    // don't switch off parser - defines may be used in later calls
    // the user has to switch it off explicitely

    // return NULL on success
    return (char*) NULL;
}

//-----------------------------------------------------------------------------//
