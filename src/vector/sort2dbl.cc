/*-----------------------------------------------------------------------------*\
| Sort a vector rearranging one other vector correspondingly        sort2dbl.cc |
|                                                                               |
| Last change: 1994                                                             |
|                                                                               |
| MatPack Library Release 1.0                                                   |
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


void sort (int n, double *v, double *w)
//
// Sorts a vector of doubles v[1..n] into ascending order using the heap sort
// algorithm, while making the corresponding rearrangements in the array w[1..n].
// This is a n*ln(n) algorithm like quicksort, but much more
// compact. It is about half as quick as quicksort.
//
{
    int i,j,k,l;
    double rv,rw;
    
    if (n < 2) return;
    
    l = (n >> 1) + 1;
    k = n;
    for (;;) {
	if (l > 1) {
	    rv = v[--l];
	    rw = w[l];
	} else {
	    rv = v[k];
	    rw = w[k];
	    v[k] = v[1];
	    w[k] = w[1];
	    if (--k == 1) {
		v[1] = rv;
		w[1] = rw;
		return;
	    }
	}
	i = l;
	j = l << 1;
	while (j <= k) {
	    if (j < k && v[j] < v[j+1]) ++j;
	    if (rv < v[j]) {
		v[i] = v[j];
		w[i] = w[j];
		j += (i = j);
	    } else 
	        j = k + 1;
	}
	v[i] = rv;
	w[i] = rw;
    }
}

//-----------------------------------------------------------------------------//
