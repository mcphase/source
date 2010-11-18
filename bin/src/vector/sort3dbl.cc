/*-----------------------------------------------------------------------------*\
| Sort a vector rearranging two other vectors correspondingly       sort3dbl.cc |
|                                                                               |
| Last change: Jan 9, 1997                                                      |
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


void sort (int n, double *v, double *v2, double *v3)
//
// Sorts a vector of doubles v[1..n] into ascending order using the heap sort
// algorithm, while making the corresponding rearrangements in the arrays 
// v2[1..n] and v3[1..n].
// This is a n*ln(n) algorithm like quicksort, but much more
// compact. It is about half as quick as quicksort.
//
{
    int i,j,k,l;
    double rv,rv2,rv3;
    
    if (n < 2) return;
    
    l = (n >> 1) + 1;
    k = n;
    for (;;) {
        if (l > 1) {
            rv = v[--l];
            rv2 = v2[l];
            rv3 = v3[l];
        } else {
            rv = v[k];
            rv2 = v2[k];
            rv3 = v3[k];
            v[k] = v[1];
            v2[k] = v2[1];
            v3[k] = v3[1];
            if (--k == 1) {
                v[1] = rv;
                v2[1] = rv2;
                v3[1] = rv3;
                return;
            }
        }
        i = l;
        j = l << 1;
        while (j <= k) {
            if (j < k && v[j] < v[j+1]) ++j;
            if (rv < v[j]) {
                v[i] = v[j];
                v2[i] = v2[j];
                v3[i] = v3[j];
                j += (i = j);
            } else 
                j = k + 1;
        }
        v[i] = rv;
        v2[i] = rv2;
        v3[i] = rv3;
    }
}

//-----------------------------------------------------------------------------//
