/*------------------------------------------------------------------------*
|  sort()                                                       sortdbl.c |
|                                                                         |
|  (C) B.M.GAMMEL, 1991                                                   |
*------------------------------------------------------------------------*/

#include "math.h"


void sort (int n, double *v)
//
// Sorts a vector of doubles v[1..n] into ascending order using the heap sort
// algorithm. This is a n*ln(n) algorithm like quicksort, but much more
// compact. It is about half as quick as quicksort.
//
{
    int i,j,k,l;
    double rv;
    
    if (n < 2) return;
    
    l = (n >> 1) + 1;
    k = n;
    for (;;) {
	if (l > 1)
	  rv = v[--l];
	else {
	    rv = v[k];
	    v[k] = v[1];
	    if (--k == 1) {
		v[1] = rv;
		return;
	    }
	}
	i = l;
	j = l << 1;
	while (j <= k) {
	    if (j < k && v[j] < v[j+1]) ++j;
	    if (rv < v[j]) {
		v[i] = v[j];
		j += (i = j);
	    } else 
	        j = k + 1;
	}
	v[i] = rv;
    }
}
