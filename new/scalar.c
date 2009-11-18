/*
 *	scalar.c
 *
 * 	Copyright 2006 Johan de Jong
 *
 *	This file is part of Frobenius
 *
 *	Frobenius is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation; either version 2 of the License, or
 *	(at your option) any later version.
 *
 *	Frobenius is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with Frobenius; if not, write to the Free Software Foundation, 
 *	Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *									*/

#include <stdio.h>
#include <stdlib.h>
#include "data.h"
#include "scalar.h"

#ifdef KIJKEN
void empty_function(void)
{
	printf("This just for breakpoints when debugging.\n");
}
#endif

/* Only called once. */
void setup_scalars(void)
{
	printf("No set up needed.\n");
}

/* Divides a by b. If b is not a unit then this assumes 	*
 * valuation(a) >= valuation(b), and the result is lifted	*
 * to an integer mod p^r.					*/
/* Does not destroy a and b.					*/
void sc_div(mscalar a, mscalar b, mscalar c)
{
	unsigned long TT[4];
	unsigned long SSS[4];
	int e;

	e = VAL4(b);
	SET4(TT, b);
	SET4(SSS, a);
	while (e >= 64) {
		TT[0] = TT[1];
		TT[1] = TT[2];
		TT[2] = TT[3];
		TT[3] = 0;
		SSS[0] = SSS[1];
		SSS[1] = SSS[2];
		SSS[2] = SSS[3];
		SSS[3] = 0;
		e = e - 64;
	}
	DIV4(TT, e);
	DIV4(SSS, e);
	INV4(TT);
	MUL4(c, SSS, TT);
}
