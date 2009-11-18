/*
 *	scalar.h
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

#include "functions_inline_asm.h"
#include "functions_C.h"

void setup_scalars(void);

#define printmscalar(a)		PRINT4(a)

#define valuation(x)		VAL4(x)

#define sc_add(a,b,c)		ADD4(c,b,a)

#define sc_mult(a,b,c)		MUL4(c,b,a)

static inline void sc_imult(int a, mscalar b, mscalar c)
{
	unsigned long TTTT[4];
	I_TO_4(TTTT, a);
	MUL4(c, b, TTTT);
}

#define sc_inv(a,b)		{SET4(b,a);INV4(b);}

void sc_div(mscalar a, mscalar b, mscalar c);

static inline void div_p(mscalar a, int k)
{
	int eeee = k;
	while (eeee >= 64) {
		a[0] = a[1];
		a[1] = a[2];
		a[2] = a[3];
		a[3] = 0;
		eeee = eeee - 64;
	}
	DIV4(a, eeee);
}

static inline void sc_add_replace(mscalar a, mscalar b)
{
	unsigned long TTTT[4];
	ADD4(TTTT, a, b);
	SET4(b, TTTT);
}

static inline void sc_mult_replace(mscalar a, mscalar b)
{
	unsigned long TTTT[4];
	MUL4(TTTT, a, b);
	SET4(b, TTTT);
}

static inline void sc_imult_replace(int a, mscalar b)
{
	unsigned long TTTT[4], SSSS[4];
	I_TO_4(TTTT, a);
	MUL4(SSSS, TTTT, b);
	SET4(b, SSSS);
}

#define sc_zero(a)		{a[0]=0;a[1]=0;a[2]=0;a[3]=0;}

#define sc_one(a)		{a[0]=1;a[1]=0;a[2]=0;a[3]=0;}

#define sc_copy(a,b)		SET4(b,a)

#define sc_negate(a)		NEG4(a);

#define ito_sc(a,b)		I_TO_4(b,a)

#define sc_is_zero(a)		ISZERO4(a)
