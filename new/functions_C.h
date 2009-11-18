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

#define DIV4(a,b)	div4_inline_asm_C(a,b)
#define ADD4(a,b,c)	add4_inline_asm_C(a,b,c)
#define MUL4(a,b,c)	mul4_inline_asm_C(a,b,c)
#define NEG4(a)		neg4_inline_asm_C(a)
#define VAL4(a)		val4_C(a)
#define	SET4(a,b)	{a[0]=b[0];a[1]=b[1];a[2]=b[2];a[3]=b[3];}
#define	EQUAL4(a,b)	((a[0]==b[0])&&(a[1]==b[1])&&(a[2]==b[2])&&(a[3]==b[3]))
#define ISZERO4(a)	((a[0]==0)&&(a[1]==0)&&(a[2]==0)&&(a[3]==0))
#define INV4(a)		inv4_C(a)
#define	I_TO_4(a,k)	i_to_4_C(a,k)
#define	PRINT4(a)	printf("%lu + %lu*2^64 + %lu*2^128 + %lu*2^192",\
							a[0],a[1],a[2],a[3])

static inline void i_to_4_C(unsigned long *A, int k)
{
	A[1] = 0;
	A[2] = 0;
	A[3] = 0;
	if (k >= 0) {
		A[0] = k;
	} else {
		A[0] = -k;
		NEG4(A);
	}
}

static inline int val4_C(unsigned long *A)
{
	unsigned long a, b;

        if (A[0] == 0) {
		if (A[1] == 0) {
			if (A[2] == 0) {
				if (A[3] == 0) {
					return(256);
				} else {
					a = 192;
					b = A[3];
				}
			} else {
				a = 128;
				b = A[2];
			}
		} else {
			a = 64;
			b = A[1];
		}
	} else {
		a = 0;
		b = A[0];
	}
        __asm__ ("bsfq %1,%0" : "=r" (b) : "r" (b));
        return((int) a + b);
}

static inline void inv4_C(unsigned long *A)
{
	int e;
	unsigned long B[4], C[4], D[4];

	SET4(C, A);
	A[0]--; /* No overflow as A is unit 1+2... */
	e = VAL4(A); /* C is 1+2^e x */
	NEG4(A);
	A[0]++; /* No overflow as A[0] is even. Now A = 1 - 2^e x */
	MUL4(D, C, A);
	SET4(C, D);
	e = 2*e;

	while (e < 256) {
		SET4(B, C); /* What is left over */
		B[0]--;
		NEG4(B);
		B[0]++; /* Now B = 1 - 2^{2e} y */ 
		MUL4(D, C, B);
		SET4(C, D); /* C is updated. */
		MUL4(D, A, B);
		SET4(A, D); /* A is updated. */
		e = 2*e;
	}
}
