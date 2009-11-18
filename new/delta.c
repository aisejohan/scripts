/*
 *	delta.c
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

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "helper.h"
#include "grobner.h"
#include "compute.h"
#include "delta.h"

/* Computes p*Delta.							*
 * This will only be run once!						*/
struct polynomial compute_delta(void)
{
	int i;
	struct polynomial A, B, C;
	A.leading = NULL;
	B.leading = NULL;
	C.leading = NULL;

	A = copy_pol(myf);
	B = copy_pol(myf);
	for (i = 2; i <= p; i++)
	{
		C = pol_mult(A, B);
		free_tail(B.leading);
		B = C;
		C.leading = NULL;
		C.degree = 0;
	}
	free_tail(A.leading);
	A.leading = NULL;
	A.degree = 0;
	A = frobenius(myf);

	/* Replace A by negative. */
	times_int(-1, &A);

	/* Add -F(f) + f^p */
	C = pol_add(A, B);

	free_tail(A.leading);
	free_tail(B.leading);
	
	return(C);
}

/* Scalar multiple of a split polynomial. */
struct polynomial **copy_pol_star(mscalar c, struct polynomial **bb)
{
	struct polynomial **uit;
	int i, len;
	len = 1 + bb[0]->degree/d;
	uit = (struct polynomial **)
		malloc(len*sizeof(struct polynomial *));
	if (!uit) {
		perror("Malloc failed! (Again?)");
		exit(1);
	}
	for (i = 0; i + 1 <= len; i++) {
		uit[i] = NULL;
		make_pol(&uit[i]);
		*uit[i] = copy_pol(*bb[i]);
		times_scalar(c, uit[i]);
	}
	return(uit);
}

/* Scalar multiple of a split polynomial. */
void free_star(struct polynomial **bb)
{
	int i, len;

	len = 1 + bb[0]->degree/d;
	for (i = 0; i + 1 <= len; i++) {
		free_tail(bb[i]->leading);
		free(bb[i]);
	}
	return;
}

/* Splits up a polynomial into pieces.				*
 * Removes the tail of f, and sets f.leading=NULL		*/
struct polynomial **split_up(struct polynomial *f)
{
	int i, count;
	struct polynomial **uit;
#ifdef OLD_GROBNER
	struct polynomial *tussen;
	struct polynomial **aa;
#else
	struct polynomial *aa;
#endif

	count = 1 + f->degree/d;
	uit = (struct polynomial **)malloc(count*sizeof(struct polynomial *));
	if (!uit) {
		perror("Malloc failed!");
		exit(1);
	}
	for (i = 0; i + 1 <= count; i++) {
		uit[i] = NULL;
		make_pol(&uit[i]);
		uit[i]->degree = f->degree - i*d;
	}

	/* Zero case still has to work.	*/
	if (f->leading == NULL) {
		return(uit);
	}

#ifdef OLD_GROBNER
	tussen = &myf;
#endif
	uit[0]->leading = f->leading;
	f->leading = NULL;
	for (i = 0; i + 1 + 1 <= count; i++) {
#ifdef OLD_GROBNER
		aa = gen_division(uit[i], 1, &tussen);
		times_int(-1, aa[0]);
		uit[i + 1]->leading = aa[0]->leading;
		uit[i + 1]->degree = aa[0]->degree;
		free(aa[0]);
#else
		aa = myf_division(uit[i]);
		times_int(-1, aa);
		uit[i + 1]->leading = aa->leading;
		uit[i + 1]->degree = aa->degree;
		free(aa);
#endif
	}
	return(uit);
}

/* Replaces f by f+g. Destroys the contents of g.
 * It could happen that the result has leading term 0. */
void merge_add_split(struct polynomial ***f, struct polynomial **g)
{
	int i, flen, glen;
	struct polynomial **tussen;

	flen = 1 + (*f)[0]->degree/d;
	glen = 1 + g[0]->degree/d;

	if (flen < glen) {
		tussen = *f;
		*f = g;
		g = tussen;
		i = glen;
		glen = flen;
		flen = i;
	}

	for (i = 0; i + 1 <= glen; i++) {
		merge_add((*f)[i + flen - glen], *g[i]);
		free(g[i]);
	}
	free(g);

	return;
}


/* Returns the product.						*
 * Does not modify f or g. 					*
 * The main term ``f[0]*g[0] mod myf'' is NOT zero since	*
 * neither f[0] nor g[0] is divisible by myf.			*/
struct polynomial **mult_split(struct polynomial **f, struct polynomial **g)
{
	int i, j, k, uitlen, flen, glen;
	struct polynomial tmp1, tmp2;
	struct polynomial **uit, **aa;
	tmp1.leading = NULL;
	tmp2.leading = NULL;

	flen = 1 + f[0]->degree/d;
	glen = 1 + g[0]->degree/d;
	uitlen = 1 + (f[0]->degree + g[0]->degree)/d;
	uit = (struct polynomial **)malloc(uitlen*sizeof(struct polynomial *));
	if (!uit) {
		perror("Malloc failed!");
		exit(1);
	}
	for (i = 0; i + 1 <= uitlen; i++) {
		uit[i] = NULL;
		make_pol(&uit[i]);
		uit[i]->degree = f[0]->degree + g[0]->degree - i*d; 
	}
	for (i = 0; i + 1 <= uitlen; i++) {
		tmp2.degree = uit[i]->degree;
		tmp2.leading = NULL;
		j= (glen < i + 1) ? (i + 1 - glen) : 0;
		while ((j <= i) && (j + 1 <= flen)) {
			tmp1 = pol_mult(*f[j], *g[i - j]);
			merge_add(&tmp2, tmp1);
			j++;
		}
		aa = split_up(&tmp2);
		/* The number of terms of aa will be uitlen-i.	*/
		/* We should also free aa again. 		*/
		for (k = 0; k + 1 <= uitlen - i; k++) {
			merge_add(uit[i + k], *aa[k]);
			free(aa[k]);
		}
		free(aa);
	}	
	return(uit);
}


#ifdef KIJKEN
/* Warning: destroys aa! */
void test_split(struct polynomial **aa, struct polynomial orig)
{
	int i, j, aalen;
	struct polynomial tmp;
	tmp.leading = NULL;
	
	if (orig.degree != aa[0]->degree) {
		printf("Wrong degrees! Stop.");
		exit(1);
	}
	
	aalen = 1 + aa[0]->degree/d;
	
	for (i = 1; i + 1 <= aalen; i++) {
		printf(" %d \n", i);
		tmp = pol_mult(myf, *aa[i]);
		free_tail(aa[i]->leading);
		merge_add(aa[0], tmp);
		for (j = i + 1; j + 1 <= aalen; j++) {
			tmp = pol_mult(myf, *aa[j]);
			free_tail(aa[j]->leading);
			*aa[j] = tmp;
		}
	}
	times_int(-1, aa[0]);
	printf("Here is the test result: \n");
	print_pol(pol_add(*aa[0], orig));
	return;
}
#endif
