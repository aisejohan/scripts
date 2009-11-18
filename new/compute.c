/*
 *	compute.c
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
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "grobner.h"
#include "helper.h"
#include "compute.h"
#include "delta.h"
#include "reduce.h"

/* The rule is that 0<=i<j<=G.len always. 	*
 * Since maxlength is hopefully never bigger	*
 * than 2^16 short should be enough.		*/
struct pair {
	unsigned short int i;
	unsigned short int j;
};

/* Variable used outside this file as well. */
struct lijst G;
struct polynomial myf;

/* Variables only used in this file.			*/
static unsigned char **V;
static struct pair *M;
static struct pair *Mold;
static struct pair *Mnew;
static int M_len = 0;
static int G_len = 0;

static void allocate_GVMnew(int at_least)
{
	int old, i;

	if (G_len > at_least) return;

	if (G_len == 0) {
		G.BC = NULL;
		G.ff = NULL;
		G.ee = NULL;
		V = NULL;
		Mnew = NULL;
		old = 0;
		G_len = maxlength;
	} else {
		old = G_len;
		G_len = 2*G_len;
	}
	if (G_len < at_least + 2) G_len = at_least + 1;
	printf("Allocate G_len %d elements and old = %d.\n", G_len, old);

	/* Allocate memory for G */
	G.BC = (struct base_change **)
			realloc(G.BC, G_len*sizeof(struct base_change *));
	if (!G.BC) {
		perror("Malloc failed!");
		exit(1);
	}
	G.ff = (struct polynomial **)
			realloc(G.ff, G_len*sizeof(struct polynomial *));
	if (!G.ff) {
		perror("Malloc failed!");
		exit(1);
	}
	G.ee = (struct exponents **)
			realloc(G.ee, G_len*sizeof(struct exponents *));
	if (!G.ee) {
		perror("Malloc failed!");
		exit(1);
	}
	for (i = old; i < G_len; i++) {
		G.BC[i] = (struct base_change *)
				malloc(sizeof(struct base_change));
		if (!G.BC[i]) {
			perror("Malloc failed!");
			exit(1);
		}
		G.ff[i] = NULL;
		make_pol(&G.ff[i]);
		G.ee[i] = (struct exponents *)
				malloc(sizeof(struct exponents));
		if (!G.ee[i]) {
			perror("Malloc failed!");
			exit(1);
		}
	}
	V = (unsigned char **) realloc(V, G_len*sizeof(unsigned char *));
	if (!V) {
		perror("Malloc failed.\n");
		exit(1);
	}
	for (i = 0; i < old; i++) {
		V[i] = realloc(V[i], G_len*sizeof(unsigned char));
	}
	for (i = old; i < G_len; i++) {
		V[i] = malloc(G_len*sizeof(unsigned char));
		if (!V[i]) {
			perror("Malloc failed.\n");
			exit(1);
		}
	}
	Mnew = (struct pair *) realloc(Mnew, G_len*sizeof(struct pair));
}

static void allocate_MMold(int at_least)
{
	if (M_len > at_least) return;

	if (M_len == 0) {
		M_len = maxlength*maxlength;
		M = NULL;
		Mold = NULL;
	} else {
		M_len = 2*M_len;
	}
	if (M_len < at_least + 2) M_len = at_least + 1;
	printf("Allocate M_len = %d elements.\n", M_len);
	M = (struct pair *) realloc(M, M_len*sizeof(struct pair));
	Mold = (struct pair *) realloc(Mold, M_len*sizeof(struct pair));
}

void deallocate_GVMnewMMold(void )
{
	int i;

	free(Mnew);
	for (i = 0; i < G_len; i++) free(V[i]);
	free(V);
	for (i = 0; i < G_len; i++) {
		free(G.BC[i]);
		free(G.ff[i]);
		free(G.ee[i]);
	}
	free(G.ee);
	free(G.ff);
	free(G.BC);
	G.len = 0;
	free(Mold);
	free(M);
	M_len = 0;
}

/* Note that this produces a segfault or hangs if either	*
 * f.leading is NULL or if f.leading->c == 0.			*/
static struct exponents take_exponents(struct polynomial f)
{
	struct exponents uit;
	uit.e1 = f.leading->n1;
	uit.e2 = f.leading->n2;
	uit.e3 = f.leading->n3;
	uit.e4 = f.leading->n4;
	uit.e5 = (unsigned int) valuation(f.leading->c);
	return(uit);
}

/* Least common multiple.					*/
static struct exponents lcm(struct exponents *mon1, struct exponents *mon2)
{
	struct exponents uit;
	uit.e1 = (mon1->e1 > mon2->e1) ? mon1->e1 : mon2->e1;
	uit.e2 = (mon1->e2 > mon2->e2) ? mon1->e2 : mon2->e2;
	uit.e3 = (mon1->e3 > mon2->e3) ? mon1->e3 : mon2->e3;
	uit.e4 = (mon1->e4 > mon2->e4) ? mon1->e4 : mon2->e4;
	uit.e5 = (mon1->e5 > mon2->e5) ? mon1->e5 : mon2->e5;
	return(uit);
}

/* Rarely the case.							*/
static unsigned int rel_prime(struct exponents *mon1, struct exponents *mon2)
{
	if ((mon1->e1 > 0) && (mon2->e1 > 0)) return(0);
	if ((mon1->e2 > 0) && (mon2->e2 > 0)) return(0);
	if ((mon1->e3 > 0) && (mon2->e3 > 0)) return(0);
	if ((mon1->e4 > 0) && (mon2->e4 > 0)) return(0);
	if ((mon1->e5 > 0) && (mon2->e5 > 0)) return(0);
	return(1);
}

static unsigned int divides(struct exponents *mon1, struct exponents *mon2)
{
	return((mon1->e1 <= mon2->e1) && (mon1->e2 <= mon2->e2) && 
	(mon1->e3 <= mon2->e3) && (mon1->e4 <= mon2->e4) && 
	(mon1->e5 <= mon2->e5));
}

/* Smaller degree means smaller. Otherwise:				*
 * Make sure the ordering on the first 4 is the same as in the 		*
 * function kleiner, and finally if these are the same, then the	*
 * valuation of the coefficients being smaller means smaller.		*/
static unsigned int smaller(struct exponents mon1, struct exponents mon2)
{
	if (d1*mon1.e1 + d2*mon1.e2 + d3*mon1.e3 + d4*mon1.e4 !=
	d1*mon2.e1 + d2*mon2.e2 + d3*mon2.e3 + d4*mon2.e4) return((
		d1*mon1.e1 + d2*mon1.e2 + d3*mon1.e3 + d4*mon1.e4 < 
		d1*mon2.e1 + d2*mon2.e2 + d3*mon2.e3 + d4*mon2.e4));
	/* Same as in kleiner...				*/
#ifdef REVLEX_ORDER
	if (mon1.e4 != mon2.e4) return((mon1.e4 > mon2.e4));
	if (mon1.e3 != mon2.e3) return((mon1.e3 > mon2.e3));
	if (mon1.e2 != mon2.e2) return((mon1.e2 > mon2.e2));
#endif
#ifdef LEX_ORDER
	if (mon1.e1 != mon2.e1) return((mon1.e1 < mon2.e1));
	if (mon1.e2 != mon2.e2) return((mon1.e2 < mon2.e2));
	if (mon1.e3 != mon2.e3) return((mon1.e3 < mon2.e3));
#endif
	/* Extra measuring valuation. 				*/
	if (mon1.e5 != mon2.e5) return((mon1.e5 < mon2.e5));
	/* Means equal so not smaller. 				*/
	return(0);
}

/* Computes the coefficient terms needed to make the s_pol.	*/
static void s_pol_terms(struct term *a, struct term *b, struct term *fterm, struct term *gterm)
{
	if (fterm->n1 > gterm->n1) {
		a->n1 = 0;
		b->n1 = fterm->n1 - gterm->n1;
	} else {
		a->n1 = gterm->n1 - fterm->n1;
		b->n1 = 0;
	}
	if (fterm->n2 > gterm->n2) {
		a->n2 = 0;
		b->n2 = fterm->n2 - gterm->n2;
	} else {
		a->n2 = gterm->n2 - fterm->n2;
		b->n2 = 0;
	}
	if (fterm->n3 > gterm->n3) {
		a->n3 = 0;
		b->n3 = fterm->n3 - gterm->n3;
	} else {
		a->n3 = gterm->n3 - fterm->n3;
		b->n3 = 0;
	}
	if (fterm->n4 > gterm->n4) {
		a->n4 = 0;
		b->n4 = fterm->n4 - gterm->n4;
	} else {
		a->n4 = gterm->n4 - fterm->n4;
		b->n4 = 0;
	}
	sc_copy(gterm->c, a->c);
	sc_copy(fterm->c, b->c);
	/* Note sign. */
	sc_negate(b->c);
	while ((valuation(a->c) > 0) && (valuation(b->c) > 0)) {
		div_p(a->c, 1);
		div_p(b->c, 1);
	}
	return;
}

/* Computes the s_pol.						*/
static struct polynomial s_pol(struct polynomial f, struct polynomial g)
{
	struct term a, b;
	struct polynomial A, B;
	A.leading = NULL;
	B.leading = NULL;

	s_pol_terms(&a, &b, f.leading, g.leading);
	A = make_times_term(a, f);
	clean_pol(&A);
	B = make_times_term(b, g);
	merge_add(&A, B);
	return(A);
}

/* The copy paste method of programming.			*
 * Computes the base change vector of the s_pol.		*/
static struct base_change s_pol_BC(unsigned int i, unsigned int j)
{
	struct base_change uit;
	struct polynomial A, B;
	struct term a, b;
	A.leading = NULL;
	B.leading = NULL;

	s_pol_terms(&a, &b, G.ff[i]->leading, G.ff[j]->leading);
	/* Do the same onto BC as you do onto G.ff.	*/
	if ((G.BC[i]->bc1.leading) && (G.BC[j]->bc1.leading)) {
		A = make_times_term(a, G.BC[i]->bc1);
		clean_pol(&A);
		B = make_times_term(b, G.BC[j]->bc1);
		merge_add(&A, B);
		uit.bc1 = A;
	} else if (G.BC[i]->bc1.leading) {
		A = make_times_term(a, G.BC[i]->bc1);
		clean_pol(&A);
		uit.bc1 = A;
	} else if (G.BC[j]->bc1.leading) {
		B = make_times_term(b, G.BC[j]->bc1);
		clean_pol(&B);
		uit.bc1 = B;
	} else {
		uit.bc1.degree = 0; /* Not correct! */
		uit.bc1.leading = NULL;
	}
	if ((G.BC[i]->bc2.leading) && (G.BC[j]->bc2.leading)) {
		A = make_times_term(a, G.BC[i]->bc2);
		clean_pol(&A);
		B = make_times_term(b, G.BC[j]->bc2);
		merge_add(&A, B);
		uit.bc2 = A;
	} else if (G.BC[i]->bc2.leading) {
		A = make_times_term(a, G.BC[i]->bc2);
		clean_pol(&A);
		uit.bc2 = A;
	} else if (G.BC[j]->bc2.leading) {
		B = make_times_term(b, G.BC[j]->bc2);
		clean_pol(&B);
		uit.bc2 = B;
	} else {
		uit.bc2.degree = 0; /* Not correct! */
		uit.bc2.leading = NULL;
	}
	if ((G.BC[i]->bc3.leading) && (G.BC[j]->bc3.leading)) {
		A = make_times_term(a, G.BC[i]->bc3);
		clean_pol(&A);
		B = make_times_term(b, G.BC[j]->bc3);
		merge_add(&A, B);
		uit.bc3 = A;
	} else if (G.BC[i]->bc3.leading) {
		A = make_times_term(a, G.BC[i]->bc3);
		clean_pol(&A);
		uit.bc3 = A;
	} else if (G.BC[j]->bc3.leading) {
		B = make_times_term(b, G.BC[j]->bc3);
		clean_pol(&B);
		uit.bc3 = B;
	} else {
		uit.bc3.degree = 0; /* Not correct! */
		uit.bc3.leading = NULL;
	}
	if ((G.BC[i]->bc4.leading) && (G.BC[j]->bc4.leading)) {
		A = make_times_term(a, G.BC[i]->bc4);
		clean_pol(&A);
		B = make_times_term(b, G.BC[j]->bc4);
		merge_add(&A, B);
		uit.bc4 = A;
	} else if (G.BC[i]->bc4.leading) {
		A = make_times_term(a, G.BC[i]->bc4);
		clean_pol(&A);
		uit.bc4 = A;
	} else if (G.BC[j]->bc4.leading) {
		B = make_times_term(b, G.BC[j]->bc4);
		clean_pol(&B);
		uit.bc4 = B;
	} else {
		uit.bc4.degree = 0; /* Not correct! */
		uit.bc4.leading = NULL;
	}
	if ((G.BC[i]->bc5.leading) && (G.BC[j]->bc5.leading)) {
		A = make_times_term(a, G.BC[i]->bc5);
		clean_pol(&A);
		B = make_times_term(b, G.BC[j]->bc5);
		merge_add(&A, B);
		uit.bc5 = A;
	} else if (G.BC[i]->bc5.leading) {
		A = make_times_term(a, G.BC[i]->bc5);
		clean_pol(&A);
		uit.bc5 = A;
	} else if (G.BC[j]->bc5.leading) {
		B = make_times_term(b, G.BC[j]->bc5);
		clean_pol(&B);
		uit.bc5 = B;
	} else {
		uit.bc5.degree = 0; /* Not correct! */
		uit.bc5.leading = NULL;
	}
		
	return(uit);
}

#ifdef KIJKEN
/* Test function. 						*/
static void test_base_change(
	struct base_change B, struct polynomial new,
	struct polynomial myf1,
	struct polynomial myf2,
	struct polynomial myf3,
	struct polynomial myf4)
{
	unsigned int degree, i;
	struct polynomial lijst[10];
	degree = new.degree;
	for (i = 0; i <= 4; i++) {
		lijst[i].degree = degree; /* Just for this test function. */
		lijst[i].leading = NULL;
	}
	if (B.bc1.leading) lijst[0] = pol_mult(B.bc1, myf1);
	if (B.bc2.leading) lijst[1] = pol_mult(B.bc2, myf2);
	if (B.bc3.leading) lijst[2] = pol_mult(B.bc3, myf3);
	if (B.bc4.leading) lijst[3] = pol_mult(B.bc4, myf4);
	if (B.bc5.leading) lijst[4] = pol_mult(B.bc5, myf);
	lijst[5] = pol_add(lijst[0], lijst[1]);
	lijst[6] = pol_add(lijst[2], lijst[3]);
	lijst[7] = pol_add(lijst[5], lijst[6]);
	lijst[8] = pol_add(lijst[7], lijst[4]);
	times_int(-1, &lijst[8]);
	lijst[9] = pol_add(new, lijst[8]);
	if (lijst[9].leading) {
		print_pol(lijst[8]);
		printf("? =========================================== ?\n");
		print_pol(new);
		exit(1);
	}
	for (i = 0; i <= 9; i++) {
		free_tail(lijst[i].leading);
	}
	return;
}

/* Outputs M.							*/
static void print_M(unsigned int mm, struct pair *MM)
{
	int i;
	struct exponents tmp;
	for (i = 0; i + 1 <= mm; i++) {
		printf("[%d, %d] ", MM[i].i, MM[i].j);
		tmp = lcm(G.ee[MM[i].i], G.ee[MM[i].j]);
		printf("[%d, %d, %d, %d, %d] ", tmp.e1, tmp.e2, tmp.e3, tmp.e4,
					tmp.e5);
		printf("%d ", d1*tmp.e1 + d2*tmp.e2 + d3*tmp.e3 + d4*tmp.e4);
		printf("\n");
	}
	return;
}

/* Outputs V.							*/
static void print_V(unsigned int mm)
{
	int i, j;

	for (i = 0; i + 1 <= mm; i++) {
		for (j = 0; j + 1 <= mm; j++) {
			printf("%d ", V[i][j]);
		}
		printf("\n");
	}
	return;
}
#endif

/* Outputs G.							*/
static unsigned int print_G(void)
{
	int i, s1 = 0, s2 = 0, s3 = 0, s4 = 0, success;
	struct exponents tmp;

	for (i = 0; i + 1 <= G.len; i++) {
		tmp = *G.ee[i];
		printf("[%d, %d, %d, %d, %d]  \t",
			tmp.e1, tmp.e2, tmp.e3, tmp.e4, tmp.e5);
		printf("%d\t", d1*tmp.e1 + d2*tmp.e2 + d3*tmp.e3 + d4*tmp.e4);
		printf("%d ", number_terms(*G.ff[i]));
		if ((!tmp.e5) && ((tmp.e1 + tmp.e2 + tmp.e3) == 0)) {
			printf(" <--- 4");
			s1 = 1;
		}
		if ((!tmp.e5) && ((tmp.e1 + tmp.e2 + tmp.e4) == 0)) {
			printf(" <--- 3");
			s2 = 1;
		}
		if ((!tmp.e5) && ((tmp.e1 + tmp.e3 + tmp.e4) == 0)) {
			printf(" <--- 2");
			s3 = 1;
		}
		if ((!tmp.e5) && ((tmp.e2 + tmp.e3 + tmp.e4) == 0)) {
			printf(" <--- 1");
			s4 = 1;
		}
#ifdef KIJKEN
		if (
		(tmp.e1 != G.ff[i]->leading->n1) ||
		(tmp.e2 != G.ff[i]->leading->n2) ||
		(tmp.e3 != G.ff[i]->leading->n3) ||
		(tmp.e4 != G.ff[i]->leading->n4) ||
		(tmp.e5 != valuation(G.ff[i]->leading->c))) {
			printf("Wrong exponents!\n");
			exit(1);
		}
#endif
		printf("\n");
	}
	success = s1 + s2 + s3 + s4;
	return(success);
}

static unsigned int test_G(void)
{
	int i, s1 = 0, s2 = 0, s3 = 0, s4 = 0, success;
	struct exponents tmp;

	for (i = 0; i + 1 <= G.len; i++) {
		tmp = *G.ee[i];
		if ((!tmp.e5) && ((tmp.e1 + tmp.e2 + tmp.e3) == 0)) {
			s1 = 1;
		}
		if ((!tmp.e5) && ((tmp.e1 + tmp.e2 + tmp.e4) == 0)) {
			s2 = 1;
		}
		if ((!tmp.e5) && ((tmp.e1 + tmp.e3 + tmp.e4) == 0)) {
			s3 = 1;
		}
		if ((!tmp.e5) && ((tmp.e2 + tmp.e3 + tmp.e4) == 0)) {
			s4 = 1;
		}
#ifdef KIJKEN
		if (
		(tmp.e1 != G.ff[i]->leading->n1) ||
		(tmp.e2 != G.ff[i]->leading->n2) || 
		(tmp.e3 != G.ff[i]->leading->n3) || 
		(tmp.e4 != G.ff[i]->leading->n4) || 
		(tmp.e5 != valuation(G.ff[i]->leading->c))) {
			printf("Wrong exponents!\n");
			exit(1);
		}
#endif
	}
	success = s1 + s2 + s3 + s4;
	return(success);
}

/* Silly sort should be OK since the length of G is at most maxlength.	*
 * We sort the basis so that all the elements with high power of p	*
 * in the leading coefficient come last.				*/
static void sort_G(void)
{
	int i, j;
	struct exponents *s_ee;
	struct base_change *s_bc;
	struct polynomial *s_ff;

	for (i = 0; i+1 <= G.len; i++) {
		for (j = i+1; j+1 <= G.len; j++) {
			if (
			/* Test for p-adic valuation of leading coefficient. */
			(G.ee[j]->e5 < G.ee[i]->e5) ||
			/* Test for degree of leading term. */
			((G.ee[j]->e5 == G.ee[i]->e5) &&
				(G.ff[j]->degree < G.ff[i]->degree)) ||
			/* Test for ordering. */
			((G.ee[j]->e5 == G.ee[i]->e5) &&
				(G.ff[j]->degree == G.ff[i]->degree) &&
				((smaller(*G.ee[i],*G.ee[j]))))) {
					s_ee = G.ee[i];
					s_bc = G.BC[i];
					s_ff = G.ff[i];
					G.ee[i] = G.ee[j];
					G.BC[i] = G.BC[j];
					G.ff[i] = G.ff[j];
					G.ee[j] = s_ee;
					G.BC[j] = s_bc;
					G.ff[j] = s_ff;
			}
		}
	}
}


static unsigned int test_skip(struct pair try, struct exponents least)
{
	int k;

	for (k = 0; k + 1 <= try.i; k++) {
		if ((!V[k][try.i]) && (!V[k][try.j]) &&
					divides(G.ee[k],&least)) {
			return(1);
		}
	}
	for (k = try.i + 1; k + 1 <= try.j; k++) {
		if ((!V[try.i][k]) && (!V[k][try.j]) &&
					divides(G.ee[k],&least)) {
			return(1);
		}
	}
	for (k = try.j + 1; k + 1 <= G.len; k++) {
		if ((!V[try.i][k]) && (!V[try.j][k]) &&
					divides(G.ee[k],&least)) {
			return(1);
		}
	}
	return(0);
}

int setup(int silent)
{
	int i, j, k, ii, jj, old, new, check, epsilon;
	mscalar c;
	struct pair tmppair;
	struct polynomial EEN, NIKS, SS, T;
	struct polynomial *Tff;
	struct exponents *Tee;
	struct base_change *TBC;
	struct polynomial **aa, **bb;
	struct polynomial myf1, myf2, myf3, myf4;
	unsigned int m, mold, mnew;
	struct exponents lcm_new, lcm_old;
	EEN.leading = NULL;
	NIKS.leading = NULL;
	SS.leading = NULL;
	T.leading = NULL;

	/* Unit polynomial */
	EEN.degree = 0;
	make_term(&EEN.leading);
	sc_one(EEN.leading->c);
	EEN.leading->n1 = 0;
	EEN.leading->n2 = 0;
	EEN.leading->n3 = 0;
	EEN.leading->n4 = 0;
	
	/* Zero polynomial of degree 0. */
	NIKS.degree = 0;

	/* Initialize myf,myf1,myf2,myf3,myf4 */
	if (!silent) {
		printf("\n");
		printf("Here is the polynomial we're using this time:\n");
		printf("\n");
		print_pol(myf);
	}
	if (!myf.leading) {
		if (!silent) printf("Polynomial is zero!\n");
		free_tail(EEN.leading);
		return(1);
	}
	if (!silent) printf("\n");
	
	myf1 = deriv(myf, 1);
	if (!myf1.leading) {
		if (!silent) printf("Polynomial does not depend on x!\n");
		free_tail(EEN.leading);
		free_tail(myf.leading);
		return(1);
	}
	myf2 = deriv(myf, 2);
	if (!myf2.leading) {
		if (!silent) printf("Polynomial does not depend on y!\n");
		free_tail(EEN.leading);
		free_tail(myf.leading);
		free_tail(myf1.leading);
		return(1);
	}
	myf3 = deriv(myf, 3);
	if (!myf3.leading) {
		if (!silent) printf("Polynomial does not depend on z!\n");
		free_tail(EEN.leading);
		free_tail(myf.leading);
		free_tail(myf1.leading);
		free_tail(myf2.leading);
		return(1);
	}
	myf4 = deriv(myf, 4);
	if (!myf4.leading) {
		if (!silent) printf("Polynomial does not depend on w!\n");
		free_tail(EEN.leading);
		free_tail(myf.leading);
		free_tail(myf1.leading);
		free_tail(myf2.leading);
		free_tail(myf3.leading);
		return(1);
	}
	if (!silent) printf("\n");
	
	/* Allocate memory for G, V, M, Mold, Mnew */
	allocate_GVMnew(10);
	allocate_MMold(100);
	
	/* Initialize G */
	*G.ff[0] = copy_pol(myf4);
	G.BC[0]->bc1 = copy_pol(NIKS);
	G.BC[0]->bc2 = copy_pol(NIKS);
	G.BC[0]->bc3 = copy_pol(NIKS);
	G.BC[0]->bc4 = copy_pol(EEN);
	G.BC[0]->bc5 = copy_pol(NIKS);
	*G.ee[0] = take_exponents(myf4);

	*G.ff[1] = copy_pol(myf3);
	G.BC[1]->bc1 = copy_pol(NIKS);
	G.BC[1]->bc2 = copy_pol(NIKS);
	G.BC[1]->bc3 = copy_pol(EEN);
	G.BC[1]->bc4 = copy_pol(NIKS);
	G.BC[1]->bc5 = copy_pol(NIKS);
	*G.ee[1] = take_exponents(myf3);

	*G.ff[2] = copy_pol(myf2);
	G.BC[2]->bc1 = copy_pol(NIKS);
	G.BC[2]->bc2 = copy_pol(EEN);
	G.BC[2]->bc3 = copy_pol(NIKS);
	G.BC[2]->bc4 = copy_pol(NIKS);
	G.BC[2]->bc5 = copy_pol(NIKS);
	*G.ee[2] = take_exponents(myf2);

	*G.ff[3] = copy_pol(myf1);
	G.BC[3]->bc1 = copy_pol(EEN);
	G.BC[3]->bc2 = copy_pol(NIKS);
	G.BC[3]->bc3 = copy_pol(NIKS);
	G.BC[3]->bc4 = copy_pol(NIKS);
	G.BC[3]->bc5 = copy_pol(NIKS);
	*G.ee[3] = take_exponents(myf1);

	*G.ff[4] = copy_pol(myf);
	G.BC[4]->bc1 = copy_pol(NIKS);
	G.BC[4]->bc2 = copy_pol(NIKS);
	G.BC[4]->bc3 = copy_pol(NIKS);
	G.BC[4]->bc4 = copy_pol(NIKS);
	G.BC[4]->bc5 = copy_pol(EEN);
	*G.ee[4] = take_exponents(myf);

	G.len = 5;

	/* Deal with leading coefficients being divisible by p! */
	i = 0;
	while (i + 1 <= G.len) {
		if (G.ee[i]->e5 > 0) {
			T = copy_pol(*G.ff[i]);
			sc_one(c);
			for (j = 1; j <= r-G.ee[i]->e5; j++) 
					sc_imult_replace(p, c);
			times_scalar(c, &T);
			if (T.leading) {
				G.len++;
				*G.ff[G.len - 1] = T;
				*G.ee[G.len - 1] = take_exponents(T);
				T = copy_pol(G.BC[i]->bc1);
				times_scalar(c, &T);
				G.BC[G.len - 1]->bc1 = T;
				T = copy_pol(G.BC[i]->bc2);
				times_scalar(c, &T);
				G.BC[G.len - 1]->bc2 = T;
				T = copy_pol(G.BC[i]->bc3);
				times_scalar(c, &T);
				G.BC[G.len - 1]->bc3 = T;
				T = copy_pol(G.BC[i]->bc4);
				times_scalar(c, &T);
				G.BC[G.len - 1]->bc4 = T;
				T = copy_pol(G.BC[i]->bc5);
				times_scalar(c, &T);
				G.BC[G.len-1]->bc5 = T;
			}
		}
		i++;
	}
	check = 0;
	
	/* Initialize V */
	for (i = 0; i + 1 <= G.len; i++) {
		for (j = 0; j + 1 <= G.len; j++) {
			V[i][j] = 0;
		}
	}

	/* Initialize M and m. */
	m = 0;
	for (i = 0; i + 1 <= G.len; i++) {
		for (j = i + 1; j + 1 <= G.len; j++) {
			if (!rel_prime(G.ee[i], G.ee[j])) {
				m = m + 1;
				M[m - 1].i = i;
				M[m - 1].j = j;
				V[i][j] = 1;
			}
		}
	}

	/* Order entries in M such that the smallest one comes last! */
	for (i = 0; i + 1 <= m; i++) {
		for (j = i + 1; j + 1 <= m; j++) {
			if (smaller(lcm(G.ee[M[i].i], G.ee[M[i].j]),
					lcm(G.ee[M[j].i], G.ee[M[j].j]))) {
				tmppair = M[i];
				M[i] = M[j];
				M[j] = tmppair;
			}
		}
	}


/* Loop for computing the Grobner basis.	*
 * Needlessly complicated.			*/
while ((m > 0) || (check == 1)) {
	if (check) {
		/* Here we multiply the previous one found by a power	*
		 * of p to cancel off the leading term and we see if	*
		 * there is anything left.				*/
		sc_one(c);
		for (i = 1; i <= r-G.ee[G.len - 1]->e5; i++)
			sc_imult_replace(p, c);
		SS = copy_pol(*G.ff[G.len - 1]);
		/* Multiply by power of p. */
		times_scalar(c, &SS);
		if ((SS.leading) && (!zero_on_division(SS, G.len, G.ff))) {
			G.len++;
			if (G.len > G_len) {
				allocate_GVMnew(G.len);
			}
			T = copy_pol(G.BC[G.len - 1 - 1]->bc1);
			times_scalar(c, &T);
			G.BC[G.len - 1]->bc1 = T;
			T = copy_pol(G.BC[G.len - 1 - 1]->bc2);
			times_scalar(c, &T);
			G.BC[G.len - 1]->bc2 = T;
			T = copy_pol(G.BC[G.len - 1 - 1]->bc3);
			times_scalar(c, &T);
			G.BC[G.len - 1]->bc3 = T;
			T = copy_pol(G.BC[G.len - 1 - 1]->bc4);
			times_scalar(c, &T);
			G.BC[G.len - 1]->bc4 = T;
			T = copy_pol(G.BC[G.len - 1 - 1]->bc5);
			times_scalar(c, &T);
			G.BC[G.len - 1]->bc5 = T;
			check = 2; /* success. */
		} else {
			free_tail(SS.leading);
			check = 0; /* go around. */
		}
	} else {
		/* Here we take the last pair from M and we reduce	*
		 * it and we see if there is anything left.		*/
		ii = M[m - 1].i;
		jj = M[m - 1].j;
		V[ii][jj] = 0; 		/* Update V. */
		m = m - 1;		/* Update M. */
		if (!test_skip(M[m], lcm(G.ee[ii], G.ee[jj]))) {
			/* Make S-pol. */
			SS = s_pol(*G.ff[ii], *G.ff[jj]); 
			if ((SS.leading) &&
					(!zero_on_division(SS, G.len, G.ff))) {
				G.len++;		
				if (G.len > G_len) {
					allocate_GVMnew(G.len);
				}
				*G.BC[G.len - 1] = s_pol_BC(ii, jj);
				check = 2; /* success. */
			} else {
				free_tail(SS.leading);
				check = 0;
			}
		} else {
			check = 0;
		}
	}
	
	/* Common code for the two cases of adding a member of G. */
	if (check == 2) {
		/* Already increased G.len so have to substract one here. */
		aa = gen_division(&SS, G.len - 1, G.ff);
		*G.ff[G.len - 1] = SS;
		for (i = 0; i + 1 + 1 <= G.len; i++)  {
			if ((G.BC[i]->bc1.leading) && 
					(G.BC[G.len - 1]->bc1.leading)) {
				T = pol_mult(*aa[i], G.BC[i]->bc1);
				merge_add(&(G.BC[G.len-1]->bc1), T);
			} else if (G.BC[i]->bc1.leading) {
				G.BC[G.len - 1]->bc1 = 
					pol_mult(*aa[i], G.BC[i]->bc1);
			}
			if ((G.BC[i]->bc2.leading) && 
					(G.BC[G.len-1]->bc2.leading)) {
				T = pol_mult(*aa[i], G.BC[i]->bc2);
				merge_add(&(G.BC[G.len - 1]->bc2), T);
			} else if (G.BC[i]->bc2.leading) {
				G.BC[G.len - 1]->bc2 = 
					pol_mult(*aa[i], G.BC[i]->bc2);
			}
			if ((G.BC[i]->bc3.leading) && 
					(G.BC[G.len - 1]->bc3.leading)) {
				T = pol_mult(*aa[i], G.BC[i]->bc3);
				merge_add(&(G.BC[G.len - 1]->bc3), T);
			} else if (G.BC[i]->bc3.leading) {
				G.BC[G.len - 1]->bc3 = 
					pol_mult(*aa[i], G.BC[i]->bc3);
			}
			if ((G.BC[i]->bc4.leading) && 
					(G.BC[G.len - 1]->bc4.leading)) {
				T = pol_mult(*aa[i], G.BC[i]->bc4);
				merge_add(&(G.BC[G.len - 1]->bc4), T);
			} else if (G.BC[i]->bc4.leading) {
				G.BC[G.len - 1]->bc4 = 
					pol_mult(*aa[i], G.BC[i]->bc4);
			}
			if ((G.BC[i]->bc5.leading) && 
					(G.BC[G.len - 1]->bc5.leading)) {
				T = pol_mult(*aa[i], G.BC[i]->bc5);
				merge_add(&(G.BC[G.len - 1]->bc5), T);
			} else if (G.BC[i]->bc5.leading) {
				G.BC[G.len - 1]->bc5 = 
					pol_mult(*aa[i], G.BC[i]->bc5);
			}
		}
		*G.ee[G.len - 1] = take_exponents(SS); /* Done updating G. */

#ifdef KIJKEN
		test_base_change(*G.BC[G.len - 1], SS, myf1, myf2, myf3, myf4);
#endif

		/* Frees space allocated for aa. */
		for (i = 0; i + 1 + 1 <= G.len; i++) {
			free_tail(aa[i]->leading);
			free(aa[i]);
		}
		free(aa); 

		/* Update M and V. */
		for (i = 0; i < G.len; i++) {
			V[G.len - 1][i] = 0;
			V[i][G.len - 1] = 0;
		}

		/* List the new pairs in order in Mnew. */
		mnew = 0;
		for (i = 0; i<= (G.len - 1) - 1; i++) {
			if (!rel_prime(G.ee[i], G.ee[G.len - 1])) {
				lcm_new = lcm(G.ee[i], G.ee[G.len - 1]);
				j = 0;
				while ((j + 1 <= mnew) &&
					(smaller(lcm_new, lcm(G.ee[Mnew[j].i],
						G.ee[Mnew[j].j])))) j++;
				if (j == mnew) {
					mnew = mnew + 1;
					Mnew[mnew - 1].i = i;
					Mnew[mnew - 1].j = G.len - 1;
					V[i][G.len - 1] = 1;
				} else {
					for (k = mnew; k >= j + 1; k--) {
						Mnew[k] = Mnew[k - 1];
					}
					mnew = mnew + 1;
					Mnew[j].i = i;
					Mnew[j].j = G.len - 1;
					V[i][G.len - 1] = 1;
				}
			}
		}

		if (mnew > 0) {
			/* Save the M we have sofar into Mold. */
			for (i = 0; i + 1 <= m; i++) {
				Mold[i] = M[i];
			}
			mold = m;

			/* Merge old and new into M. */
			old = 0;
			lcm_old = lcm(G.ee[Mold[old].i], G.ee[Mold[old].j]);
			new = 0;
			lcm_new = lcm(G.ee[Mnew[new].i], G.ee[Mnew[new].j]);
			m = mold + mnew;
			allocate_MMold(m);
			i = 0;
			while ((new + 1 <= mnew) && (old + 1 <= mold)) {
				if (smaller(lcm_new, lcm_old)) {
					M[i] = Mold[old];
					i = i + 1;
					old = old + 1;
					if (old + 1 <= mold) lcm_old =
						lcm(G.ee[Mold[old].i],
						G.ee[Mold[old].j]);
				} else {
					M[i] = Mnew[new];
					i = i + 1;
					new = new + 1;
					if (new + 1 <= mnew) lcm_new = 
						lcm(G.ee[Mnew[new].i],
						G.ee[Mnew[new].j]);
				}
			}
			while (old + 1 <= mold) {
				M[i] = Mold[old];
				i = i + 1;
				old = old + 1;
			}
			while (new + 1 <= mnew) {
				M[i] = Mnew[new];
				i = i + 1;
				new = new + 1;
			}
		}

		/* Do it again for the powers of p times the new one!	*
		 * This is needed to deal with something like		*
		 * p*x + y.						*/
		if (G.ee[G.len - 1]->e5 > 0) {
			check = 1;
		} else {
			check = 0;
		}
	}
} /* End loop computing Grobner basis. */

	i = test_G();
	if (i < 4) {
		/* These are safe to free. */
		free_tail(EEN.leading);
		free_tail(myf.leading);
		free_tail(myf1.leading);
		free_tail(myf2.leading);
		free_tail(myf3.leading);
		free_tail(myf4.leading);
		for (i = 0; i + 1 <= G.len; i++) {
			free_tail(G.BC[i]->bc1.leading);
			free_tail(G.BC[i]->bc2.leading);
			free_tail(G.BC[i]->bc3.leading);
			free_tail(G.BC[i]->bc4.leading);
			free_tail(G.BC[i]->bc5.leading);
			free_tail(G.ff[i]->leading);
		}
		if (!silent) printf("Not smooth!\n");
		return(1);
	}

#ifdef KIJKEN
	printf("The initial length of G is %d.\n", G.len);
	print_G();
	printf("\n");
#endif
	
	/* Weed out G. 						*
	 * Do not need to update M,m,V since they are no	*
	 * longer used.						*
	 * Below we order the elements of G.			*/
	old = G.len; /* Remember for freeing bb later. */
	bb = (struct polynomial **)malloc(G.len*sizeof(struct polynomial *));
	if (!bb) {
		perror("Malloc failed!");
		exit(1);
	}
	for (i = 0; i + 1 <= G.len; i++) {
		bb[i] = NULL;
		make_pol(&bb[i]);
	}
	i = 0;
	while (i + 1 <= G.len) {
		for (j = 0; j + 1 <= G.len; j++) {
			if (j != i) {
				epsilon = (j>i) ? 1 : 0;
				bb[j-epsilon]->degree = G.ff[j]->degree;
				bb[j-epsilon]->leading = G.ff[j]->leading;
			}
		}
		
		new = G.len - 1; /* Remember for freeing aa later. */
		aa = gen_division(G.ff[i], G.len - 1, bb);

#ifdef KIJKEN
		if ((G.ff[i]->leading) && 
		((G.ff[i]->leading->n1 != G.ee[i]->e1) || 
		 (G.ff[i]->leading->n2 != G.ee[i]->e2) || 
		 (G.ff[i]->leading->n3 != G.ee[i]->e3) || 
		 (G.ff[i]->leading->n4 != G.ee[i]->e4))) {
			printf("The following should have been zero: ");
			print_pol(*G.ff[i]);
			exit(1);
		}
#endif
		
		/* Either omit G[i] or replace it. */
		if (!G.ff[i]->leading) {
			Tff = G.ff[i];
			Tee = G.ee[i];
			TBC = G.BC[i];
			for (j = i; j + 1 + 1 <= G.len; j++) {
				G.ff[j] = G.ff[j + 1];
				G.ee[j] = G.ee[j + 1];
				G.BC[j] = G.BC[j + 1];
			}
			G.ff[G.len - 1] = Tff;
			G.ee[G.len - 1] = Tee;
			free_tail(TBC->bc1.leading);
			free_tail(TBC->bc2.leading);
			free_tail(TBC->bc3.leading);
			free_tail(TBC->bc4.leading);
			free_tail(TBC->bc5.leading);
			G.BC[G.len - 1] = TBC;
			G.len--;
		} else {
			/* Note that we do not have to update G.len	*
			 * or G.ff[i] since it has already been		*
			 * changed.					*/
			for (j = 0; j + 1 <= G.len; j++)  {
				/* Do not allow j=i, and don't do it
				 * when the coefficient is zero. */
				epsilon = (j > i) ? 1 : 0;
				if ((j != i) && (aa[j - epsilon]->leading)) {
					if ((G.BC[j]->bc1.leading) &&
					(G.BC[i]->bc1.leading)) {
						T = pol_mult(*aa[j - epsilon],
							G.BC[j]->bc1);
						merge_add(&(G.BC[i]->bc1), T);
					} else if (G.BC[j]->bc1.leading) {
						G.BC[i]->bc1 = pol_mult(
							*aa[j - epsilon],
							G.BC[j]->bc1);
					}
					if ((G.BC[j]->bc2.leading) &&
					(G.BC[i]->bc2.leading)) {
						T = pol_mult(*aa[j - epsilon],
							G.BC[j]->bc2);
						merge_add(&(G.BC[i]->bc2), T);
					} else if (G.BC[j]->bc2.leading) {
						G.BC[i]->bc2 = pol_mult(
							*aa[j - epsilon],
							G.BC[j]->bc2);
					}
					if ((G.BC[j]->bc3.leading) &&
					(G.BC[i]->bc3.leading)) {
						T = pol_mult(*aa[j - epsilon],
							G.BC[j]->bc3);
						merge_add(&(G.BC[i]->bc3), T);
					} else if (G.BC[j]->bc3.leading) {
						G.BC[i]->bc3 = pol_mult(
							*aa[j - epsilon],
							G.BC[j]->bc3);
					}
					if ((G.BC[j]->bc4.leading) &&
						(G.BC[i]->bc4.leading)) {
						T = pol_mult(*aa[j - epsilon],
							G.BC[j]->bc4);
						merge_add(&(G.BC[i]->bc4), T);
					} else if (G.BC[j]->bc4.leading) {
						G.BC[i]->bc4 = pol_mult(
							*aa[j - epsilon],
							G.BC[j]->bc4);
					}
					if ((G.BC[j]->bc5.leading) &&
					(G.BC[i]->bc5.leading)) {
						T = pol_mult(*aa[j - epsilon],
							G.BC[j]->bc5);
						merge_add(&(G.BC[i]->bc5), T);
					} else if (G.BC[j]->bc5.leading) {
						G.BC[i]->bc5 = pol_mult(
							*aa[j - epsilon],
							G.BC[j]->bc5);
					}
				}
			}

#ifdef KIJKEN
			/* Test.					*/
			test_base_change(*G.BC[i], *G.ff[i],
						myf1, myf2, myf3, myf4);
#endif

			/* This should not be necessary. */
			*G.ee[i] = take_exponents(*G.ff[i]); 
			/* Done updating G. */

			/* Only in this case do we update i! */
			i++;
		}

		/* Free aa. */
		for (j = 0; j + 1 <= new; j++) {
			free_tail(aa[j]->leading);
			free(aa[j]);
		}
		free(aa);
	}
	
	/* Free bb. */
	for (j = 0; j + 1 <= old; j++) {
		bb[j]->leading = NULL;
		free(bb[j]);
	}
	free(bb);

	sort_G();

	if (!silent) {
		printf("The final length of G is %d\n", G.len);
		print_G();
		printf("------\n");
	}

#ifdef KIJKEN
	/* Sanity Check! Takes some time... */
	printf("Checking p-powers!\n");
	for (i = 0; i <= G.len - 1; i++) {
		if (G.ee[i]->e5 > 0) {
			sc_one(c);
			for (j = 1; j <= r - G.ee[i]->e5; j++)
				sc_imult_replace(p, c);
			T = copy_pol(*G.ff[i]);
			times_scalar(c, &T);
			bb = gen_division(&T, G.len, G.ff);
			if (T.leading) {
				printf("NOT GOOD.\n");
				print_pol(T);
				printf("\n");
				exit(1);
			}
			for (j = 0; j + 1 <= G.len; j++) {
				free_tail(bb[j]->leading);
				free(bb[j]);
			}
			free(bb);
		}
	}
	/* Recheck all S-pols reduce to zero! */
	printf("Checking S-pols.\n");
	for (i = 0; i <= G.len - 1; i++) {
		for (j = i + 1; j <= G.len - 1; j++) {
			SS = s_pol(*G.ff[i], *G.ff[j]);
			if (zero_on_division(SS, G.len, G.ff)) {
				free_tail(SS.leading);
			} else {
				printf("Not OK at the following"
				" indices: %d and %d\n", i, j);
				print_pol(*G.ff[i]);
				printf("\n");
				print_pol(*G.ff[j]);
				printf("\n");
				print_pol(SS);
				exit(1);
			}
		}
	}
#endif /* KIJKEN */

	/* These are not used outside this file... */
	free_tail(EEN.leading);
	free_tail(myf1.leading);
	free_tail(myf2.leading);
	free_tail(myf3.leading);
	free_tail(myf4.leading);

	/* Success. */
	return(0);
}
