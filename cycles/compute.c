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
static int s1, s2, s3, s4;
static int p12, p13, p14, p23, p24, p34;
static int t123, t124, t134, t234;

void allocate_GVMnew(int at_least)
{
	int old, i;

	if (G_len > at_least) return;

	if (G_len == 0) {
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

	/* Allocate memory for G */
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

void allocate_MMold(int at_least)
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
		free(G.ff[i]);
		free(G.ee[i]);
	}
	free(G.ee);
	free(G.ff);
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
	return(uit);
}

/* Rarely the case.							*/
static unsigned int rel_prime(struct exponents *mon1, struct exponents *mon2)
{
	if ((mon1->e1 > 0) && (mon2->e1 > 0)) return(0);
	if ((mon1->e2 > 0) && (mon2->e2 > 0)) return(0);
	if ((mon1->e3 > 0) && (mon2->e3 > 0)) return(0);
	if ((mon1->e4 > 0) && (mon2->e4 > 0)) return(0);
	return(1);
}

static unsigned int divides(struct exponents *mon1, struct exponents *mon2)
{
	return((mon1->e1 <= mon2->e1) && (mon1->e2 <= mon2->e2) && 
	(mon1->e3 <= mon2->e3) && (mon1->e4 <= mon2->e4));
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
	a->c = gterm->c;
	b->c = fterm->c;
	/* Note sign. */
	b->c = prime - b->c;
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


/* Outputs G.							*/
static unsigned int print_G(void)
{
	int i, s1 = 0, s2 = 0, s3 = 0, s4 = 0, success;
	struct exponents tmp;

	for (i = 0; i + 1 <= G.len; i++) {
		tmp = *G.ee[i];
		printf("[%d, %d, %d, %d]  \t",
			tmp.e1, tmp.e2, tmp.e3, tmp.e4);
		printf("%d\t", d1*tmp.e1 + d2*tmp.e2 + d3*tmp.e3 + d4*tmp.e4);
		printf("%d ", number_terms(*G.ff[i]));
		if (tmp.e1 + tmp.e2 + tmp.e3 == 0) {
			printf(" <--- 4");
			s4 = 1;
		}
		if (tmp.e1 + tmp.e2 + tmp.e4 == 0) {
			printf(" <--- 3");
			s3 = 1;
		}
		if (tmp.e1 + tmp.e3 + tmp.e4 == 0) {
			printf(" <--- 2");
			s2 = 1;
		}
		if (tmp.e2 + tmp.e3 + tmp.e4 == 0) {
			printf(" <--- 1");
			s1 = 1;
		}
		printf("\n");
	}
	success = s1 + s2 + s3 + s4;
	return(success);
}

static void adjust(struct exponents tmp)
{
	if (tmp.e1 + tmp.e2 + tmp.e3 == 0) t123 = 0;
	if (tmp.e1 + tmp.e2 + tmp.e4 == 0) t124 = 0;
	if (tmp.e1 + tmp.e3 + tmp.e4 == 0) t134 = 0;
	if (tmp.e2 + tmp.e3 + tmp.e4 == 0) t234 = 0;
	if (tmp.e1 + tmp.e2 == 0) p12 = 0;
	if (tmp.e1 + tmp.e2 == 0) p13 = 0;
	if (tmp.e1 + tmp.e2 == 0) p14 = 0;
	if (tmp.e2 + tmp.e3 == 0) p24 = 0;
	if (tmp.e2 + tmp.e4 == 0) p23 = 0;
	if (tmp.e3 + tmp.e4 == 0) p34 = 0;
	if (tmp.e1 == 0) s1 = 0;
	if (tmp.e2 == 0) s2 = 0;
	if (tmp.e3 == 0) s3 = 0;
	if (tmp.e4 == 0) s4 = 0;
}

static unsigned int compute_dim_G(void)
{
	if ((s1) || (s2) || (s3) || (s4)) return 3;
	if ((p12) || (p13) || (p14) || (p23) || (p24) || (p34)) return 2;
	if ((t123) || (t124) || (t134) || (t234)) return 1;
	return 0;
}

static unsigned int dim_G(void)
{
	int i;

	s1 = 1;
	s2 = 1;
	s3 = 1;
	s4 = 1;
	p12 = 1;
	p13 = 1;
	p14 = 1;
	p23 = 1;
	p24 = 1;
	p34 = 1;
	t123 = 1;
	t124 = 1;
	t134 = 1;
	t234 = 1;

	if (G.len == 0) return 4;
	for (i = 0; i + 1 <= G.len; i++) {
		adjust(*G.ee[i]);
	}
	return compute_dim_G();
}

/* Silly sort should be OK since the length of G is at most maxlength.	*
 * We sort the basis so that all the elements with high power of p	*
 * in the leading coefficient come last.				*/
static void sort_G(void)
{
	int i, j;
	struct exponents *s_ee;
	struct polynomial *s_ff;

	for (i = 0; i+1 <= G.len; i++) {
		for (j = i+1; j+1 <= G.len; j++) {
			if (!smaller(*G.ee[i],*G.ee[j])) {
					s_ee = G.ee[i];
					s_ff = G.ff[i];
					G.ee[i] = G.ee[j];
					G.ff[i] = G.ff[j];
					G.ee[j] = s_ee;
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

/* bound means dim >= bound */
int setup(struct polynomial A, struct polynomial B, struct polynomial C, struct polynomial D, int bound, int silent)
{
	int i, j, k, ii, jj, old, new, check, dimension;
	struct pair tmppair;
	struct polynomial SS, T;
	unsigned int m, mold, mnew;
	struct exponents lcm_new, lcm_old;
	SS.leading = NULL;
	T.leading = NULL;

	s1 = 1;
	s2 = 1;
	s3 = 1;
	s4 = 1;
	p12 = 1;
	p13 = 1;
	p14 = 1;
	p23 = 1;
	p24 = 1;
	p34 = 1;
	t123 = 1;
	t124 = 1;
	t134 = 1;
	t234 = 1;

	/* Allocate memory for G, V, M, Mold, Mnew */
	allocate_GVMnew(10);
	allocate_MMold(100);
	
	/* Initialize G */
	*G.ff[0] = copy_pol(A);
	*G.ee[0] = take_exponents(A);
	G.len = 1;
	adjust(*G.ee[0]);
	dimension = compute_dim_G();
	if (dimension < bound) goto uit;
	*G.ff[1] = copy_pol(B);
	*G.ee[1] = take_exponents(B);
	G.len = 2;
	adjust(*G.ee[1]);
	dimension = compute_dim_G();
	if (dimension < bound) goto uit;
	if (C.leading) {
		*G.ff[2] = copy_pol(C);
		*G.ee[2] = take_exponents(C);
		G.len = 3;
		adjust(*G.ee[2]);
		dimension = compute_dim_G();
		if (dimension < bound) goto uit;
		if (D.leading) {
			*G.ff[3] = copy_pol(D);
			*G.ee[3] = take_exponents(D);
			G.len = 4;
			adjust(*G.ee[3]);
			dimension = compute_dim_G();
			if (dimension < bound) goto uit;
		}
	}

	sort_G();

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
	while (m > 0) {
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
				check = 2; /* success. */
			} else {
				free_tail(SS.leading);
				check = 0;
			}
		} else {
			check = 0;
		}
	
		if (check == 2) {
			check = 0;
			/* Already increased G.len so have to substract one here. */
			gen_division(&SS, G.len - 1, G.ff);
			*G.ff[G.len - 1] = SS;
			*G.ee[G.len - 1] = take_exponents(SS); /* Done updating G. */

			adjust(*G.ee[G.len - 1]);
			dimension = compute_dim_G();
			if (dimension < bound) goto uit;

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
						(smaller(lcm_new,
							lcm(G.ee[Mnew[j].i],
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
				if (old < m) lcm_old = lcm(G.ee[Mold[old].i], G.ee[Mold[old].j]);
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
		}
	}/* End loop computing Grobner basis. */

	if (!silent) {
		sort_G();
		printf("The final length of G is %d\n", G.len);
		print_G();
		printf("------\n");
	}

	dimension = compute_dim_G();

uit:
	for (i = 0; i + 1 <= G.len; i++) {
		free_tail(G.ff[i]->leading);
	}

	return dimension;
}
