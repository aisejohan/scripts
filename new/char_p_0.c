/*
 *	char_p_0.c
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
#include "reduce.h"

int *find_gap(void )
{
	int i, ma, a[2], mb, *b;
	int powers[r];
	
	for (i = 0; i + 1 <= r; i++) powers[i] = 0;
	for (i = 0; i + 1 <= G.len; i++) powers[G.ee[i]->e5] = 1;
	ma = 0;
	a[0] = 0;
	a[1] = 0;
	mb = 0;
	b = (int *) malloc(2*sizeof(int));
	if (!b) {
		printf("Malloc failed!");
		exit(1);
	}
	b[0] = 0;
	b[1] = 0;
	for (i = 0; i + 2 <= r; i++) {
		if ((powers[i] == 1) && (powers[i + 1] == 0)) {
			a[0] = i + 1;
			ma = 0;
		}
		if (powers[i] == 0) ma++;
		if ((powers[i] == 0) && (powers[i + 1] == 1)) a[1] = i;
		if ((i + 2 == r) && (powers[i] == 0) && (powers[i + 1] == 0)) {
			a[1] = i + 1;
			ma++;
		}
		if ((a[0] < a[1]) && (ma > mb)) {
			mb = ma;
			b[0] = a[0];
			b[1] = a[1];
		}
	}
	if (b[0] == 0) {
		printf("No gap found! You need to increase the integer r.\n");
		exit(1);
	}
	printf("The gap is between %d and %d.\n", b[0], b[1]);
	return(b);
}

int char_0(unsigned int degree, int *gap)
{
	int i, test;
	int count;
	unsigned int a1, a2, a3, a4;

	count = 0;

	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
		test = 0;
		for (i = 0; i + 1 <= G.len; i++) {
			if ((G.ee[i]->e1 <= a1) &&
			(G.ee[i]->e2 <= a2) &&
			(G.ee[i]->e3 <= a3) &&
			(G.ee[i]->e4 <= a4) &&
			(G.ee[i]->e5 < gap[0])) test = 1;
		}
		if (!test) count++;
	      }
	    }
	  }
	}
	return(count);
}

/* Order the list so the largest is first.	*
 * This is stupid sorting so hopefully		*
 * the list is not too long!			*/
static void sort_terms(struct term **tt, int blen)
{
	int i, j;
	struct term tmp;

	for (i = 0; i < blen; i++) {
		for (j = i + 1; j < blen; j++) {
			if (kleiner(tt[i], tt[j]) == KLEINER) {
				copy_term(tt[i], &tmp);
				copy_term(tt[j], tt[i]);
				copy_term(&tmp, tt[j]);
			}
		}
	}
}

void print_terms(struct term **tt, int blen)
{
	int i;
	struct polynomial T;

	for (i = 0; i < blen; i++) {
		T.degree = tt[i]->n1*d1 + tt[i]->n2*d2 + tt[i]->n3*d3
								+ tt[i]->n4*d4;
		T.leading = tt[i];
		print_pol(T);
		T.leading = NULL;
	}
	printf("\n");
}

void free_list_terms(struct term **tt, int blen)
{
	int i;

	for (i = 0; i < blen; i++) {
		free_term(tt[i]);
	}
	free(tt);
}

/* Finds the char 0 basis of terms in degree degree.		*
 * This function assumes the function char_0 has been run	*
 * previously and has returned the nonnegative integer blen.	*/
struct term **char_0_basis(unsigned int degree, int blen, int *gap)
{
	int i, a1, a2, a3, a4, count, test;
	struct term tmp;
	struct term **tt;

	tt = (struct term **)malloc(blen*sizeof(struct term *));
	if (!tt) {
		perror("Malloc failed!");
		exit(1);
	}
	for (i = 0; i + 1 <= blen; i++) {
		tt[i] = NULL;
		make_term(&tt[i]);
	}
	
	count = 0;
	sc_one(tmp.c);
	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
		test = 0;
		for (i = 0; i + 1 <= G.len; i++) {
			if ((G.ee[i]->e1 <= a1) &&
			(G.ee[i]->e2 <= a2) &&
			(G.ee[i]->e3 <= a3) &&
			(G.ee[i]->e4 <= a4) &&
			(G.ee[i]->e5 < gap[0])) test = 1;
		}
		if (!test) {
			count++;
			if (count > blen) {
				printf("Wrong length basis!");
				exit(1);
			}
			/* tmp.c = 1 */
			tmp.n1 = a1;
			tmp.n2 = a2;
			tmp.n3 = a3;
			tmp.n4 = a4;
			copy_term(&tmp, tt[count - 1]);
			tt[count - 1]->next = NULL;
		}
	      }
	    }
	  }
	}

	if (count < blen) {
		printf("Wrong length basis; too short!\n");
		exit(1);
	}

	sort_terms(tt, blen);

	return(tt);
}

int char_p(unsigned int degree)
{
	int i, test;
	int count;
	unsigned int a1, a2, a3, a4;
	
	count = 0;

	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
		test = 0;
		for (i = 0; i + 1 <= G.len; i++) {
			if ((G.ee[i]->e1 <= a1) &&
			(G.ee[i]->e2 <= a2) &&
			(G.ee[i]->e3 <= a3) &&
			(G.ee[i]->e4 <= a4) &&
			(G.ee[i]->e5 == 0)) test = 1;
		}
		if (!test) count++;
	      }
	    }
	  }
	}
	return(count);
}

int __extra(unsigned int degree, int *gap)
{
	int i, extra, e;
	unsigned int a1, a2, a3, a4;
	
	extra = 0;

	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
		e = -1;
		for (i = 0; i + 1 <= G.len; i++) {
			if ((G.ee[i]->e1 <= a1) &&
			(G.ee[i]->e2 <= a2) &&
			(G.ee[i]->e3 <= a3) &&
			(G.ee[i]->e4 <= a4) &&
			(G.ee[i]->e5 < gap[0])) {
				if ((e < 0) || (e > G.ee[i]->e5))
							e = G.ee[i]->e5;
			}
		}
		if (e > 0) extra = extra + e;
	      }
	    }
	  }
	}
	return(extra + 1);
}

/* Finds the char p generators of terms in degree degree.	*
 * This function assumes the function char_p has been		*
 * run previously and has returned a nonnegative integer blen.	*/
struct term **char_p_generators(unsigned int degree, int blen)
{
	int i, a1, a2, a3, a4, count, test;
	struct term tmp;
	struct term **tt;

	tt = (struct term **)malloc(blen*sizeof(struct term *));
	if (!tt) {
		perror("Malloc failed!");
		exit(1);
	}
	for (i = 0; i + 1 <= blen; i++) {
		tt[i] = NULL;
		make_term(&tt[i]);
	}
	
	count = 0;
	sc_one(tmp.c);
	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
		test = 0;
		for (i = 0; i + 1 <= G.len; i++) {
			if ((G.ee[i]->e1 <= a1) &&
			(G.ee[i]->e2 <= a2) &&
			(G.ee[i]->e3 <= a3) &&
			(G.ee[i]->e4 <= a4) &&
			(G.ee[i]->e5 == 0)) test = 1;
		}
		if (!test) {
			count++;
			if (count > blen) {
				printf("Wrong length basis!");
				exit(1);
			}
			/* tmp.c = 1 */
			tmp.n1 = a1;
			tmp.n2 = a2;
			tmp.n3 = a3;
			tmp.n4 = a4;
			copy_term(&tmp, tt[count - 1]);
			tt[count - 1]->next = NULL;
		}
	      }
	    }
	  }
	}

	if (count < blen) {
		printf("Wrong length basis; too short!\n");
		exit(1);
	}

	sort_terms(tt, blen);

	return(tt);
}

/* Takes the coefficients if aa completely reduces
 * otherwise returns NULL. */
mscalar *coefficients(
	struct polynomial **aa,
	int blen1, struct term **basis1,
	int blen2, struct term **basis2,
	int blen3, struct term **basis3)
{
	int row, fail;
	mscalar *column;
	struct term *aaterm;

	column = (mscalar *)
		malloc((blen1 + blen2 + blen3)*sizeof(mscalar));
	if (!column) {
		perror("Not again!\n");
		exit(1);
	}

	for (row = 0; row < blen1 + blen2 + blen3; row++) {
		sc_zero(column[row]);
	}

	fail = 0;
	row = 0;
	aaterm = aa[0]->leading;
	while ((aaterm) && (row < blen3)) {
		if (
		(aaterm->n1 == basis3[row]->n1) &&
		(aaterm->n2 == basis3[row]->n2) &&
		(aaterm->n3 == basis3[row]->n3) &&
		(aaterm->n4 == basis3[row]->n4)) {
			sc_copy(aaterm->c, column[row]);
			row++;
			aaterm = aaterm->next;
		} else {
			row++;
		}
	}
	if (aaterm) fail = 1;
	row = blen3;
	aaterm = aa[1]->leading;
	while ((aaterm) && (row < blen3 + blen2)) {
		if (
		(aaterm->n1 == basis2[row - blen3]->n1) &&
		(aaterm->n2 == basis2[row - blen3]->n2) &&
		(aaterm->n3 == basis2[row - blen3]->n3) &&
		(aaterm->n4 == basis2[row - blen3]->n4)) {
			sc_copy(aaterm->c, column[row]);
			row++;
			aaterm = aaterm->next;
		} else {
			row++;
		}
	}
	if (aaterm) fail = 1;
	row = blen3 + blen2;
	aaterm = aa[2]->leading;
	while ((aaterm) && (row < blen3 + blen2 + blen1)) {
		if (
		(aaterm->n1 == basis1[row - blen3 - blen2]->n1) &&
		(aaterm->n2 == basis1[row - blen3 - blen2]->n2) &&
		(aaterm->n3 == basis1[row - blen3 - blen2]->n3) &&
		(aaterm->n4 == basis1[row - blen3 - blen2]->n4)) {
			sc_copy(aaterm->c, column[row]);
			row++;
			aaterm = aaterm->next;
		} else {
			row++;
		}
	}
	if (aaterm) fail = 1;
	if (fail) {
		free(column);
		column = NULL;
	}
	free_tail(aa[0]->leading);
	free(aa[0]);
	free_tail(aa[1]->leading);
	free(aa[1]);
	free_tail(aa[2]->leading);
	free(aa[2]);
	free(aa);
	return(column);
}

static struct polynomial **split_up_term(struct term *tt)
{
	int degree;
	struct polynomial **aa;

	aa = (struct polynomial **)malloc(3*sizeof(struct polynomial *));
	if (!aa) {
		perror("Malloc failed!");
		exit(1);
	}
	aa[0] = NULL;
	aa[1] = NULL;
	aa[2] = NULL;
	make_pol(&aa[0]);
	make_pol(&aa[1]);
	make_pol(&aa[2]);

	aa[0]->degree = 3*d - d1 - d2 - d3 - d4;
	aa[1]->degree = 2*d - d1 - d2 - d3 - d4;
	aa[2]->degree = d - d1 - d2 - d3 - d4;

	if (sc_is_zero(tt->c)) return(aa);

	degree = tt->n1*d1 + tt->n2*d2 + tt->n3*d3 + tt->n4*d4;
	if (degree == aa[2]->degree) {
		make_term(&aa[2]->leading);
		copy_term(tt, aa[2]->leading);
	}
	if (degree == aa[1]->degree) {
		make_term(&aa[1]->leading);
		copy_term(tt, aa[1]->leading);
	}
	if (degree == aa[0]->degree) {
		make_term(&aa[0]->leading);
		copy_term(tt, aa[0]->leading);
	}

	return(all_the_way_split(aa));
}

mscalar **gens_to_basis(
	int blen1, struct term **basis1,
	int blen2, struct term **basis2,
	int blen3, struct term **basis3,
	int glen1, struct term **gens1,
	int glen2, struct term **gens2,
	int glen3, struct term **gens3,
	int *e)
{
	mscalar **matrix;
	mscalar *column;
	struct polynomial **aa;
	int i, j, k;
	mscalar c;

	matrix = (mscalar **)malloc((glen1 + glen2 + glen3)*sizeof(mscalar *));
	if (!matrix) {
		perror("Malloc failed!\n");
		exit(1);
	}

	sc_one(c);
	for (i = 1; i <= *e; i++) sc_imult_replace(p, c);
	i = 0;
	while (i < glen1 + glen2 + glen3) {
		if (i < glen3) {
			sc_copy(c, gens3[i]->c);
			aa = split_up_term(gens3[i]);
		} else if (i < glen3 + glen2) {
			sc_copy(c, gens2[i - glen3]->c);
			aa = split_up_term(gens2[i - glen3]);
		} else {
			sc_copy(c, gens1[i - glen3 - glen2]->c);
			aa = split_up_term(gens1[i - glen3 - glen2]);
		}
	
		column = coefficients(aa, blen1, basis1, blen2, basis2,
								blen3, basis3);
		if (!column) {
			(*e)++;
			sc_imult_replace(p, c);
			for (j = 0; j < i; j++) {
				for (k = 0; k < blen1 + blen2 + blen3; k++) {
					sc_imult_replace(p, matrix[j][k]);
				}
			}
		} else {
			matrix[i] = column;
			column = NULL;
			i++;
		}
	}

	return(matrix);
}

mscalar **prod_matrix(int n, int m, int l, mscalar **A, mscalar **B)
{
	int i, j, k;
	mscalar c;
	mscalar **C;

	C = (mscalar **)malloc(n*sizeof(mscalar *));
	if (!C) {
		perror("Malloc failed!\n");
		exit(1);
	}
	for (i = 0; i < n; i++) {
		C[i] = (mscalar *)malloc(l*sizeof(mscalar));
		if (!C[i]) {
			perror("Malloc failed!\n");
			exit(1);
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < l; j++) {
			sc_zero(C[i][j]);
			for (k = 0; k < m; k++) {
				sc_mult(A[i][k], B[k][j], c);
				sc_add_replace(c, C[i][j]);
			}
		}
	}
	return(C);
}			

void print_matrix(int rows, int columns, mscalar **matrix)
{
	int i, j;

	printf("[");
	for (i = 0; i < columns; i++) {
		for (j = 0; j < rows; j++) {
			printmscalar(matrix[j][i]);
			if (j + 1 < rows) printf(",");
		}
		if (i + 1 < columns) printf(";\\\n");
	}
	printf("]\n");
	return;
}

int clean_matrix(int rows, int columns, mscalar **matrix)
{
	int i, j, e, v;

	e = -1;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			if (!sc_is_zero(matrix[i][j])) {
				v = valuation(matrix[i][j]);
				if ((e < 0) || (e > v)) {
					e = v;
				}
			}
		}
	}
	if ((e == 0) || (e < 0)) return(0);
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			div_p(matrix[i][j], e);
		}
	}
	return(e);
}

void free_matrix(int rows, int columns, mscalar **matrix)
{
	int i;

	for (i = 0; i < rows; i++) {
		free(matrix[i]);
	}
	free(matrix);
	return;
}
