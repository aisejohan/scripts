/*
 *	helper.c
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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "data.h"
#include "scalar.h"
#include "pol.h"
#include "helper.h"

unsigned int ivaluation(int a)
{
	unsigned int uit = 0;
	while (!(a % p)) {
		a = a / p;
		uit++;
	}
	return(uit);
}

void set_seed(unsigned int zaadje)
{
	int fd, uit;
	unsigned int willekeurig;
	if (zaadje) {
		srand(zaadje);
		return;
	}
	fd = open("/dev/urandom", O_RDONLY);
	if (fd < 0) {
		printf("Unable to open /dev/urandom. So seed=666.\n");
		willekeurig = 666;
	} else {
		uit = read(fd, &willekeurig, sizeof(willekeurig));
		if (uit <= 0) {
			printf("Failure reading /dev/urandom. So Seed=666.\n");
			willekeurig = 666;
		}
	}
	srand(willekeurig);
	uit = close(fd);
	return;
}

unsigned int count_sum(int degree)
{
	unsigned int count, a1, a2, a3;

	if (degree < 0) return(0);

	count = 0;
	for (a1 = 0; (d1*a1 <= degree); a1++) {
		for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
			for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
				if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4
									== 0) {
					count++;
				}
			}
		}
	}
	return(count);
}

int hilbert(int degree)
{
	int goodcount;

	goodcount = count_sum(degree);
	goodcount -= count_sum(degree - d + d1);
	goodcount -= count_sum(degree - d + d2);
	goodcount -= count_sum(degree - d + d3);
	goodcount -= count_sum(degree - d + d4);
	goodcount += count_sum(degree - 2*d + (d1 + d2));
	goodcount += count_sum(degree - 2*d + (d1 + d3));
	goodcount += count_sum(degree - 2*d + (d1 + d4));
	goodcount += count_sum(degree - 2*d + (d2 + d3));
	goodcount += count_sum(degree - 2*d + (d2 + d4));
	goodcount += count_sum(degree - 2*d + (d3 + d4));
	goodcount -= count_sum(degree - 3*d + (d1 + d2 + d3));
	goodcount -= count_sum(degree - 3*d + (d1 + d2 + d4));
	goodcount -= count_sum(degree - 3*d + (d1 + d3 + d4));
	goodcount -= count_sum(degree - 3*d + (d2 + d3 + d4));
	goodcount += count_sum(degree - 4*d + (d1 + d2 + d3 + d4));

	return(goodcount);
}

#if 0
static void print_sum(unsigned int degree)
{
	unsigned int a1, a2, a3;
	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		printf("[%d, %d, %d, %d]\n", a1, a2, a3,
				(degree - (a1*d1 + a2*d2 + a3*d3))/d4);
	      }
	    }
	  }
	}
	return;
}
#endif

/* Makes a polynomial of degree degree with all terms		*
 * having coefficient 1.					*/
static struct polynomial make_full(unsigned int degree)
{
	unsigned int a1, a2, a3, a4;
	struct polynomial uit;
	struct term *uitterm;
	struct term **ptrterm;
	uitterm = NULL;
	uit.degree = degree;
	uit.leading = NULL;

	for (a1 = 0; (d1*a1 <= degree); a1++) {
	  for (a2 = 0; (d1*a1 + d2*a2 <= degree); a2++) {
	    for (a3 = 0; (d1*a1 + d2*a2 + d3*a3 <= degree); a3++) {
	      if ((degree - (a1*d1 + a2*d2 + a3*d3)) % d4 == 0) {
		a4 = (degree - (a1*d1 + a2*d2 + a3*d3))/d4;
		make_term(&uitterm);
		uitterm->n1 = a1;
		uitterm->n2 = a2;
		uitterm->n3 = a3;
		uitterm->n4 = a4;
		ito_sc(1, uitterm->c);
		ptrterm = &(uit.leading);
		while ((*ptrterm) && (kleiner(uitterm, *ptrterm) == KLEINER)) {
			ptrterm = &((*ptrterm)->next);
		}
		uitterm->next = *ptrterm;
		*ptrterm = uitterm;
		uitterm=NULL;
	      }
	    }
	  }
	}
	return(uit);
}

/* Makes a random polynomial of degree degree.			*
 * The result may terms which have coefficient zero.		*/
struct polynomial make_random(unsigned int degree)
{
	int c;
	struct polynomial uit;
	struct term *uitterm;
	uitterm = NULL;

	uit = make_full(degree);

	uitterm = uit.leading;

	while (uitterm) {
		c = rand() % p;
		/* Change for compatibility with previous version.	*
		 * This does not make any difference to the pol mod p.	*/
		if (c < -(p - 1)/2) c += p; /* OK, this never happens. */
		if (c > (p - 1)/2) c -= p; /* This does happen. */
		ito_sc(c, uitterm->c);
		uitterm = uitterm->next;
	}
	return(uit);
}

static void list_print(struct polynomial f)
{
	int c;
	unsigned int a1, a2, a3, a4;
	struct term *uitterm;

	uitterm = f.leading;
	while (uitterm) {
		a1 = uitterm->n1;
		a2 = uitterm->n2;
		a3 = uitterm->n3;
		a4 = uitterm->n4;
		c = 0;
		printf("Coefficient of   ");
		if (a1) {
			printf("x^%d", a1);
			c++;
		}
		if ((a1) && (a2+a3+a4)) {
			printf(" * ");
			c++;
		}
		if (a2) {
			printf("y^%d", a2);
			c++;
		}
		if ((a2) && (a3 + a4)) {
			printf(" * ");
			c++;
		}
		if (a3) {
			printf("z^%d", a3);
			c++;
		}
		if ((a3) && (a4)) {
			printf(" * ");
			c++;
		}
		if (a4) {
			printf("w^%d", a4);
			c++;
		}
		while (8 - c) {
			printf("   ");
			c++;
		}
		printf("= ");
#ifdef OUTPUT_LIST
 		printf("\n");
#else
 #ifdef INPUT_F
		scanf("%d", &c);
		ito_sc(c, uitterm->c);
 #else
		printmscalar(uitterm->c);
 		printf("\n");
 #endif
#endif
		uitterm = uitterm->next;
	}
}

struct polynomial get_f(void )
{
	struct polynomial uit;
#ifdef OUTPUT_LIST
	uit = make_full(d);
	list_print(uit);
	exit(0);
#else
 #ifdef INPUT_F
	uit = make_full(d);
	printf("\n");
	printf("Please input coefficients below.\n");
	list_print(uit);
	clean_pol(&uit);
	return(uit);
 #else
	uit = make_random(d);
	list_print(uit);
	clean_pol(&uit);
	return(uit);
 #endif
#endif
}


/* Counts the number of terms and verifies f is a polynomial.	*
 * Allows f to be the zero polynomial, warns if degree not 0.	*/
unsigned int number_terms(struct polynomial f)
{
	int count = 0;
	struct term *fterm;

	if (!f.leading) {
		if (!f.degree) {
			printf("Zero but not degree 0.");
		}
		return(0);
	}
	fterm = f.leading;
	while (fterm) {
		count++;
		fterm = fterm->next;
	}
	return(count);
}

/* Frobenius operation by raising varaibles tot the	*
 * pth power.						*/
struct polynomial frobenius(struct polynomial f)
{
	struct polynomial uit;
	struct term *fterm;
	struct term **ptrterm;
	uit.leading = NULL;
	ptrterm = NULL;

	uit.degree = p * f.degree;
	ptrterm = &(uit.leading);
	fterm = f.leading;
	while (fterm) {
		make_term(ptrterm);
		sc_copy(fterm->c, (*ptrterm)->c);
		(*ptrterm)->n1 = p*fterm->n1;
		(*ptrterm)->n2 = p*fterm->n2;
		(*ptrterm)->n3 = p*fterm->n3;
		(*ptrterm)->n4 = p*fterm->n4;
		ptrterm = &((*ptrterm)->next);
		fterm = fterm->next;
	}
	return(uit);
}

/* Ouputs the derivative. 					*
 * The result is nonsense if the degree of f is too low.	*/
struct polynomial deriv(struct polynomial f, unsigned int i)
{
	mscalar c;
	struct polynomial uit;
	struct term *fterm;
	struct term **ptrterm;
	uit.leading = NULL;
	ptrterm = NULL;

	fterm = f.leading;
	switch (i) {
		case 1: 
		uit.degree = (f.degree > d1) ? (f.degree - d1) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n1, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1 - 1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		return(uit);

		case 2: 
		uit.degree = (f.degree > d2) ? (f.degree - d2) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n2, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2 - 1;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		return(uit);

		case 3:
		uit.degree = (f.degree > d3) ? (f.degree - d3) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n3, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3 - 1;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		return(uit);
		
		case 4:
		uit.degree = (f.degree > d4) ? (f.degree - d4) : 0;
		ptrterm = &(uit.leading);
		while (fterm) {
			sc_imult(fterm->n4, fterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4 - 1;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
		}
		return(uit);
		
		default:
		printf("Wrong again honey!");
		exit(1);
	}
	exit(1);
}

/* Same as above. Replaces f by the derivative. 		*
 * The result is nonsense if the degree of f is too low.	*/
void rep_deriv(struct polynomial *f, unsigned int i)
{
	struct term *fterm;
	struct term **ptrterm;

	fterm = f->leading;
	ptrterm = &(f->leading);
	switch (i) {
		case 1: 
		f->degree = (f->degree > d1) ? (f->degree - d1) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n1, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n1 = fterm->n1 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		case 2: 
		f->degree = (f->degree > d2) ? (f->degree - d2) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n2, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n2 = fterm->n2 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		case 3: 
		f->degree = (f->degree > d3) ? (f->degree - d3) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n3, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n3 = fterm->n3 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		case 4:
		f->degree = (f->degree > d4) ? (f->degree - d4) : 0;
		while (fterm) {
			sc_imult_replace(fterm->n4, fterm->c);
			if (!sc_is_zero(fterm->c)) {
				fterm->n4 = fterm->n4 - 1;
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			} else {
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			}
		}
		return;

		default:
		printf("Wrong again honey!");
		exit(1);
	}
	exit(1);
}
