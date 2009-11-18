/*
 *	pol.c
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

#include "data.h"
#include "scalar.h"
#include "pol.h"

#define NR_RESERVE	1024
static struct term *reserve = NULL;
static struct term *lists = NULL;

/* Makes a term.				*/
void make_term(struct term **mon)
{
#ifdef KIJKEN
	if (*mon) {
		printf("Term already exists!\n");
		exit(1);
	}
#endif
	if (!reserve) {
		int i;
		struct term *list;

		list = (struct term *)
				malloc((NR_RESERVE + 1)*sizeof(struct term));
		if (!list) {
			perror("Malloc failed.");
			exit(1);
		}
		reserve = &(list[0]);
		for (i = 0; i < NR_RESERVE - 1; i++) {
			list[i].next = &(list[i+1]);
		}
		list[NR_RESERVE - 1].next = NULL;
		list[NR_RESERVE].next = lists;
		lists = list;
	}
	*mon = reserve;
	reserve = reserve->next;
	(*mon)->next = NULL;
}

void make_pol(struct polynomial **f)
{
#ifdef KIJKEN
	if (*f) {
		printf("Polynomial already exists!\n");
		exit(1);
	}
#endif
	*f = (struct polynomial *) malloc (sizeof (struct polynomial));
	if (!*f) {
		perror("Malloc failed.");
		exit(1);
	}
	(*f)->leading=NULL;
}

/* Frees a term. */
void free_term(struct term *mon)
{
#ifdef KIJKEN
	if (!mon) {
		printf("free_term: term does not exist!\n");
		exit(1);
	}
#endif
	mon->next = reserve;
	reserve = mon;
	mon = NULL;
}

void free_reserves(void )
{
	int nr;
	struct term *list;

	nr = 0;
	while (lists) {
		nr++;
		list = lists;
		lists = list[NR_RESERVE].next;
		free(list);
	}
	printf("/* Freed %d lists. */\n", nr);
}

/* Copies data not pointer.						*/
/* sc_copy does allocation.						*/
void copy_term(struct term *mon1, struct term *mon2)
{
	sc_copy(mon1->c, mon2->c);
	mon2->n1 = mon1->n1;
	mon2->n2 = mon1->n2;
	mon2->n3 = mon1->n3;
	mon2->n4 = mon1->n4;
}

/* This frees memory starting with mon.					*
 * It does not set mon=NULL. This means you will/might			*
 * get a segfault if you use memory management wrong.			*/
void free_tail(struct term *mon)
{
	struct term *tmp;
	while (mon) {
		tmp = mon->next;
		free_term(mon);
		mon = tmp;
	}
	return;
}

/* Remove terms with zero coefficient from putative pol.	*
 * The result is a polynomial of the same degree as the input	*
 * which may have ``degree'' NOT zero! This is usefull in case	*
 * where we directly after using this add this to another pol	*
 * of the same degree as the input.				*/
void clean_pol(struct polynomial *pol)
{
	struct term *tmp;
	struct term **ptrterm;

	ptrterm = &(pol->leading);
	while (*ptrterm) {
		if (sc_is_zero((*ptrterm)->c)) {
			tmp = *ptrterm;
			*ptrterm = (*ptrterm)->next;
			free_term(tmp);
		} else {
			ptrterm = &((*ptrterm)->next);
		}
	}
	return;
}

/* Copies tail of pol.					*/
void copy_tail(struct term *mon, struct term **ptrterm)
{
	while (mon) {
		make_term(ptrterm);
		copy_term(mon, *ptrterm);
		ptrterm = &((*ptrterm)->next);
		mon = mon->next;
	}
}


/* Copies a polynomial.					*/
struct polynomial copy_pol(struct polynomial f)
{
	struct polynomial uit;
	uit.leading = NULL;
	
	uit.degree = f.degree;
	copy_tail(f.leading, &(uit.leading));
	return(uit);
}



/* Prints a polynomial. 				*/
void print_pol(struct polynomial f)
{
	struct term *fterm;

	fterm = f.leading;
	if (!fterm) {
		printf("0\n");
		return;
	}
	while (fterm) {
		if (!sc_is_zero(fterm->c)) {
			printmscalar(fterm->c);
			if (fterm->n1) printf("*x^%d", fterm->n1);
			if (fterm->n2) printf("*y^%d", fterm->n2);
			if (fterm->n3) printf("*z^%d", fterm->n3);
			if (fterm->n4) printf("*w^%d", fterm->n4);
			if (fterm->next) printf(" + ");
			else printf("\n");
		} else {
			printf("ZERO COEFFICIENT! FIXME!\n");
			exit(1);
		}
		fterm = fterm->next;
	}
	return;
}

/* This function adds polynomials of the same degree	*
 * into a new polynomial structure: 			*
 * The pols MUST have the same degree.			*
 * The pols may be zero.				*
 * The pols MUST NOT have terms which are zero.		*/
struct polynomial pol_add(struct polynomial f, struct polynomial g)
{
	mscalar c;
	struct polynomial uit;
	struct term *fterm, *gterm;
	struct term **ptrterm;
	uit.leading = NULL;

#ifdef KIJKEN
	if (f.degree != g.degree) {
		printf("Can't add these! Wrong degrees in pol_add.\n");
		exit(1);
	}
#endif

	uit.degree = f.degree;
	ptrterm = &(uit.leading); /* uit.leading will be set later.*/
	fterm=f.leading;
	gterm=g.leading;
	while (1) {
		if (!fterm) {
			*ptrterm = NULL;
			copy_tail(gterm, ptrterm);
			return(uit);
		}
		if (!gterm) {
			copy_tail(fterm, ptrterm);
			return(uit);
		}

		switch (kleiner(fterm, gterm)) {

			case GROTER:
			make_term(ptrterm);
			copy_term(fterm, *ptrterm);
			ptrterm = &((*ptrterm)->next);
			fterm = fterm->next;
			break;

			case KLEINER:
			make_term(ptrterm);
			copy_term(gterm, *ptrterm);
			ptrterm = &((*ptrterm)->next);
			gterm = gterm->next;
			break;

			default:
			/* GELIJK */
			sc_add(fterm->c, gterm->c, c);
			if (!sc_is_zero(c)) {
				make_term(ptrterm);
				sc_copy(c, (*ptrterm)->c);
				(*ptrterm)->n1 = fterm->n1;
				(*ptrterm)->n2 = fterm->n2;
				(*ptrterm)->n3 = fterm->n3;
				(*ptrterm)->n4 = fterm->n4;
				ptrterm = &((*ptrterm)->next);
			}
			fterm = fterm->next;
			gterm = gterm->next;
			break;
		}
	}
}

/* Almost the same as rep_pol_add:			*
 * 	replace f by (f+g)				*
 * 	empty out g					*
 * Here g may have terms that are zero but not f.	*/
void merge_add(struct polynomial *f, struct polynomial g)
{
	struct term *fterm, *gterm;
	struct term **ptrterm;
#ifdef KIJKEN
	if ((f->degree != g.degree) && (g.leading)) {
		printf("Can't add these! Wrong degrees in merge_add.\n");
		exit(1);
	}
#endif
	ptrterm = &(f->leading);
	fterm = f->leading;
	*ptrterm = NULL;
	gterm = g.leading;
	g.leading = NULL;
	while (1) {
		if (!fterm) {
			/* if (!gterm) return; */
			*ptrterm = gterm;
			while (*ptrterm)  {
				if (sc_is_zero((*ptrterm)->c)) {
					gterm = *ptrterm;
					*ptrterm = (*ptrterm)->next;
					free_term(gterm);
				} else {
					ptrterm = &((*ptrterm)->next);
				}
			}
			return;
		}

		if (!gterm) {
			*ptrterm = fterm;
			return;
		}

		switch (kleiner(fterm, gterm)) {

			case GROTER:
			*ptrterm = fterm;
			ptrterm = &(fterm->next);
			fterm = fterm->next;
			break;

			case KLEINER:
			/* Check for zero in g. */
			if (!sc_is_zero(gterm->c)) {
				*ptrterm = gterm;
				ptrterm = &(gterm->next);
				gterm = gterm->next;
			} else {
				*ptrterm = gterm->next;
				free_term(gterm);
				gterm = *ptrterm;
			}
			break;

			default:
			/* GELIJK */
			sc_add_replace(gterm->c, fterm->c);
			if (sc_is_zero(fterm->c)) {
				/* Here we use *ptrterm as temp	*
				 * storage. A little ugly.	*/
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			} else {
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			}
			*ptrterm = gterm->next;
			free_term(gterm);
			gterm = *ptrterm;
			break;
		}
		*ptrterm = NULL;
	}
}


/* Same as above but replace f by (f+g).		*
 * Here g may have terms that are zero but not f.	*/
void rep_pol_add(struct polynomial *f, struct polynomial g)
{
	struct term *fterm, *gterm;
	struct term **ptrterm;

#ifdef KIJKEN
	if (f->degree != g.degree) {
		printf("Can't add these! Wrong degrees in rep_pol_add.\n");
		printf("Degree f is %d and degree g is %d\n",
				f->degree, g.degree);
		print_pol(*f);
		printf("\n");
		print_pol(g);
		exit(1);
	}
#endif

	ptrterm = &(f->leading);
	fterm=f->leading;
	*ptrterm = NULL;
	gterm=g.leading;
	while (1) {
		if (!fterm) {
			/* if (!gterm) return; */
			/* copy_tail with check for zero */
			while (gterm)  {
				if (!sc_is_zero(gterm->c)) {
					make_term(ptrterm);
					copy_term(gterm, *ptrterm);
					ptrterm = &((*ptrterm)->next);
				}
				gterm = gterm->next;
			}
			return;
		}

		if (!gterm) {
			*ptrterm = fterm;
			return;
		}

		switch (kleiner(fterm, gterm)) {

			case GROTER:
			*ptrterm = fterm;
			ptrterm = &(fterm->next);
			fterm = fterm->next;
			*ptrterm = NULL;
			break;

			case KLEINER:
			/* Check for zero in g. */
			if (!sc_is_zero(gterm->c)) {
				make_term(ptrterm);
				copy_term(gterm, *ptrterm);
				ptrterm = &((*ptrterm)->next);
			}
			gterm = gterm->next;
			break;

			default:
			/* GELIJK */
			sc_add_replace(gterm->c, fterm->c);
			if (sc_is_zero(fterm->c)) {
				/* Here we use *ptrterm as temp	*
				 * storage. A little ugly.	*/
				*ptrterm = fterm->next;
				free_term(fterm);
				fterm = *ptrterm;
			} else {
				*ptrterm = fterm;
				ptrterm = &(fterm->next);
				fterm = fterm->next;
			}
			*ptrterm = NULL;
			gterm=gterm->next;
			break;
		}
	}
}

/* Replaces f by the product (possibly null).		*/
void times_int(int c, struct polynomial *f)
{
	struct term *tmp;
	struct term **ptrterm;

	ptrterm = &(f->leading);
	while (*ptrterm) {
		sc_imult_replace(c, (*ptrterm)->c);
		if (sc_is_zero((*ptrterm)->c)) {
			tmp = *ptrterm;
			*ptrterm = (*ptrterm)->next;
			free_term(tmp);
		} else {
			ptrterm = &((*ptrterm)->next);
		}
	}
	return;
}


/* Replaces f by the product (possibly null).		*/
void times_scalar(mscalar c, struct polynomial *f)
{
	struct term *tmp;
	struct term **ptrterm;

	ptrterm = &(f->leading);
	while (*ptrterm) {
		sc_mult_replace(c, (*ptrterm)->c);
		if (sc_is_zero((*ptrterm)->c)) {
			tmp = *ptrterm;
			*ptrterm = (*ptrterm)->next;
			free_term(tmp);
		} else {
			ptrterm = &((*ptrterm)->next);
		}
	}
	return;
}

/* Divides by p^k if possible. */
void div_p_pol(int k, struct polynomial *f)
{
	struct term *aaterm;

	aaterm = f->leading;
	while (aaterm) {
		if (valuation(aaterm->c) >= k) {
			div_p(aaterm->c, k);
			aaterm = aaterm->next;
		} else {
			printf("FIXME!\n");
			exit(1);
		}
	}
	return;
}

/* This function assumes the data structures f and g	*
 * have the same length, and f is initialized, but	*
 * g need not be. The result t*f is put into g, but	*
 * it may not be a polynomial since some terms may be	*
 * zero.						*/
void times_term(struct term t, struct polynomial f, struct polynomial *g)
{
	struct term *fterm, *gterm;

	g->degree = f.degree + d1*t.n1 + d2*t.n2 + d3*t.n3 + d4*t.n4;
	gterm = g->leading;
	fterm = f.leading;
	while (fterm) {
		sc_mult(t.c, fterm->c, gterm->c);
		gterm->n1 = t.n1 + fterm->n1;
		gterm->n2 = t.n2 + fterm->n2;
		gterm->n3 = t.n3 + fterm->n3;
		gterm->n4 = t.n4 + fterm->n4;
		fterm = fterm->next;
		gterm = gterm->next;
	}
	return;
}


/* Same as above but it creates and outputs the product. 	*
 * Same caveats as above.					*/
struct polynomial
make_times_term(struct term t, struct polynomial f)
{
	struct term *fterm;
	struct term **ptrterm;
	struct polynomial uit;
	uit.leading = NULL;

	uit.degree = f.degree + d1*t.n1 + d2*t.n2 + d3*t.n3 + d4*t.n4;
	ptrterm = &(uit.leading);
	fterm = f.leading;
	while (fterm) {
		make_term(ptrterm);
		sc_mult(t.c, fterm->c, (*ptrterm)->c);
		(*ptrterm)->n1 = t.n1 + fterm->n1;
		(*ptrterm)->n2 = t.n2 + fterm->n2;
		(*ptrterm)->n3 = t.n3 + fterm->n3;
		(*ptrterm)->n4 = t.n4 + fterm->n4;
		fterm = fterm->next;
		ptrterm = &((*ptrterm)->next);
	}
	return(uit);
}

/* Only clean up and do modulo modulus at the very end. */
static struct polynomial __pol_mult(struct term *tt, struct polynomial g)
{
	struct polynomial uit, tmppol;

	uit = make_times_term(*tt, g);
	tt = tt->next;

	tmppol = make_times_term(*tt, g);
	rep_pol_add(&uit, tmppol);
	tt = tt->next;

	while (tt) {
		times_term(*tt, g, &tmppol);
		rep_pol_add(&uit, tmppol);
		tt = tt->next;
	}

	free_tail(tmppol.leading);

	return(uit);
}

struct polynomial pol_mult(struct polynomial f, struct polynomial g)
{
	struct term *tf, *tg;
	struct polynomial uit;

	if ((!f.leading) || (!g.leading)) {
		uit.leading = NULL;
		uit.degree = f.degree + g.degree;
		return uit;
	}

	if (!f.leading->next) {
		uit = make_times_term(*f.leading, g);
		clean_pol(&uit);
		return uit;
	}

	if (!g.leading->next) {
		uit = make_times_term(*g.leading, f);
		clean_pol(&uit);
		return uit;
	}
	
	tf = f.leading->next;
	tg = g.leading->next;

	while ((tf) && (tg)) {
		tf = tf->next;
		tg = tg->next;
	}

	if (tf) {
		return __pol_mult(g.leading, f);
	} else {
		return __pol_mult(f.leading, g);
	}
}
