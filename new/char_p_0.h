/*
 *	char_p_0.h
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


int *find_gap(void );
int char_0(unsigned int degree, int *gap);
void print_terms(struct term **tt, int blen);
void free_list_terms(struct term **tt, int blen);
struct term **char_0_basis(unsigned int degree, int blen, int *gap);
int char_p(unsigned int degree);

int __extra(unsigned int degree, int *gap);

struct term **char_p_generators(unsigned int degree, int blen);

mscalar *coefficients(
	struct polynomial **aa,
	int blen1, struct term **basis1,
	int blen2, struct term **basis2,
	int blen3, struct term **basis3);
mscalar **gens_to_basis(
	int blen1, struct term **basis1,
	int blen2, struct term **basis2,
	int blen3, struct term **basis3,
	int glen1, struct term **gens1,
	int glen2, struct term **gens2,
	int glen3, struct term **gens3,
	int *e);
void print_matrix(int rows, int columns, mscalar **matrix);
int clean_matrix(int rows, int columns, mscalar **matrix);
void free_matrix(int rows, int columns, mscalar **matrix);
mscalar **prod_matrix(int n, int m, int l, mscalar **A, mscalar **B);


