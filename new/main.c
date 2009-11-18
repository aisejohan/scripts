/*
 *	basis.c
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
#include "char_p_0.h"
#include "basis.h"

int main(void )
{
	int i, retry;

	printf("/* Comments begin... \n");
#ifdef KIJKEN
	printf("Debug is set! To unset do not define KIJKEN.\n");
#endif

	/* Setup the scalars. */
	setup_scalars();

	/* Seed the randomness. */
	set_seed(0);

	/* Get myf. */
	retry = 1;
	while (retry == 1) {
		myf = get_f();
		retry = setup(0);
	}

	/* Compute the frobenius matrix. */
	compute_frobenius_matrix();

	/* Free G and myf. */
	free_tail(myf.leading);
	for (i = 0; i + 1 <= G.len; i++) {
		free_tail(G.BC[i]->bc1.leading);
		free_tail(G.BC[i]->bc2.leading);
		free_tail(G.BC[i]->bc3.leading);
		free_tail(G.BC[i]->bc4.leading);
		free_tail(G.BC[i]->bc5.leading);
		free_tail(G.ff[i]->leading);
	}
	deallocate_GVMnewMMold();
	free_reserves();

	return(0);
}
