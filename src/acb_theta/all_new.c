/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_all_new(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
	int sqr, slong prec);
{
    slong g = acb_mat_nrows(tau);
	slong n = 1 << g;
    fmpz_mat_t mat;
	acb_mat_t w, c, cinv;
	acb_ptr new_zs;
	slong j;

	FLINT_ASSERT(nb >= 0);
	if (nb == 0) return;

	fmpz_mat_init(mat, g, g);
	acb_mat_init(w, g, g);
	acb_mat_init(c, g, g);
	acb_mat_init(cinv, g, g);
	new_zs = _acb_vec_init(nb * g);

    acb_siegel_reduce(mat, tau, prec);
	acb_siegel_transform_cocycle_inv(w, c, cinv, mat, tau, prec);
	if (acb_siegel_is_reduced(w, -10, prec))
	{
		acb_mat_transpose(cinv, cinv);
		for (j = 0; j < nb; j++)
		{
			acb_mat_vector_mul_col(new_zs + j * g, cinv, zs + j * g, prec);
		}
		acb_theta_all_notransform(th, new_zs, nb, w, sqr, prec);

		/* transformation formula... */
	}
	else
	{
		_acb_vec_indeterminate(th, n * nb);
	}
}
