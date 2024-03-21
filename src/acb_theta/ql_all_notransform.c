/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* compute c0, c1, c2 */

int acb_theta_ql_all_midpoint(acb_ptr th, acb_srcptr zs, slong nb,
	const acb_mat_t tau, int sqr, const arb_t c, const arb_t rho, slong prec)
{
	slong g = acb_mat_nrows(tau);
	acb_ptr zs_mid;
	acb_mat_t tau_mid;
	slong j, k;
	int res;

	FLINT_ASSERT(nb >= 0);
	if (nb == 0) return 1;

	for (j = 0; j < nb * g; j++)
	{
		acb_get_mid(&zs_mid[j], &zs[j]);
	}
	for (j = 0; j < g; j++)
	{
		for (k = j; k < g; k++)
		{
			acb_get_mid(acb_mat_entry(tau_mid, j, k), acb_mat_entry(tau, j, k));
			acb_set(acb_mat_entry(tau_mid, k, j), acb_mat_entry(tau_mid, j, k));
		}
	}

	/* todo: compute error bound */
	if (sqr)
	{
		/* todo: could adjust prec in terms of the error bound we just computed. */
		res = acb_theta_ql_all_sqr_exact(th, zs_mid, nb, tau_mid, prec);
	}
	else
	{
		res = acb_theta_ql_all_exact(th, zs_mid, nb, tau_mid, prec);
	}

	for (j = 0; j < n * nb; j++)
	{
		acb_add_error_arb();
	}
}

