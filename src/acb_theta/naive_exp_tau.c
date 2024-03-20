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
acb_theta_naive_exp(acb_mat_t exp_tau, const acb_mat_t tau, slong prec)
{
	slong g = acb_mat_nrows(tau);
	slong j, k;

	/* Exponentials of tau */
	/* These are the square roots of the entries of exp_tau in naive_worker. */
	for (j = 0; j < g; j++)
	{
		for (k = j; k < g; k++)
		{
			acb_set_round(acb_mat_entry(exp_tau, j, k), acb_mat_entry(tau, j, k), prec);
			if (k == j)
			{
				acb_mul_2exp_si(acb_mat_entry(ctx->exp_tau, j, k), acb_mat_entry(ctx->exp_tau, j, k), -1);
			}
			acb_exp_pi_i(acb_mat_entry(ctx->exp_tau, j, k), acb_mat_entry(ctx->exp_tau, j, k));
		}
	}
}
