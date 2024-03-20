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
acb_theta_ql_ctx_set(acb_theta_ql_ctx_t ctx, acb_srcptr z, const acb_mat_t tau, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	int res;
	slong j;

	acb_mat_get_imag(ctx->Y, tau);

	if (g >= 2)
	{
		acb_siegel_yinv(ctx->Yinv, tau, prec);
		acb_siegel_cho(ctx->C, tau, prec);
		res = acb_mat_inv(ctx->Cinv, ctx->C, prec);
		if (!res)
		{
			acb_mat_indeterminate(ctx->Cinv);
		}
	}

	/* t = 0 at initialization */
	ctx->z_is_real = _acb_vec_is_real(z, g);
	ctx->z_is_zero = _acb_vec_is_zero(z, g);
	ctx->t_is_zero = 1;

	/* Exponentials */
	for (j = 0; j < g; j++)
	{
		acb_one(&exp_zs[j]);
	}
	if (!ctx->z_is_zero)
	{
		acb_theta_naive_exp_z(ctx->exp_zs + g, z, g, prec);
	}
	acb_theta_naive_exp_tau(ctx->exp_tau, tau, prec);

	/* Center of ellipsoid and distances */
	if (g >= 2)
	{
		arb_ptr w;
		slong lp = ACB_THETA_LOW_PREC;
		slong a;

		w = _arb_vec_init(g);

		_acb_vec_zero(ctx->vs, g);
		_acb_vec_get_imag(ctx->vs + g, z, g);
		_arb_vec_neg(ctx->vs + g, ctx->vs + g, g);
		arb_mat_vector_mul_col(ctx->vs + g, ctx->Yinv, ctx->vs + g, lp);

		for (a = 0; a < n; a++)
		{
			acb_theta_char_get_arb(w, a, g);
			_arb_vec_add(w, ctx->vs, w, g, lp);
			arb_mat_vector_mul_col(w, ctx->C, w, lp);
			acb_theta_dist_lat(ctx->dists_a0 + a, w, ctx->C, lp);
		}
		if (!ctx->z_is_zero)
		{
			for (a = 0; a < n; a++)
			{
				acb_theta_char_get_arb(w, a, g);
				_arb_vec_add(w, ctx->vs + g, w, g, lp);
				arb_mat_vector_mul_col(w, ctx->C, w, lp);
				acb_theta_dist_lat(ctx->dists_a0 + n + a, w, ctx->C, lp);
			}
		}

		_acb_vec_clear(w, g);
	}
}
