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

	acb_mat_get_imag(acb_theta_ql_ctx_y(ctx), tau);

	if (g >= 2)
	{
		acb_siegel_yinv(acb_theta_ql_ctx_yinv(ctx), tau, prec);
		acb_siegel_cho(acb_theta_ql_ctx_cho(ctx), tau, prec);
		res = acb_mat_inv(acb_theta_ql_ctx_choinv(ctx), acb_theta_ql_ctx_cho(ctx), prec);
		if (!res)
		{
			acb_mat_indeterminate(acb_theta_ql_ctx_choinv(ctx));
		}
	}

	/* t = 0 at initialization */
	ctx->z_is_real = _acb_vec_is_real(z, g);
	ctx->z_is_zero = _acb_vec_is_zero(z, g);
	ctx->t_is_zero = 1;

	/* Exponentials */
	for (j = 0; j < g; j++)
	{
		acb_one(&acb_theta_ql_ctx_exp_zs(ctx)[j]);
	}
	if (!ctx->z_is_zero)
	{
		acb_theta_naive_exp_z(acb_theta_ql_ctx_exp_zs(ctx) + g, z, g, prec);
	}
	acb_theta_naive_exp_tau(acb_theta_ql_ctx_exp_tau(ctx), tau, prec);

	/* Center of ellipsoid, distances, cofactor */
	if (g == 1 && !ctx->z_is_zero)
	{
		arb_t x;
		arb_init(x);

		arb_sqr(x, acb_imagref(z));
		arb_div(x, x, acb_mat_entry(acb_theta_ql_ctx_y(ctx), 0, 0));
		arb_neg(x, x);
		arb_const_pi(acb_theta_ql_ctx_c(ctx), prec);
		arb_mul(acb_theta_ql_ctx_c(ctx), acb_theta_ql_ctx_c(ctx), x, prec);
		arb_exp(acb_theta_ql_ctx_c(ctx), acb_theta_ql_ctx_c(ctx), prec);

		arb_clear(x);
	}
	else if (g >= 2)
	{
		arb_ptr w, y;
		arb_t x;
		slong lp = ACB_THETA_LOW_PREC;
		slong a;

		w = _arb_vec_init(g);
		y = _arb_vec_init(g);
		arb_init(x);

		_acb_vec_get_imag(acb_theta_ql_ctx_v(ctx), z, g);
		_arb_vec_neg(acb_theta_ql_ctx_v(ctx), acb_theta_ql_ctx_v(ctx), g);
		arb_mat_vector_mul_col(acb_theta_ql_ctx_v(ctx), acb_theta_ql_ctx_yinv(ctx),
			acb_theta_ql_ctx_v(ctx), lp);

		for (a = 0; a < n; a++)
		{
			acb_theta_char_get_arb(w, a, g);
			_arb_vec_add(w, acb_theta_ql_ctx_v(ctx), w, g, lp);
			arb_mat_vector_mul_col(w, acb_theta_ql_ctx_cho(ctx), w, lp);
			acb_theta_dist_lat(&acb_theta_ql_ctx_dists(ctx)[a], w, acb_theta_ql_ctx_cho(ctx), lp);
		}
		if (!ctx->z_is_zero)
		{
			for (a = 0; a < n; a++)
			{
				acb_theta_char_get_arb(w, a, g);
				_arb_vec_add(w, acb_theta_ql_ctx_v(ctx), w, g, lp);
				arb_mat_vector_mul_col(w, acb_theta_ql_ctx_cho(ctx), w, lp);
				acb_theta_dist_lat(&acb_theta_ql_ctx_dists(ctx)[n + a], acb_theta_ql_ctx_cho(ctx), lp);
			}

			_acb_vec_get_imag(y, z, g);
			arb_mat_vector_mul_col(w, acb_theta_ql_ctx_yinv(ctx), y, prec);
			arb_dot(x, y, 1, w, 1, g, prec);
			arb_neg(x, x);
			arb_const_pi(acb_theta_ql_ctx_c(ctx), prec);
			arb_mul(acb_theta_ql_ctx_c(ctx), acb_theta_ql_ctx_c(ctx), x, prec);
			arb_exp(acb_theta_ql_ctx_c(ctx), acb_theta_ql_ctx_c(ctx), prec);
		}

		_arb_vec_clear(w, g);
		_arb_vec_clear(y, g);
		arb_clear(x);
	}
}
