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
	slong j, k, a;
	slong lp = ACB_THETA_LOW_PREC;
	arb_ptr w;

	if (g >= 2)
	{
		acb_mat_get_imag(ctx->Y, tau);
		acb_siegel_yinv(ctx->Yinv, tau, prec);
		acb_siegel_cho(ctx->C, tau, prec);
		res = acb_mat_inv(ctx->Cinv, ctx->C, prec);
		if (!res)
		{
			acb_mat_indeterminate(ctx->Cinv);
		}
	}

	/* Exponentials of tau */
	/* These are the square roots of the entries of exp_tau in naive_worker. */
	for (j = 0; j < g; j++)
	{
		for (k = j; k < g; k++)
		{
			acb_set(acb_mat_entry(ctx->exp_tau, j, k), acb_mat_entry(tau, j, k));
			if (k == j)
			{
				acb_mul_2exp_si(acb_mat_entry(ctx->exp_tau, j, k), acb_mat_entry(ctx->exp_tau, j, k), -1);
			}
			acb_exp_pi_i(acb_mat_entry(ctx->exp_tau, j, k), acb_mat_entry(ctx->exp_tau, j, k));

			if (arb_is_zero(acb_imagref(acb_mat_entry(tau, j, k))))
			{
				acb_conj(acb_mat_entry(ctx->exp_tau_inv, j, k), acb_mat_entry(ctx->exp_tau, j, k));
			}
			else
			{
				acb_inv(acb_mat_entry(ctx->exp_tau_inv, j, k), acb_mat_entry(ctx->exp_tau, j, k));
			}
		}
	}

	/* t = 0 at initialization */
	ctx->z_is_real = _acb_vec_is_real(z, g);
	ctx->z_is_zero = _acb_vec_is_zero(z, g);
	ctx->t_is_zero = 1;

	/* Exponentials of 0 and z */
	/* These are the entries of exp_z as in naive_worker for the two vectors 0, z. */
	for (j = 0; j < g; j++)
	{
		acb_one(&exp_zs[j]);
		acb_one(&exp_zs_inv[j]);

		if (!ctx->z_is_zero)
		{
			acb_mul_2exp_si(&exp_zs[j + g], &z[j], 1);
			acb_exp_pi_i(&exp_zs[j + g], &exp_zs[j + g], prec);
			if (arb_is_zero(acb_imagref(&z[j])))
			{
				acb_conj(&exp_zs_inv[j + g], &exp_zs[j + g], prec);
			}
			else
			{
				acb_inv(&exp_zs_inv[j + g], &exp_zs[j + g], prec);
			}
		}
	}

	/* Center of ellipsoid and distances */
	if (g >= 2)
	{
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
