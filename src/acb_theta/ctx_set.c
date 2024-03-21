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
acb_theta_ctx_set(acb_theta_ctx_t ctx, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
	slong g = acb_theta_ctx_g(ctx);
	acb_t x;
	arb_t u;
	arb_ptr y;
	slong j, k;
	int res;

	FLINT_ASSERT(nb <= acb_theta_ctx_nb(ctx));

	acb_init(x);
	arb_init(u);
	y = _arb_vec_init(g);

	/* Set Y, exp_tau, exp_z, exp_z_inv */

	for (j = 0; j < nb * g; j++)
	{
		acb_set_round(x, &zs[j], prec);
		acb_mul_2exp_si(x, x, 1);
		acb_exp_pi_i(&acb_theta_ctx_exp_zs(ctx)[j], x, prec);
		if (acb_is_real(&zs[j]))
		{
			acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[j], &acb_theta_ctx_exp_zs(ctx)[j]);
		}
		else
		{
			acb_inv(&acb_theta_ctx_exp_zs_inv(ctx)[j],  &acb_theta_ctx_exp_zs(ctx)[j], prec);
		}
	}

	/* Set exp_tau_inv, cs, vs */
	if (g == 1)
	{
		for (j = 0; j < nb; j++)
		{
			arb_sqr(u, acb_imagref(zs + j * g), prec);
			arb_neg(u, u);
			arb_const_pi(&acb_theta_ctx_exp_cs(ctx)[j], prec);
			arb_mul(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], u, prec);
			arb_exp(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], prec);
		}
	}
	else
	{
		}
		for (j = 0; j < nb; j++)
		{
			_acb_vec_get_imag(y, zs + j * g, g);
			_arb_vec_neg(acb_theta_ctx_vs(ctx) + j * g, y, g);
			arb_mat_vector_mul_col(acb_theta_ctx_vs(ctx) + j * g, acb_theta_ctx_yinv(ctx),
				acb_theta_ctx_vs(ctx) + j * g, prec);
			if (j == 0)
			{
				_arb_vec_set(acb_theta_ctx_v(ctx), acb_theta_ctx_vs(ctx), g);
			}
			else
			{
				_arb_vec_union(acb_theta_ctx_v(ctx), acb_theta_ctx_v(ctx),
					acb_theta_ctx_vs(ctx) + j * g, g, prec);
			}
			arb_dot(u, y, 1, acb_theta_ctx_vs(ctx) + j * g, 1, g, prec);
			arb_const_pi(&acb_theta_ctx_exp_cs(ctx)[j], prec);
			arb_mul(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], u, prec);
			arb_exp(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], prec);
		}
	}

	acb_clear(x);
	arb_clear(u);
	_arb_vec_clear(y, g);
}
