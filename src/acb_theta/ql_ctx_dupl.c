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
acb_theta_ql_ctx_dupl(acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong nbt = 1;
	slong j, k;

	arb_mat_scalar_mul_2exp_si(acb_theta_ql_ctx_y(ctx), acb_theta_ql_ctx_y(ctx), 1);
	arb_sqr(acb_theta_ql_ctx_c(ctx), acb_theta_ql_ctx_c(ctx), prec);

	/* Update exponentials using squarings */
	for (j = 0; j < g; j++)
	{
		for (k = j; k < g; k++)
		{
			acb_sqr(acb_mat_entry(acb_theta_ql_ctx_exp_tau(ctx), j, k),
				acb_mat_entry(acb_theta_ql_ctx_exp_tau(ctx), j, k), prec);
		}
	}

	if (!ctx->t_is_zero)
	{
		_acb_vec_set(acb_theta_ql_ctx_exp_zs(ctx) + g, acb_theta_ql_ctx_exp_zs(ctx) + 2 * g, g);
		for (j = 0; j < g; j++)
		{
			acb_sqr(&acb_theta_ql_ctx_exp_zs(ctx)[2 * g + j],
				&acb_theta_ql_ctx_exp_zs(ctx)[g + j], prec);
		}
		nbt = 3;
	}

	if (!ctx->z_is_zero)
	{
		for (j = nbt * g; j < 2 * nbt * g; j++)
		{
			acb_sqr(&acb_theta_ql_ctx_exp_zs(ctx)[j], &acb_theta_ql_ctx_exp_zs(ctx)[j], prec);
		}
	}

	/* Higher genus things */
	if (g >= 2)
	{
		arb_t sqrt2;
		arb_init(sqrt2);

		arb_set_si(sqrt2, 2);
		arb_sqrt(sqrt2, sqrt2, prec);

		arb_mat_scalar_mul_2exp_si(acb_theta_ql_ctx_yinv(ctx), acb_theta_ql_ctx_yinv(ctx), -1);
		arb_mat_scalar_mul(acb_theta_ql_ctx_cho(ctx), acb_theta_ql_ctx_cho(ctx), sqrt2, prec);
		arb_mat_scalar_div(acb_theta_ql_ctx_choinv(ctx), acb_theta_ql_ctx_choinv(ctx), sqrt2, prec);

		arb_clear(sqrt2);
	}
}
