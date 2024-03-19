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
acb_theta_ql_ctx_dupl(acb_theta_ql_ctx_t ctx2, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong nbt = 1;
	slong j, k;

	/* Update exp_tau using squares (todo: conj for real entries?) */
	for (j = 0; j < g; j++)
	{
		for (k = j; k < g; k++)
		{
			acb_sqr(acb_mat_entry(ctx2->exp_tau, j, k), acb_mat_entry(ctx->exp_tau, j, k), prec);
			acb_sqr(acb_mat_entry(ctx2->exp_tau_inv, j, k), acb_mat_entry(ctx->exp_tau_inv, j, k), prec);
		}
	}

	/* Zero information */
	ctx2->t_is_zero = ctx->z_is_zero;
	ctx2->z_is_zero = ctx->z_is_zero;
	ctx2->z_is_real = ctx->z_is_real;

	/* Ones */
	for (j = 0; j < g; j++)
	{
		acb_one(&ctx2->exp_zs[j]);
		acb_one(&ctx2->exp_zs_inv[j]);
	}

	/* Update exp_t if present */
	if (!ctx->t_is_zero)
	{
		_acb_vec_set(ctx2->exp_zs + g, ctx->exp_zs + 2 * g, g);
		_acb_vec_set(ctx2->exp_zs_inv + g, ctx->exp_zs_inv + 2 * g, g);
		for (j = 0; j < g; j++)
		{
			acb_sqr(&ctx2->exp_zs[2 * g + j], &ctx2->exp_zs[g + j], prec);
			acb_conj(&ctx2->exp_zs_inv[2 * g + j], &ctx2->exp_zs[2 * g + j], prec);
		}
		nbt = 3;
	}

	/* Update exp_z using squares (todo: conj for real entries?) */
	if (!ctx->z_is_zero)
	{
		for (j = nbt * g; j < 2 * nbt * g; j++)
		{
			acb_sqr(&ctx2->exp_zs[j], &ctx->exp_zs[j], prec);
			acb_sqr(&ctx2->exp_zs_inv[j], &ctx->exp_zs_inv[j], prec);
		}
	}

	/* Higher genus stuff */
	if (g >= 2)
	{
		arb_t c;
		arb_init(c);

		arb_set_si(c, 2);
		arb_sqrt(c, c, prec);

		arb_mat_scalar_mul_2exp_si(ctx2->Y, ctx->Y, 1);
		arb_mat_scalar_mul_2exp_si(ctx2->Yinv, ctx->Yinv, -1);
		arb_mat_scalar_mul(ctx2->C, ctx->C, c, prec);
		arb_mat_scalar_div(ctx2->Cinv, ctx->Cinv, c, prec);

		arb_clear(c);
	}


}
