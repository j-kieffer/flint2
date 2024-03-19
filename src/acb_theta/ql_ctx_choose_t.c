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
acb_theta_ql_ctx_choose_t(acb_theta_ql_ctx_t ctx, const acb_ptr t, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong j;

	if (ctx->t_is_zero) /* put exponentials of z as vector no 3*/
	{
		_acb_vec_set(ctx->exp_zs + 3 * g, ctx->exp_zs + g, g);
		_acb_vec_set(ctx->exp_zs_inv + 3 * g, ctx->exp_zs_inv + g, g);
	}
	ctx->t_is_zero = 0;

	/* Update exponentials of 0, t, 2t, z, z + t, z + 2t */
	for (j = 0; j < g; j++)
	{
		acb_mul_2exp_si(&ctx->exp_zs[g + j], &t[j], 1);
		acb_exp_pi_i(&ctx->exp_zs[g + j], &ctx->exp_zs[g + j], prec);
		acb_conj(&ctx->exp_zs_inv[g + j], &ctx->exp_zs[g + j]);

		acb_sqr(&ctx->exp_zs[2 * g + j], &ctx->exp_zs[g + j], prec);
		acb_conj(&ctx->exp_zs_inv[2 * g + j], &ctx->exp_zs[2 * g + j]);

		acb_mul(&ctx->exp_zs[4 * g + j], &ctx->exp_zs[g + j], &ctx->exp_zs[3 * g + j], prec);
		acb_mul(&ctx->exp_zs[5 * g + j], &ctx->exp_zs[2 * g + j], &ctx->exp_zs[3 * g + j], prec);
		acb_mul(&ctx->exp_zs_inv[4 * g + j], &ctx->exp_zs_inv[g + j], &ctx->exp_zs_inv[3 * g + j], prec);
		acb_mul(&ctx->exp_zs_inv[5 * g + j], &ctx->exp_zs_inv[2 * g + j], &ctx->exp_zs_inv[3 * g + j], prec);
	}
}
