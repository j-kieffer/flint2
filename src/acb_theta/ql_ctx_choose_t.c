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
		_acb_vec_set(acb_theta_ql_ctx_exp_zs(ctx) + 3 * g, acb_theta_ql_ctx_exp_zs(ctx) + g, g);
	}
	ctx->t_is_zero = 0;

	/* Update exponentials of 0, t, 2t, z, z + t, z + 2t */
	acb_theta_naive_exp_z(acb_theta_ql_ctx_exp_zs(ctx) + g, t, g, prec);
	for (j = 0; j < g; j++)
	{
		acb_sqr(&acb_theta_ql_ctx_exp_zs(ctx)[2 * g + j],
			&acb_theta_ql_ctx_exp_zs(ctx)[g + j], prec);
		acb_mul(&acb_theta_ql_ctx_exp_zs(ctx)[4 * g + j],
			&acb_theta_ql_ctx_exp_zs(ctx)[g + j], &acb_theta_ql_ctx_exp_zs(ctx)[3 * g + j], prec);
		acb_mul(&acb_theta_ql_ctx_exp_zs(ctx)[5 * g + j],
			&acb_theta_ql_ctx_exp_zs(ctx)[2 * g + j], &acb_theta_ql_ctx_exp_zs(ctx)[3 * g + j], prec);
	}
}
