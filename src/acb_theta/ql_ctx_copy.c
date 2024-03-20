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
acb_theta_ql_ctx_copy(acb_theta_ql_ctx_t ctx2, const acb_theta_ql_ctx_t ctx);
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;

	acb_mat_set(acb_theta_ql_ctx_y(ctx2), acb_theta_ql_ctx_y(ctx));
	acb_mat_set(acb_theta_ql_ctx_exp_tau(ctx2), acb_theta_ql_ctx_exp_tau(ctx));

	ctx2->z_is_real = ctx->z_is_real;
	ctx2->z_is_zero = ctx->z_is_zero;
	ctx2->t_is_zero = ctx->t_is_zero;

	_acb_vec_set(acb_theta_ql_ctx_exp_zs(ctx2), acb_theta_ql_ctx_exp_zs(ctx), 6 * g);
	acb_set(acb_theta_ql_ctx_c(ctx2), acb_theta_ql_ctx_c(ctx));

	if (g >= 2)
	{
		arb_mat_set(acb_theta_ql_ctx_yinv(ctx2), acb_theta_ql_ctx_yinv(ctx));
		arb_mat_set(acb_theta_ql_ctx_cho(ctx2), acb_theta_ql_ctx_cho(ctx));
		arb_mat_set(acb_theta_ql_ctx_choinv(ctx2), acb_theta_ql_ctx_choinv(ctx));
		_arb_vec_set(acb_theta_ql_ctx_v(ctx2), acb_theta_ql_ctx_v(ctx));
		_acb_vec_set(acb_theta_ql_ctx_dists(ctx2), acb_theta_ql_ctx_dists(ctx), 2 * n);
	}
}
