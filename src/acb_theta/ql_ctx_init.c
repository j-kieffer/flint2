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
acb_theta_ql_ctx_init(acb_theta_ql_ctx_t ctx, slong g)
{
	slong n = 1 << g;

	acb_mat_init(acb_theta_ql_ctx_exp_tau(ctx), g, g);
	acb_theta_ql_ctx_exp_zs(ctx) = _acb_vec_init(6 * g);
	arb_mat_init(acb_theta_ql_ctx_y(ctx), g, g);
	acb_init(acb_theta_ql_ctx_c(ctx));

	ctx->z_is_zero = 1;
	ctx->z_is_real = 1;
	ctx->t_is_zero = 1;

	if (g >= 2)
	{
		arb_mat_init(acb_theta_ql_ctx_yinv(ctx), g, g);
		arb_mat_init(acb_theta_ql_ctx_cho(ctx), g, g);
		arb_mat_init(acb_theta_ql_ctx_choinv(ctx), g, g);
		acb_theta_ql_ctx_v(ctx) = _acb_vec_init(g);
		acb_theta_ql_ctx_dists(ctx) = _acb_vec_init(2 * n);
	}
}
