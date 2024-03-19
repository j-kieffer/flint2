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

	acb_mat_init(ctx->exp_tau, g, g);
	acb_mat_init(ctx->exp_tau_inv, g, g);
	ctx->exp_zs = _acb_vec_init(6 * g);
	ctx->exp_zs_inv = _acb_vec_init(6 * g);

	if (g >= 2)
	{
		arb_mat_init(ctx->Y, g, g);
		arb_mat_init(ctx->Yinv, g, g);
		arb_mat_init(ctx->C, g, g);
		arb_mat_init(ctx->Cinv, g, g);
		ctx->vs = _acb_vec_init(2 * g);
		ctx->dists_a0 = _acb_vec_init(2 * n);
	}
}
