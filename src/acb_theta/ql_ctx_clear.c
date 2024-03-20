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
acb_theta_ql_ctx_clear(acb_theta_ql_ctx_t ctx)
{
	slong g = acb_mat_nrows(ctx->exp_tau);
	slong n = 1 << g;

	acb_mat_clear(ctx->exp_tau);
	_acb_vec_clear(ctx->exp_zs, 6 * g);
	arb_mat_clear(ctx->Y);

	if (g >= 2)
	{
		arb_mat_clear(ctx->Yinv);
		arb_mat_clear(ctx->C);
		arb_mat_clear(ctx->Cinv);
		_arb_vec_clear(vs, 2 * g);
		_arb_vec_clear(dists, 2 * n);
	}
}
