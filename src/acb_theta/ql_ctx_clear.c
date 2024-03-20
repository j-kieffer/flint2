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
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;

	acb_mat_clear(acb_theta_ql_ctx_exp_tau(ctx));
	_acb_vec_clear(acb_theta_ql_ctx_exp_zs(ctx), 6 * g);
	arb_mat_clear(acb_theta_ql_ctx_y(ctx));
	acb_clear(acb_theta_ql_ctx_c(ctx));

	if (g >= 2)
	{
		arb_mat_clear(acb_theta_ql_ctx_yinv(ctx));
		arb_mat_clear(acb_theta_ql_ctx_cho(ctx));
		arb_mat_clear(acb_theta_ql_ctx_choinv(ctx));
		_arb_vec_clear(acb_theta_ql_ctx_v(ctx), 2 * g);
		_arb_vec_clear(acb_theta_ql_ctx_dists(ctx), 2 * n);
	}
}
