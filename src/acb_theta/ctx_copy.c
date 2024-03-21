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
acb_theta_ctx_copy(acb_theta_ctx_t ctx2, const acb_theta_ctx_t ctx);
{
	slong g = acb_theta_ctx_g(ctx);
	slong nb = acb_theta_ctx_nb(ctx);
	slong n = 1 << g;

	FLINT_ASSERT(g == acb_theta_ctx_g(ctx2));
	FLINT_ASSERT(nb == acb_theta_ctx_nb(ctx2));

	acb_mat_set(acb_theta_ctx_y(ctx2), acb_theta_ctx_y(ctx));
	acb_mat_set(acb_theta_ctx_yinv(ctx2), acb_theta_ctx_yinv(ctx));
	acb_mat_set(acb_theta_ctx_exp_tau_div_4(ctx2), acb_theta_ctx_exp_tau_div_4(ctx));
	acb_mat_set(acb_theta_ctx_exp_tau_div_2(ctx2), acb_theta_ctx_exp_tau_div_2(ctx));
	acb_mat_set(acb_theta_ctx_exp_tau(ctx2), acb_theta_ctx_exp_tau(ctx));
	_acb_vec_set(acb_theta_ctx_exp_zs(ctx2), acb_theta_ctx_exp_zs(ctx), nb * g);
	_acb_vec_set(acb_theta_ctx_exp_zs_inv(ctx2), acb_theta_ctx_exp_zs_inv(ctx), nb * g);
	_arb_vec_set(acb_theta_ctx_cs(ctx2), acb_theta_ctx_cs(ctx), nb);

	ctx2->z_is_real = ctx->z_is_real;
	ctx2->z_is_zero = ctx->z_is_zero;
	ctx2->t_is_zero = ctx->t_is_zero;

	if (g >= 2)
	{
		arb_mat_set(acb_theta_ctx_cho(ctx2), acb_theta_ctx_cho(ctx));
		arb_mat_set(acb_theta_ctx_choinv(ctx2), acb_theta_ctx_choinv(ctx));
		acb_mat_set(acb_theta_ctx_exp_tau_inv(ctx2), acb_theta_ctx_exp_tau_inv(ctx));
		_arb_vec_set(acb_theta_ctx_vs(ctx2), acb_theta_ctx_vs(ctx), nb * g);
		_arb_vec_set(acb_theta_ctx_d0(ctx2), acb_theta_ctx_d0(ctx), n);
		_arb_vec_set(acb_theta_ctx_d(ctx2), acb_theta_ctx_d(ctx), n);
	}
}
