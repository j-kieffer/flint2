/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_ql_a0_naive_from_ctx_gen(cb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	arb_ptr new_v;
	acb_theta_eld_t E;
	arf_t R2, eps;
	acb_ptr aux, new_exp_zs, cs;
	slong new_prec;
	slong a, j;
	int success = 1;

	new_v = _arb_vec_init(g);
	arf_init(R2);
	arf_init(eps);
	new_exp_zs = _acb_vec_init(nbt * g);
	aux = _acb_vec_init(nbt * n); /* todo: could do 6 at once if z is real */
	cs = _acb_vec_init(nbt);

	/* Get theta_a0 at 0, t, 2t */
	for (a = 0; (a < n) && success; a++)
	{
		acb_theta_char_get_arb(new_v, a, g);
		_arb_vec_neg(new_v, new_v, g);
		_arb_vec_add(new_v, new_v, ctx->vs, g, prec);
		new_prec = prec + acb_theta_dist_addprec(&ctx->dists_a0[a]);
		acb_theta_naive_radius(R2, eps, &ctx->C, 0, new_prec);
		success = acb_theta_eld_set(E, &ctx->C, R2, new_v);

		if (!success) break;

		for (j = 0; j < nbt; j++)
		{
			acb_theta_naive_exp_translate(&cs[j], new_exp_zs + j * g, ctx->exp_zs + j * g,
				&ctx->exp_tau, a, prec);
		}

		acb_theta_naive_worker(aux, n, new_exp_zs,  nbt,
			&ctx->exp_tau, E, 0, new_prec, acb_theta_naive_0b_worker);
		for (j = 0; j < nbt; j++)
		{
			_acb_vec_scalar_mul(th + j * n * n + a * n, aux + j * n, &cs[j], n, prec);
		}
	}

	/* Get theta_a0 at z, z + t, z + 2t */
	if (!ctx->z_is_zero)
	{
		for (a = 0; (a < n) && success; a++)
		{
			acb_theta_char_get_arb(new_v, a, g);
			_arb_vec_neg(new_v, new_v, g);
			_arb_vec_add(new_v, new_v, ctx->vs + g, g, prec);
			new_prec = prec + acb_theta_dist_addprec(&ctx->dists_a0[n + a]);
			acb_theta_naive_radius(R2, eps, &ctx->C, 0, new_prec);
			success = acb_theta_eld_set(E, &ctx->C, R2, new_v);

			if (!success) break;

			acb_theta_naive_worker(aux, n, new_exp_zs,  nbt,
				&ctx->exp_tau, E, 0, new_prec, acb_theta_naive_worker_0b);
			for (j = 0; j < nbt; j++)
			{
				_acb_vec_scalar_mul(th + (nbt + j) * n * n + a * n, aux + j * n, &cs[j], n, prec);
			}
		}
	}

	_arb_vec_clear(new_v, g);
	arf_clear(R2);
	arf_clear(eps);
	_acb_vec_clear(new_exp_zs, nbt * g);
	_acb_vec_clear(aux, nbt * n);
	_acb_vec_clear(cs, nbt);
}

static void
acb_theta_ql_all_naive_from_ctx_g1(acb_ptr th, const acb_theta_ql_ctx ctx, slong prec)
{
	arb_t x, y;
	acb_t th1, th2, th3, th4;
	slong lp = ACB_THETA_LOW_PREC;
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	slong j;

	arb_init(x);
	arb_init(y);
	acb_init(th1);
	acb_init(th2);
	acb_init(th3);
	acb_init(th4);

	/* taking the maximal possible prec is fine for g=1 */
	arb_set(y, acb_imagref(acb_mat_entry(ctx->Y, 0, 0)));
	arb_const_pi(x, lp);
	arb_mul(x, x, y, lp);
	arb_mul_2exp_si(x, x, -2);
	newprec = prec + acb_theta_dist_addprec(x);

	/* todo: we don't need the theta constants to extract square roots... */
	for (j = 0; j < nbt; j++)
	{
		if (j == 0)
		{
			acb_modular_theta_const_sum(th2, th3, th4, acb_mat_entry(&ctx->exp_tau, 0, 0), newprec);
		}
		else
		{
			acb_modular_theta_sum(th1, th2, th3, th4, &ctx->exp_zs[j], 1,
				acb_mat_entry(&ctx->exp_tau, 0, 0), 1, newprec);
		}
		acb_set(&th[4 * j], th3);
		acb_set(&th[4 * j + 1], th4);
		acb_set(&th[4 * j + 2], th2);
		acb_neg(&th[4 * j + 3], th1);
	}

	if (!ctx->z_is_zero)
	{
		for (j = nbt; j < 2 * nbt; j++)
		{
			acb_modular_theta_sum(th1, th2, th3, th4, &ctx->exp_zs[j],
				ctx->z_is_real, acb_mat_entry(&ctx->exp_tau, 0, 0), 1, newprec);
			acb_set(&th[4 * j], th3);
			acb_set(&th[4 * j + 1], th4);
			acb_set(&th[4 * j + 2], th2);
			acb_neg(&th[4 * j + 3], th1);
		}
	}

	arb_clear(x);
	arb_clear(y);
	acb_clear(th1);
	acb_clear(th2);
	acb_clear(th3);
	acb_clear(th4);
}

void
acb_theta_ql_all_naive_from_ctx(acb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);

	if (g == 1)
	{
		acb_theta_ql_all_naive_from_ctx_g1(th, ctx, prec);
	}
	else
	{
		acb_theta_ql_a0_naive_from_ctx_gen(th, ctx, prec);
	}
}
