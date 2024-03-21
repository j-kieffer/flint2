/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"


static int
acb_theta_ql_naive_from_ctx_gen(acb_ptr th, const acb_theta_ql_ctx_t ctx,
	int all, int only_roots, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	acb_ptr exp_zs;
	slong nb = 1;
	int success = 1;

	if (!ctx->t_is_zero && only_roots)
	{
		nb = 2;
	}
	else if (!ctx->t_is_zero)
	{
		nb = 3;
	}
	if (!ctx->z_is_zero && ctx->z_is_real)
	{
		nb *= 2;
	}

	_acb_vec_init(exp_zs, nb);

	/* Copy the right exponentials into exp_zs */
	if (ctx->t_is_zero || !only_roots)
	{
		_acb_vec_set(exp_zs, acb_theta_ql_ctx_exp_zs(ctx), nb * g);
	}
	else
	{
		_acb_vec_set(exp_zs, acb_theta_ql_ctx_exp_zs(ctx) + g, 2 * g);
		if (!ctx->z_is_zero && ctx->z_is_real)
		{
			_acb_vec_set(exp_zs + 2 * g, acb_theta_ql_ctx_exp_zs(ctx) + 4 * g, 2 * g);
		}
	}

	/* Call worker */
	acb_theta_ql_naive_worker_gen(th, exp_zs, nb, acb_theta_ql_ctx_exp_tau(ctx),
		acb_theta_ql_ctx_v(ctx), acb_theta_ql_ctx_dists(ctx), acb_theta_ql_ctx_cho(ctx),
		all, prec);

	if (!ctx->z_is_zero && !ctx->z_is_real)
	{
		/* Do it again with another ellipsoid. nb is 1, 2 or 3. */
		if (!ctx->t_is_zero && only_roots)
		{
			_acb_vec_set(exp_zs, acb_theta_ql_ctx_exp_zs(ctx) + 4 * g, 2 * g);
		}
		else
		{
			_acb_vec_set(exp_zs, acb_theta_ql_ctx_exp_zs(ctx) + 3 * g, nb * g);
		}
		acb_theta_ql_naive_worker_gen(th + nb * (all ? n : 1), exp_zs, nb,
			acb_theta_ql_ctx_exp_tau(ctx), acb_theta_ql_ctx_v(ctx) + g,
			acb_theta_ql_ctx_dists(ctx) + n, acb_theta_ql_ctx_cho(ctx), all, prec);
	}

	_acb_vec_clear(exp_zs, nb);
	return success;
}

static void
acb_theta_ql_naive_worker_g1(acb_ptr th, const acb_t exp_z, int w_is_unit,
	const acb_mat_t exp_tau, int all, slong prec)
{
	/* todo: need to multiply by q1/4 */
	if (all)
	{
		acb_modular_theta_sum(&th[3], &th[2], &th[0], &th[1], exp_z, w_is_unit,
			acb_mat_entry(exp_tau, 0, 0), 1, prec);
		acb_neg(&th[3], &th[3]);
	}
}

static void
acb_theta_ql_naive_from_ctx_g1(acb_ptr th, const acb_theta_ql_ctx_t ctx,
	int all, int only_roots, slong prec)
{
	arb_t x;
	slong lp = ACB_THETA_LOW_PREC;
	slong nbth = (all ? 4 : 2);
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	slong j;

	arb_init(x);

	/* taking the maximal possible prec is fine for g=1 */
	arb_const_pi(x, lp);
	arb_mul(x, x, arb_mat_entry(acb_theta_ql_ctx_y(ctx), 0, 0), lp);
	arb_mul_2exp_si(x, x, -2);
	newprec = prec + acb_theta_dist_addprec(x);

	/* we don't care if z is real or not, except for the w_is_unit argument. */
	if (!ctx->t_is_zero && only_roots)
	{
		acb_theta_ql_naive_worker_g1(th, acb_theta_ql_ctx_exp_zs(ctx) + 1, 1,
			acb_theta_ql_ctx_exp_tau(ctx), all, newprec);
		acb_theta_ql_naive_worker_g1(th + nbth, acb_theta_ql_ctx_exp_zs(ctx) + 2, 1,
			acb_theta_ql_ctx_exp_tau(ctx), all, newprec);
		if (!ctx->z_is_zero)
		{
			acb_theta_ql_naive_worker_g1(th + 2 * nbth, acb_theta_ql_ctx_exp_zs(ctx) + 4,
				ctx->z_is_real, acb_theta_ql_ctx_exp_tau(ctx), all, newprec);
			acb_theta_ql_naive_worker_g1(th + 3 * nbth, acb_theta_ql_ctx_exp_zs(ctx) + >5,
				ctx->z_is_real, acb_theta_ql_ctx_exp_tau(ctx), all, newprec);
		}
	}
	else
	{
		for (j = 0; j < nbt; j++)
		{
			acb_theta_ql_naive_worker_g1(th + j * nbth, acb_theta_ql_ctx_exp_zs(ctx) + j, 1,
				acb_theta_ql_ctx_exp_tau(ctx), all, newprec);
		}
		for (j = nbt; j < 2 * nbt; j++)
		{
			acb_theta_ql_naive_worker_g1(th + j * nbth, acb_theta_ql_ctx_exp_zs(ctx) + j,
				ctx->z_is_real, acb_theta_ql_ctx_exp_tau(ctx), all, newprec);
		}
	}

	acb_clear(x);
}

int
acb_theta_ql_a0_naive_from_ctx(acb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);

	if (g == 1)
	{
		acb_theta_ql_naive_from_ctx_g1(th, ctx, 0, 0, prec);
		return 1;
	}
	else
	{
		return acb_theta_ql_naive_from_ctx_gen(th, ctx, 0, 0, prec);
	}
}

int
acb_theta_ql_roots_from_ctx(acb_ptr rts, const acb_theta_ql_ctx_t ctx, slong guard)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	slong nbr = (ctx->t_is_zero ? 1 : 2) * (ctx->z_is_zero ? 1 : 2);
	slong nbt = (ctx->t_is_zero ? 1 : 2);
	slong j;
	int res;

	/* Compute roots */
	if (g == 1)
	{
		acb_theta_ql_a0_naive_from_ctx_g1(th, ctx, 0, 1, guard);
		res = 1;
	}
	else
	{
		res = acb_theta_ql_naive_from_ctx_gen(th, ctx, 0, 1, guard);
	}

	/* Rescale roots corresponding to a nonzero z */
	if (!ctx->z_is_zero)
	{
		_acb_vec_scalar_mul(th + n * nbt, th + n * nbt, n * nbt, acb_theta_ql_ctx_c(ctx), guard);
	}

	/* Check that roots don't contain zero. */
	for (j = 0; (j < n * nbr) && res; j++)
	{
		if (acb_contains_zero(&th[j]))
		{
			res = 0;
		}
	}

	return res;
}
