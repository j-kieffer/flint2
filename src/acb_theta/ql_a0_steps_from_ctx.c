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
acb_theta_ql_step_1(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;

    acb_theta_agm_mul_tight(res, th0, th, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(res, res, n, g);
    acb_theta_agm_sqrt(res, res, rts, n, prec);
}

static void
acb_theta_ql_step_2(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;

    aux = _acb_vec_init(3 * n);

    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    _acb_vec_set(res, aux, 3 * n);
    _acb_vec_clear(aux, 3 * n);
}

static void
acb_theta_ql_step_3(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr aux;
    ulong a;

    aux = _acb_vec_init(3 * n);

    acb_theta_agm_mul_tight(aux + n, th0, th + n, d0, d, g, prec);
    acb_theta_agm_mul_tight(aux + 2 * n, th0, th + 2 * n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux + n, aux + n, 2 * n, g);
    acb_theta_agm_sqrt(aux + n, aux + n, rts, 2 * n, prec);

    acb_theta_agm_mul_tight(aux, th0 + n, th + n, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(aux, aux, n, g);
    for (a = 0; a < n; a++)
    {
        acb_div(&aux[a], &aux[a], &aux[2 * n + a], prec);
    }

    _acb_vec_set(res, aux, 3 * n);
    _acb_vec_clear(aux, 3 * n);
}

static void
acb_theta_ql_a0_step(acb_ptr th, acb_srcptr rts, arb_srcptr d0, arb_srcptr d,
	const acb_theta_ql_ctx_t ctx, int last_step, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	slong nbz = (ctx->z_is_zero ? 1 : 2) * (ctx->t_is_zero ? 1 : 3);
	acb_ptr next;

	next = _acb_vec_init(nbz * n);

	if (ctx->t_is_zero)
	{
		acb_theta_ql_step_1(next, th, th, rts, d0, d0, g, prec);
		if (!ctx->z_is_zero)
		{
			acb_theta_ql_step_1(next + n, th, th + n, rts + n, d0, d, g, prec);
		}
	}
	else
	{
        acb_theta_ql_step_3(next, th, th, rts, d0, d0, g, prec);
        if (!ctx->z_is_zero)
		{
			if (last_step)
			{
				acb_theta_ql_step_3(next + 3 * n, th, th + 3 * n,
					rts + 2 * n, d0, d, g, prec);
			}
			else
			{
				acb_theta_ql_step_2(next + 3 * n, th, th + 3 * n,
					rts + 2 * n, d0, d, g, prec);
			}
		}
	}

	_acb_vec_set(th, next, nbz * n);
	_acb_vec_clear(next, nbz * n);
}

int
acb_theta_ql_a0_steps_from_ctx(acb_ptr th, const acb_theta_ql_ctx_t ctx, slong nb_steps,
	slong s, slong prec, slong guard, acb_ql_worker_from_ctx_t worker)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	slong nbz = acb_theta_ql_ctx_nbz(ctx);
	slong nbr = (ctx->t_is_zero ? 1 : 2) * (ctx->z_is_zero ? 1 : 2);
	acb_theta_ql_ctx_t new_ctx;
	acb_ptr rts;
	arb_ptr d, d0;
	int res = 1;

	acb_theta_ql_ctx_init(new_ctx, g);
	rts = _acb_vec_init(nbr * n * nb_steps);
	d = _arb_vec_init(n);
	d0 = _arb_vec_init(n);

	acb_theta_ql_ctx_copy(new_ctx, ctx);

	for (k = 0; (k < nb_steps) && res; k++)
	{
		/* Compute rts efficiently */
		res = acb_theta_ql_a0_get_roots(rts + nbr * n * k, new_ctx, guard);
		acb_theta_ql_ctx_dupl(new_ctx, new_ctx, prec);
	}

	if (res)
	{
		/* Starting values */
		if (s == 0)
		{
			res = acb_theta_ql_a0_naive_from_ctx(th, new_ctx, prec);
		}
		else
		{
			res = acb_theta_ql_a0_split_from_ctx(th, new_ctx, s, prec, worker);
		}
	}

	if (res)
	{
		/* Duplication steps */
		for (k = nb_steps - 1; k >= 0; k--)
		{
			_arb_vec_scalar_mul_2exp_si(d0, &ctx->dists_a0, n, k);
			_arb_vec_scalar_mul_2exp_si(d, &ctx->dists_a0 + n, n, k);
			acb_theta_ql_a0_step(th, rts, d0, d, ctx, k == 0, prec);
		}
		if (!ctx->z_is_real)
		{
			_acb_vec_scalar_div_arb(th + nbt * n, th + nbt * n, nbt * n,
				acb_theta_ql_ctx_c(ctx), prec);
		}
	}

	acb_theta_ql_ctx_clear(new_ctx, g);
	_acb_vec_clear(rts, nbz * n * nb_steps);
}
