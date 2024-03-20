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
acb_theta_ql_a0_naive_from_ctx_gen(cb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	arb_ptr new_v;
	acb_theta_eld_t E;
	arf_t R2, eps;
	acb_ptr new_exp_zs, aux, cs;
	slong new_prec;
	slong a, j;
	int success = 1;

	new_v = _arb_vec_init(g);
	arf_init(R2);
	arf_init(eps);
	new_exp_zs = _acb_vec_init(nbt * g);
	aux = _acb_vec_init(nbt); /* todo: could do 6 at once if z is real */
	cs = _acb_vec_init(nbt);

	/* Get theta_a0 at 0, t, 2t */
	for (a = 0; (a < n) && success; a++)
	{
		/* Compute ellipsoid */
		acb_theta_char_get_arb(new_v, a, g);
		_arb_vec_neg(new_v, new_v, g);
		_arb_vec_add(new_v, new_v, ctx->vs, g, prec);
		new_prec = prec + acb_theta_dist_addprec(&ctx->dists_a0[a]);
		acb_theta_naive_radius(R2, eps, &ctx->C, 0, new_prec);
		success = acb_theta_eld_set(E, &ctx->C, R2, new_v);

		if (!success) break;

		/* Translate exponentials */
		for (j = 0; j < nbt; j++)
		{
			acb_theta_naive_exp_translate(&cs[j], new_exp_zs + j * g, ctx->exp_zs + j * g,
				&ctx->exp_tau, a, prec);
		}

		/* Call worker */
		acb_theta_naive_worker(aux, 1, new_exp_zs, nbt, &ctx->exp_tau, E, 0, new_prec,
			acb_theta_naive_00_worker);
		for (j = 0; j < nbt; j++)
		{
			acb_mul(&th[j * n + a], &aux[j], &cs[j], prec);
		}
	}

	/* Get theta_a0 at z, z + t, z + 2t */
	if (!ctx->z_is_zero)
	{
		for (a = 0; (a < n) && success; a++)
		{
			/* Compute ellipsoid */
			acb_theta_char_get_arb(new_v, a, g);
			_arb_vec_neg(new_v, new_v, g);
			_arb_vec_add(new_v, new_v, ctx->vs + g, g, prec);
			new_prec = prec + acb_theta_dist_addprec(&ctx->dists_a0[n + a]);
			acb_theta_naive_radius(R2, eps, &ctx->C, 0, new_prec);
			success = acb_theta_eld_set(E, &ctx->C, R2, new_v);

			if (!success) break;

			/* Translate exponentials */
			for (j = 0; j < nbt; j++)
			{
				acb_theta_naive_exp_translate(&cs[j], new_exp_zs + j * g, ctx->exp_zs + (nbt + j) * g,
					exp_tau, a, prec);
			}

			acb_theta_naive_worker_from_ctx(aux, 1, new_exp_zs,  nbt,
					&ctx->exp_tau, E, 0, new_prec, acb_theta_naive_00_worker);
			for (j = 0; j < nbt; j++)
			{
				acb_mul(&th[(nbt + j) * n + a], &aux[j], &cs[j], prec);
			}
		}
	}

	_acb_vec_clear(new_v, g);
	arf_clear(R2);
	arf_clear(eps);
	_acb_vec_clear(new_exp_zs, nbt * g);
	_acb_vec_init(aux, nbt);
	_acb_vec_clear(cs, nbt);
}

void
acb_theta_ql_a0_naive_from_ctx(acb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);

	if (g == 1)
	{
		acb_ptr aux;
		slong j;
		aux = _acb_vec_init(4 * nbz);

		acb_theta_ql_all_naive_from_ctx(aux, ctx, prec);
		for (j = 0; j < nbz; j++)
		{
			acb_set(&th[2 * j], &aux[4 * j]);
			acb_set(&th[2 * j + 1], &aux[4 * j + 2]);
		}

		_acb_vec_clear(aux);
	}
	else
	{
		acb_theta_ql_a0_naive_from_ctx_gen(th, ctx, prec);
	}
}
