/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* copied from acb_theta_naive_00 */
static void
worker(acb_ptr th, acb_srcptr v1, acb_srcptr v2, const slong * precs, slong len,
    const acb_t cofactor, const slong * coords, slong ord, slong g, slong prec, slong fullprec)
{
    acb_t sum;

    acb_init(sum);

    acb_dot(sum, NULL, 0, v1, 1, v2, 1, len, prec);
    acb_addmul(th, sum, cofactor, fullprec);

    acb_clear(sum);
}

void
acb_theta_ql_a0_naive_from_ctx_gen(cb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
	slong g = acb_theta_ql_ctx_g(ctx);
	slong n = 1 << g;
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	arb_ptr new_v;
	acb_theta_eld_t E;
	arf_t R2, eps;
	acb_ptr aux, new_exp_zs, new_exp_zs_inv;
	slong new_prec;
	slong a, j;
	int success = 1;

	new_v = _arb_vec_init(g);
	arf_init(R2);
	arf_init(eps);
	new_exp_zs = _acb_vec_init(nbt * g);
	new_exp_zs_inv = _acb_vec_init(nbt * g);
	aux = _acb_vec_init(nbt); /* todo: could do 6 at once if z is real */

	/* Get theta_a0 at 0, t, 2t */
	for (a = 0; (a < n) && success; a++)
	{
		/* Figure out the right factors... */
		acb_theta_char_get_arb(new_v, a, g);
		arb_vec_add(new_v, new_v, ctx->vs, g, prec);
		new_prec = prec + acb_theta_dist_addprec(&ctx->dists_a0[a]);

		acb_theta_naive_radius(R2, eps, &ctx->C, 0, new_prec);
		success = acb_theta_eld_set(E, &ctx->C, R2, new_v);
		if (success)
		{
			acb_theta_naive_worker_from_ctx(aux, 1, new_exp_zs, new_exp_zs_inv, nbt,
				&ctx->exp_tau, &ctx->exp_tau_inv, E, 0, new_prec, worker);
			for (j = 0; j < nbt; j++)
			{
				acb_set(&th[j * n + a], &aux[j]);
				/* figure out the right cofactors... */
			}
		}
	}

	/* Get theta_a0 at z, z + t, z + 2t */
	if (!ctx->z_is_zero)
	{
		for (a = 0; (a < n) && success; a++)
		{
			acb_theta_char_get_arb(new_v, a, g);
			arb_vec_add(new_v, new_v, ctx->vs + g, g, prec);
			new_prec = prec + acb_theta_dist_addprec(&ctx->dists_a0[n + a]);

			acb_theta_naive_radius(R2, eps, &ctx->C, 0, new_prec);
			success = acb_theta_eld_set(E, &ctx->C, R2, new_v);
			if (success)
			{
				acb_theta_naive_worker_from_ctx(aux, 1, new_exp_zs, new_exp_zs_inv, nbt,
					&ctx->exp_tau, &ctx->exp_tau_inv, E, 0, new_prec, worker);
				for (j = 0; j < nbt; j++)
				{
					acb_set(&th[(nbt + j) * n + a], &aux[j]);
					/* figure out the right cofactors */
				}
			}
		}

	}
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
