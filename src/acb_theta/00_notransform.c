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
acb_theta_00_notransform_gen(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
	slong lp = ACB_THETA_LOW_PREC;
	acb_theta_ctx_t ctx;
	acb_theta_eld_t E;
	arb_ptr v, a, c;
	arf_t R2, eps;
	arb_t err;
	slong j;
	int res;

	FLINT_ASSERT(nb >= 0);
	if (nb == 0) return;

	acb_theta_ctx_init(ctx, nb, g);
	acb_theta_eld_init(E, g, g);
	v = _arb_vec_init(g);
	a = _arb_vec_init(g);
	c = _arb_vec_init(nb);
	arf_init(R2);
	arf_init(eps);
	arb_init(err);

	acb_theta_ctx_set_tau(ctx, tau, prec);
	for (j = 0; j < nb; j++)
	{
		acb_theta_ctx_set_z(ctx, a, &c[j], zs + j * g, j, prec);
	}

	_arb_vec_set(v, acb_theta_ctx_vs(ctx), g); /* since nb >= 1 */
	for (j = 1; j < nb; j++)
	{
		_arb_vec_union(v, v, acb_theta_ctx_vs(ctx), g);
	}
	acb_theta_naive_radius(R2, eps, acb_theta_ctx_cho(ctx), 0, prec);
	res = acb_theta_eld_set(E, acb_theta_ctx_cho(ctx), R2, v);

	if (res)
	{
		acb_theta_naive_worker(th, 1, acb_theta_ctx_exp_zs(ctx), acb_theta_exp_zs_inv(ctx),
			nb, acb_theta_ctx_exp_tau(ctx), acb_theta_ctx_exp_tau_inv(ctx), E, 0, prec,
			acb_theta_naive_00_worker);
		for (j = 0; j < nb; j++)
		{
			arb_set_arf(err, eps);
			arb_div(err, err, &acb_theta_ctx_cs(ctx)[j], lp);
			acb_add_error_arb(&th[j], err);
			acb_mul(&th[j], &th[j], &c[j], prec);
		}
	}
	else
	{
		for (j = 0; j < nb; j++)
        {
            acb_indeterminate(&th[j]);
        }
    }

	acb_theta_ctx_clear(ctx);
    acb_theta_eld_clear(E);
	_arb_vec_clear(v, g);
	_arb_vec_clear(a, g);
	_arb_vec_clear(c, nb);
    arf_clear(R2);
    arf_clear(eps);
	arb_clear(err);
}

static void
acb_theta_00_notransform_g1(acb_ptr th, acb_srcptr zs, slong nb, const acb_t tau, slong prec)
{
	slong g = 1;
	acb_theta_ctx_t ctx;
    acb_ptr res;
	arb_t a;
	acb_t c;
    int w_is_unit;
    slong j;

	FLINT_ASSERT(nb >= 0);
	if (nb == 0) return;

	acb_theta_ctx_init(ctx, 1, g); /* don't need to store every z */
	res = _acb_vec_init(4);
	arb_init(a);
	acb_init(c);

	acb_theta_ctx_set_tau(ctx, tau, prec);
	for (j = 0; j < nb; j++)
	{
		acb_theta_ctx_set_z(ctx, a, c, &zs[j], 1, prec);
		acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
			&acb_theta_ctx_exp_zs(ctx)[0], acb_is_real(&zs[j]),
			acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), 1, prec);
		acb_mul(&th[j], &res[2], c, prec);
	}

	acb_theta_ctx_clear(ctx);
    _acb_vec_clear(res, 4);
	arb_clear(a);
	acb_clear(c);
}

void
acb_theta_00_notransform(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
	slong g = acb_mat_nrows(tau);

	if (g == 1)
	{
		acb_theta_00_reduced_z_g1(th, zs, nb, acb_mat_entry(tau, 0, 0), prec);
	}
	else
	{
		acb_theta_00_reduced_z(th, zs, nb, tau, prec);
	}
}
