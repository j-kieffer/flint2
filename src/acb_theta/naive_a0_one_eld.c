/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int acb_theta_naive_a0_one_eld(acb_ptr th, acb_srcptr exp_zs, acb_srcptr exp_zs_inv, slong nb,
	const acb_mat_t exp_tau_div_4, const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv,
	acb_srcptr v, arb_srcptr d, const arb_mat_t C, int all, slong prec)
{
	slong g = acb_mat_nrows(exp_tau);
	slong n = 1 << g;
	arb_ptr new_v;
	acb_ptr cs, new_exp_zs, new_exp_zs_inv;
	arf_t R2, eps;
	acb_theta_eld_t E;
	slong new_prec;
	slong j, a;
	int success = 1;

	new_v = _arb_vec_init(g);
	arf_init(R2);
	arf_init(eps);
	acb_theta_eld_init(E, g);
	new_exp_zs = _acb_vec_init(nb * n * g);
	new_exp_zs_inv = _acb_vec_init(nb * n * g);

	for (a = 0; a < n; a++)
	{
		acb_theta_naive_exp_translate(&cs[a], exps_translate + a * g, exp_tau_div_4, a, prec);
	}
	for (j = 0; j < g; j++)
	{
		arb_unit_interval(&new_v[j]);
	}
	_arb_vec_neg(new_v, new_v, g);
	_arb_vec_scalar_mul_2exp_si(new_v, new_v, g, -1);
	_arb_vec_add(new_v, new_v, v, g, prec);

	new_prec = prec;
	for (a = 0; a < n; a++)
	{
		new_prec = FLINT_MAX(new_prec, prec + acb_theta_dist_addprec(&d[a]));
	}

	acb_theta_naive_radius(R2, eps, C, 0, new_prec);
	success = acb_theta_eld_set(E, C, R2, new_v);

	if (success)
	{
		if (all)
		{
			acb_theta_naive_worker(th, n, new_exp_zs, new_exp_zs_inv, nb * n, exp_tau,
				exp_tau_inv, E, 0, new_prec, acb_theta_naive_0b_worker);
			for (j = 0; j < nb; j++)
			{
				for (a = 0; a < n; a++)
				{
					_acb_vec_scalar_mul(th + j * n * n + a * n, th + j * n * n + a * n,
						n, &cs[a], prec);
				}
			}
		}
		else
		{
			acb_theta_naive_worker(th, 1, new_exp_zs, new_exp_zs_inv, nb * n, exp_tau,
				exp_tau_inv, E, 0, new_prec, acb_theta_naive_00_worker);
			for (j = 0; j < nb; j++)
			{
				for (a = 0; a < n; a++)
				{
					acb_mul(&th[j * n + a], &th[j * n + a], &cs[a], prec);
				}
			}
		}

	}

	_acb_vec_clear(new_v, g);
	arf_clear(R2);
	arf_clear(eps);
	acb_theta_eld_clear(E);
	_acb_vec_clear(new_exp_zs, nb * n * g);
	_acb_vec_clear(new_exp_zs_inv, nb * n * g);
}
