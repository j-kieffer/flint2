/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

int acb_theta_naive_a0(acb_ptr th, acb_srcptr exp_zs, acb_srcptr exp_zs_inv, slong nb,
	const acb_mat_t exp_tau_div_4, const acb_mat_t exp_tau, const acb_mat_t exp_tau_inv,
	acb_srcptr v, arb_srcptr d, const arb_mat_t C, int all, slong prec)
{
	slong g = acb_mat_nrows(exp_tau);
	slong n  = 1 << g;
	arb_ptr new_v;
	acb_ptr cs, aux, new_exp_zs, new_exp_zs_inv;
	arf_t R2, eps;
	acb_theta_eld_t E;
	slong new_prec;
	slong a;
	int success = 1;

	new_v = _arb_vec_init(g);
	arf_init(R2);
	arf_init(eps);
	acb_theta_eld_init(E, g);
	cs = _acb_vec_init(nb);
	aux = _acb_vec_init(nb * (all ? n : 1));
	new_exp_zs = _acb_vec_init(nb * g);
	new_exp_zs_inv = _acb_vec_init(nb * g);

	for (a = 0; (a < n) && success; a++)
	{
		acb_theta_char_get_arb(new_v, a, g);
		_arb_vec_neg(new_v, new_v, g);
		_arb_vec_add(new_v, new_v, v, g, prec);
		new_prec = prec + acb_theta_dist_addprec(&d[a]);
		acb_theta_naive_radius(R2, eps, C, 0, new_prec);
		success = acb_theta_eld_set(E, C, R2, new_v);

		if (!success) break;

		/* Translate exponentials: todo: a lot can be factored here */
		for (j = 0; j < nb; j++)
		{
			acb_theta_naive_exp_translate(&cs[j], new_exp_zs + j * g,
				exp_zs + j * g, exp_tau, a, prec);
		}
		for (j = 0; j < g * nb; j++)
		{
			/* Inversions are OK here. */
			acb_inv(&new_exp_zs_inv[j], &new_exp_zs[j], prec);
		}

		/* Call worker */
		if (all)
		{
			acb_theta_naive_worker(aux, n, new_exp_zs, new_exp_zs_inv, nb, exp_tau,
				exp_tau_inv, E, 0, new_prec, acb_theta_naive_0b_worker);
			for (j = 0; j < nb; j++)
			{
				_acb_vec_scalar_mul(th + j * n * n + a * n, aux + j * n, n, &cs[j], prec);
			}
		}
		else
		{
			acb_theta_naive_worker(aux, 1, new_exp_zs, new_exp_zs_inv, nb, exp_tau,
				exp_tau_inv, E, 0, new_prec, acb_theta_naive_00_worker);
			for (j = 0; j < nb; j++)
			{
				acb_mul(&th[j * n + a], &aux[j], &cs[j], prec);
			}
		}
	}
	/* Todo: add errors! */

	_acb_vec_clear(new_v, g);
	arf_clear(R2);
	arf_clear(eps);
	acb_theta_eld_clear(E);
	_acb_vec_clear(cs, nb);
	_acb_vec_clear(aux, nb * (all ? n : 1));
	_acb_vec_clear(new_exp_zs, nb * g);
	return success;
}
