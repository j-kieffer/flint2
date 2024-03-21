/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void acb_theta_naive_a0_g1(acb_ptr th, const acb_t exp_z, int w_is_unit,
	const acb_t exp_tau_div_4, const acb_t exp_tau, const arb_t y, int all, slong prec)
{
	acb_ptr aux;

	aux = _acb_vec_init(4);

	acb_modular_theta_sum(&aux[3], &aux[2], &aux[0], &aux[1], exp_z, w_is_unit,
		acb_mat_entry(exp_tau, 0, 0), 1, prec);
	acb_neg(&aux[3], &aux[3]);

	if (all)
	{
		_acb_vec_set(th, aux, 4);
		acb_mul(&th[2], &th[2], exp_tau_div_4, prec);
		acb_mul(&th[3], &th[3], exp_tau_div_4, prec);
	}
	else
	{
		acb_set(&th[0], &aux[0]);
		acb_mul(&th[1], &aux[2], exp_tau_div_4, prec);
	}

	_acb_vec_clear(aux, 4);
}
