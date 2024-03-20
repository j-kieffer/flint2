/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

void
acb_theta_naive_worker_old(acb_ptr th, slong len, acb_srcptr zs, slong nb,
    const acb_mat_t tau, const acb_theta_eld_t E, slong ord, slong prec,
    acb_theta_naive_worker_t worker)
{
    slong g = acb_theta_eld_ambient_dim(E);
	acb_mat_t exp_tau;
	acb_ptr exp_zs;
	slong j;

	acb_mat_init(exp_tau, g, g);
	exp_zs = _acb_vec_init(g * nb);

	acb_theta_naive_exp_tau(exp_tau, tau, prec);
	for (j = 0; j < nb; j++)
	{
		acb_theta_naive_exp_z(exp_zs + j * g, zs + j * g, g, prec);
	}

	acb_theta_naive_worker(th, len, exp_zs, nb, exp_tau, E, ord, prec, worker);

	acb_mat_clear(exp_tau);
	_acb_vec_clear(exp_zs, g * nb);
}
