/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

int
acb_theta_ql_all_new(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_mat_t new_tau;
    acb_ptr new_z, new_th;
    slong j, k;
    int add_zero = 0;
    int res;

    if (nb <= 0)
    {
        return 1;
    }
    else if (!_acb_vec_is_zero(zs, g))
    {
        add_zero = 1;
    }

    acb_mat_init(new_tau, g, g);
    new_z = _acb_vec_init((nb + add_zero) * g);

    /* Strip tau, z of error bounds */
    for (j = 0; j < g; j++)
    {
        for (k = j; k < g; k++)
        {
            acb_get_mid(acb_mat_entry(new_tau, j, k), acb_mat_entry(tau, j, k));
            acb_set(acb_mat_entry(new_tau, k, j), acb_mat_entry(new_tau, j, k));
        }
    }
    for (j = 0; j < nb * g; j++)
    {
        acb_get_mid(&new_z[j], &zs[g * add_zero + j]);
    }

    /* Call ql_exact, with an extra duplication step if sqr is set */
    if (sqr)
    {
        new_th = _acb_vec_init((nb + add_zero) * n);

        /* Duplication */
        acb_mat_scalar_mul_2exp_si(new_tau, new_tau, 1);
        _acb_vec_scalar_mul_2exp_si(new_z, new_z, (nb + add_zero) * g, 1);
        res = acb_theta_ql_exact(new_th, new_z, nb + add_zero, new_tau, 0, 0, prec);
        if (res)
        {

        }
        acb_mat_scalar_mul_2exp_si(new_tau, new_tau, -1);
        _acb_vec_scalar_mul_2exp_si(new_z, new_z, (nb + add_zero) * g, 1);

        _acb_vec_clear(new_th, (nb + add_zero) * n);
    }
    else
    {
        new_th = _acb_vec_init((nb + add_zero) * n * n);

        res = acb_theta_ql_exact(new_th, new_z, nb + add_zero, tau, 1, 0, prec);
        if (res)
        {
            _acb_vec_set(th, new_th + add_zero * n * n, nb * n * n);
        }

        _acb_vec_clear(new_th, (nb + add_zero) * n * n);
    }

    /* Add error bounds and clear */

    return res;
}
