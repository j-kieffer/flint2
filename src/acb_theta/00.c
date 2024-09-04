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

void
acb_theta_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    fmpz_mat_t mat;
    acb_mat_t new_tau, N;
    acb_ptr new_zs, exps;
    acb_ptr aux, units;
    acb_t s;
    slong kappa, e, ab;
    slong j;
    int res;

    if (nb <= 0)
    {
        return;
    }

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(N, g, g);
    new_zs = _acb_vec_init(nb * g);
    exps = _acb_vec_init(nb);
    aux = _acb_vec_init(nb);
    units = _acb_vec_init(8);
    acb_init(s);

    _acb_vec_unit_roots(units, 8, 8, prec);
    res = acb_theta_reduce_tau(new_zs, new_tau, mat, N, exps, zs, nb, tau, prec);

    if (res)
    {
        /* todo: reduce z here. */

        ab = acb_theta_transform_char(&e, mat, 0);
        kappa = acb_theta_transform_kappa(s, mat, new_tau, prec);

        acb_theta_one_notransform(aux, new_zs, nb, new_tau, ab, prec);

        for (j = 0; j < nb; j++)
        {
            acb_mul(&th[j], &aux[j], &exps[j], prec);
            acb_mul(&th[j], &th[j], s, prec);
            acb_mul(&th[j], &th[j], &units[(kappa + e) % 8], prec);
        }
    }
    else
    {
        /* todo: replace with upper bound */
        _acb_vec_indeterminate(th, nb);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(N);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(exps, nb);
    _acb_vec_clear(aux, nb);
    _acb_vec_clear(units, 8);
    acb_clear(s);
}
