/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_jet_00(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    slong ord, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nbth = acb_theta_jet_nb(ord, g);
    fmpz_mat_t mat;
    acb_mat_t new_tau, N, ct;
    acb_ptr new_zs, aux, exps, units;
    acb_t s, t;
    ulong image_ab;
    slong kappa, e, j;
    int res;

    if (nb <= 0)
    {
        return;
    }

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(N, g, g);
    acb_mat_init(ct, g, g);
    new_zs = _acb_vec_init(nb * g);
    exps = _acb_vec_init(nb);
    aux = _acb_vec_init(nb * nbth);
    units = _acb_vec_init(8);
    acb_init(s);
    acb_init(t);

    res = acb_theta_reduce_tau(new_zs, new_tau, mat, N, ct, exps, zs, nb, tau, prec);

    if (res)
    {
        /* todo: reduce z here */

        _acb_vec_unit_roots(units, 8, 8, prec);
        kappa = acb_theta_transform_kappa(s, mat, new_tau, prec);
        image_ab = acb_theta_transform_char(&e, mat, 0);

        acb_theta_jet_one_notransform(aux, new_zs, nb, new_tau, ord, image_ab, prec);

        for (j = 0; j < nb; j++)
        {
            acb_mul(t, s, &units[(kappa + e) % 8], prec);
            _acb_vec_scalar_mul(th + j * nbth, aux + j * nbth, nbth, t, prec);
            acb_theta_jet_compose(th + j * nbth, th + j * nbth, ct, ord, prec);
            acb_theta_jet_exp_qf(aux + j * nbth, zs + j * g, N, ord, prec);
            acb_theta_jet_mul(th + j * nbth, th + j * nbth, aux + j * nbth, ord, g, prec);
        }
    }
    else
    {
        _acb_vec_indeterminate(th, nb * nbth);
    }

    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(N);
    acb_mat_clear(ct);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(aux, nb * nbth);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    acb_clear(t);
}
