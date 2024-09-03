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
acb_theta_all(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n2 = 1 << (2 * g);
    fmpz_mat_t mat, gamma;
    acb_mat_t new_tau, c, cinv, N;
    acb_ptr new_zs, y, aux, units;
    acb_t s, t;
    ulong * image_ab;
    slong * e;
    slong kappa;
    slong j, ab;

    if (nb <= 0)
    {
        return;
    }

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(c, g, g);
    acb_mat_init(cinv, g, g);
    acb_mat_init(N, g, g);
    new_zs = _acb_vec_init(nb * g);
    y = _acb_vec_init(g);
    aux = _acb_vec_init(n2 * nb);
    units = _acb_vec_init(8);
    image_ab = flint_malloc(n2 * sizeof(ulong));
    e = flint_malloc(n2 * sizeof(slong));
    acb_init(s);
    acb_init(t);

    acb_siegel_reduce(mat, tau, prec);
    acb_siegel_transform_cocycle_inv(new_tau, c, cinv, mat, tau, prec);

    acb_mat_transpose(cinv, cinv);
    for (j = 0; j < nb; j++)
    {
        acb_mat_vector_mul_col(new_zs + j * g, cinv, zs + j * g, prec);
    }
    _acb_vec_unit_roots(units, 8, 8, prec);

    if (acb_siegel_is_reduced(new_tau, -10, prec))
    {
        sp2gz_inv(mat, mat);

        fmpz_mat_window_init(gamma, mat, g, 0, 2 * g, g);
        acb_mat_set_fmpz_mat(N, gamma);
        acb_mat_mul(N, c, N, prec);
        fmpz_mat_window_clear(gamma);

        /* todo: reduce z here */

        if (sqr)
        {
            kappa = acb_theta_transform_kappa2(mat);
            acb_mat_det(s, cinv, prec);
        }
        else
        {
            kappa = acb_theta_transform_kappa(s, mat, new_tau, prec);
        }
        for (ab = 0; ab < n2; ab++)
        {
            /* todo: can do better than this loop */
            image_ab[ab] = acb_theta_transform_char(&e[ab], mat, ab);
        }

        acb_theta_all_notransform(aux, new_zs, nb, new_tau, sqr, prec);

        for (j = 0; j < nb; j++)
        {
            acb_mat_vector_mul_col(y, N, new_zs + j * g, prec);
            acb_dot(t, NULL, 0, new_zs + j * g, 1, y, 1, g, prec);
            if (sqr)
            {
                acb_mul_2exp_si(t, t, 1);
            }
            acb_exp_pi_i(t, t, prec);
            acb_mul(t, t, s, prec);
            for (ab = 0; ab < n2; ab++)
            {
                acb_mul(&th[j * n2 + ab], &aux[j * n2 + image_ab[ab]], t, prec);
                acb_mul(&th[j * n2 + ab], &th[j * n2 + ab],
                    &units[((sqr ? 2 : 1) * (kappa + e[ab])) % 8], prec);
            }
        }
    }
    else
    {
        _acb_vec_indeterminate(th, n2 * nb);
    }


    fmpz_mat_clear(mat);
    acb_mat_clear(new_tau);
    acb_mat_clear(c);
    acb_mat_clear(cinv);
    acb_mat_clear(N);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(y, g);
    _acb_vec_clear(aux, n2 * nb);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    acb_clear(t);
    flint_free(e);
    flint_free(image_ab);
}
