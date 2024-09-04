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
    fmpz_mat_t mat;
    acb_mat_t new_tau, N, c;
    acb_ptr new_zs, exps, aux, units;
    acb_t s;
    ulong * image_ab;
    slong * e;
    slong kappa;
    slong j, ab;
    int res;

    if (nb <= 0)
    {
        return;
    }

    fmpz_mat_init(mat, 2 * g, 2 * g);
    acb_mat_init(new_tau, g, g);
    acb_mat_init(N, g, g);
    acb_mat_init(c, g, g);
    new_zs = _acb_vec_init(nb * g);
    exps = _acb_vec_init(nb);
    aux = _acb_vec_init(n2 * nb);
    units = _acb_vec_init(8);
    image_ab = flint_malloc(n2 * sizeof(ulong));
    e = flint_malloc(n2 * sizeof(slong));
    acb_init(s);

    res = acb_theta_reduce_tau(new_zs, new_tau, mat, N, exps, zs, nb, tau, prec);

    if (res)
    {
        /* todo: reduce z here */

        _acb_vec_unit_roots(units, 8, 8, prec);
        if (sqr)
        {
            kappa = acb_theta_transform_kappa2(mat);
            acb_siegel_cocycle(c, mat, new_tau, prec);
            acb_mat_det(s, c, prec);
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
            if (sqr)
            {
                acb_sqr(&exps[j], &exps[j], prec);
            }
            for (ab = 0; ab < n2; ab++)
            {
                acb_mul(&th[j * n2 + ab], &aux[j * n2 + image_ab[ab]], &exps[j], prec);
                acb_mul(&th[j * n2 + ab], &th[j * n2 + ab], s, prec);
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
    acb_mat_clear(N);
    _acb_vec_clear(new_zs, nb * g);
    _acb_vec_clear(exps, nb);
    _acb_vec_clear(aux, n2 * nb);
    _acb_vec_clear(units, 8);
    acb_clear(s);
    flint_free(e);
    flint_free(image_ab);
}
