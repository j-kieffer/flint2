/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

int acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int all, int sqr, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau;
    acb_ptr rts, th_init, t;
    arb_ptr distances;
    slong split, nb_steps, guard;
    slong * easy_steps;
    slong lp = ACB_THETA_LOW_PREC;
    slong hp;
    slong j;
    int res;

    if (nb <= 0)
    {
        return 1;
    }

    acb_theta_ctx_tau_init(ctx_tau, g);
    acb_theta_ctx_tau_set(ctx_tau, tau, lp);
    nb_steps = acb_theta_ql_nb_steps(&split, ctx_tau, lp);

    if (nb_steps == 0)
    {
        /* fall back to summation algorithm, return 0? (no error handling in that case) */
    }

    rts = _acb_vec_init(3 * n * nb * nb_steps);
    t = _acb_vec_init(g);
    th_init = _acb_vec_init(3 * n * nb);
    distances = _arb_vec_init(n * nb);
    easy_steps = flint_malloc(nb * sizeof(slong));

    /* Setup */
    for (j = 0; j < nb; j++)
    {
        acb_theta_dist_a0(distances + j * n, zs + j * n, tau, lp);
    }
    res = acb_theta_ql_setup(rts, t, &guard, easy_steps, zs, nb, tau, distances,
        nb_steps, all, sqr, prec);
    hp = prec + nb_steps * guard;

    /* Set th_init */
    if (res && (split == 0))
    {
        /* Use sum_a0_tilde */
        acb_theta_ctx_tau_t ctx;
        acb_theta_ctx_z_struct * aux;
        acb_theta_ctx_z_t ctxt;
        acb_mat_t new_tau;
        acb_ptr new_z;
        arb_ptr d;

        acb_theta_ctx_tau_init(ctx, g);
        aux = acb_theta_ctx_z_vec_init(3, g);
        acb_mat_init(new_tau, g, g);
        new_z = _acb_vec_init(g);
        d = _arb_vec_init(n);
        if (easy_steps[0] < nb_steps)
        {
            acb_theta_ctx_z_init(ctxt, g);
        }

        acb_mat_scalar_mul_2exp_si(new_tau, tau, nb_steps);
        acb_theta_ctx_tau_set(ctx_tau, new_tau, hp);
        if (easy_steps[0] < nb_steps)
        {
            _acb_vec_scalar_mul_2exp_si(new_z, t, g, nb_steps);
            acb_theta_ctx_z_set(ctxt, new_z, ctx_tau, hp);
        }

        for (j = 0; j < nb; j++)
        {
            _acb_vec_scalar_mul_2exp_si(new_z, zs + j * g, g, nb_steps);
            acb_theta_ctx_z_set(&aux[0], new_z, ctx_tau, hp);
            _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, nb_steps);
            if (easy_steps[j] == nb_steps)
            {
                acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 1, ctx_tau, d, hp);
            }
            else
            {
                acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, hp);
                acb_theta_ctx_z_add_real(&aux[2], &aux[1], ctxt, hp);
                acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 3, ctx_tau, d, hp);
            }
        }

        /* clear; don't forget ctxt */
    }
    else if (res)
    {
        /* Use splitting strategy */
    }

    if (res)
    {
        acb_theta_ql_steps(th, th_init, rts, nb, nb_steps, distances, easy_steps, g, hp);
    }

    /* clear */

    return res;
}
