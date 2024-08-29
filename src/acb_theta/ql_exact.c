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

static int
acb_theta_ql_exact_lower_dim(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    slong s, arb_srcptr distances, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong nba = 1 << (g - s);
    slong n0 = 1 << s;
    slong n = 1 << g;
    acb_ptr * new_zs;
    acb_ptr * cofactors;
    slong * nb_pts;
    slong * fullprec;
    arf_struct * err;
    acb_ptr th0, z0s;
    slong nb0;
    acb_mat_t tau0;
    ulong a, a0;
    slong j, k, index;
    int res = 1;

    new_zs = flint_malloc(nb * nba * sizeof(acb_ptr));
    cofactors = flint_malloc(nb * nba * sizeof(acb_ptr));
    nb_pts = flint_malloc(nb * nba * sizeof(slong));
    err = flint_malloc(nb * nba * sizeof(arf_struct));
    for (j = 0; j < nb * nba; j++)
    {
        arf_init(&err[j]);
    }
    fullprec = flint_malloc(nb * nba * sizeof(arf_struct));
    acb_mat_window_init(tau0, tau, 0, 0, s, s);

    nb0 = 0;
    for (j = 0; j < nb; j++)
    {
        for (a = 0; a < nba; a++)
        {
            if (res)
            {
                res = acb_theta_ql_lower_dim(&new_zs[j * nba + a], &cofactors[j * nba + a],
                    &nb_pts[j * nba + a], &err[j * nba + a], &fullprec[j * nba + a],
                    zs + j * g, tau, distances + j * n, s, a, prec);
                nb0 += nb_pts[j * nba + a];
            }
            else
            {
                /* Initialize with length 0 to be able to free later. */
                new_zs[j * nba + a] = _acb_vec_init(0);
                cofactors[j * nba + a] = _acb_vec_init(0);
                nb_pts[j * nba + a] = 0;
                fullprec[j * nba + a] = 0;
            }
        }
    }

    if (!res)
    {
        nb0 = 0;
    }
    z0s = _acb_vec_init(nb0 * s);
    th0 = _acb_vec_init(nb0 * n0);

    if (res)
    {
        /* Put everything together */
        index = 0;
        for (j = 0; j < nb; j++)
        {
            for (a = 0; a < nba; a++)
            {
                _acb_vec_set(z0s + index * s, new_zs[j * nba + a], nb_pts[j * nba + a] * s);
                index += nb_pts[j * nba + a];
            }
        }

        /* Call acb_theta_ql_exact in dimension s */
        res = acb_theta_ql_exact(th0, z0s, nb0, tau0, 0, shifted_prec, prec);
    }

    if (res)
    {
        /* Rescale using cofactors and sum into th */
        index = 0;
        _acb_vec_zero(th, nb * n);
        for (j = 0; j < nb; j++)
        {
            for (a = 0; a < nba; a++)
            {
                for (k = 0; k < nb_pts[j * nba + a]; k++)
                {
                    _acb_vec_scalar_mul(th0 + index * n0, th0 + index * n0,
                        n0, &cofactors[j * nba + a][k], prec);
                    for (a0 = 0; a0 < n0; a0++)
                    {
                        acb_add(&th[j * n + (a0 << (g - s)) + a],
                            &th[j * n + (a0 << (g - s)) + a],
                            &th0[index * n0 + a0], fullprec[j * nba + a]);
                    }
                    index += 1;
                }
                for (a0 = 0; a0 < n0; a0++)
                {
                    acb_add_error_arf(&th[j * n + (a0 << (g - s)) + a],
                        &err[j * nba + a]);
                }
            }
        }
    }

    /* Clear */
    for (j = 0; j < nb; j++)
    {
        for (a = 0; a < nba; a++)
        {
            _acb_vec_clear(new_zs[j * nba + a], nb_pts[j * nba + a] * s);
            _acb_vec_clear(cofactors[j * nba + a], nb_pts[j * nba + a]);
        }
    }
    flint_free(new_zs);
    flint_free(cofactors);
    flint_free(nb_pts);
    for (j = 0; j < nb * nba; j++)
    {
        arf_clear(&err[j]);
    }
    flint_free(err);
    flint_free(fullprec);
    acb_mat_window_clear(tau0);
    _acb_vec_clear(z0s, nb0 * s);
    _acb_vec_clear(th0, nb0 * n0);
    return res;
}

static void
acb_theta_ql_exact_sum(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    arb_srcptr distances, int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau;
    acb_theta_ctx_z_struct * vec;
    slong j;

    FLINT_ASSERT(nb >= 1);

    acb_theta_ctx_tau_init(ctx_tau, g);
    vec = acb_theta_ctx_z_vec_init(nb, g);

    acb_theta_ctx_tau_set(ctx_tau, tau, prec);
    for (j = 0; j < nb; j++)
    {
        acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
    }

    if (shifted_prec)
    {
        for (j = 0; j < nb; j++)
        {
            acb_theta_sum_a0_tilde(th + j * n, &vec[j], 1, ctx_tau, distances + j * n, prec);
        }
    }
    else
    {
        /* distances are all set to zero */
        acb_theta_sum_a0_tilde(th, vec, nb, ctx_tau, distances, prec);
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    acb_theta_ctx_z_vec_clear(vec, nb);
}

static int
acb_theta_ql_exact_steps(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    slong split, arb_srcptr distances, slong nb_steps, int all, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_ptr rts, th_init, t;
    slong * easy_steps;
    acb_mat_t new_tau;
    slong guard, hp;
    slong j, k;
    int res;

    rts = _acb_vec_init(3 * n * nb * nb_steps);
    t = _acb_vec_init(g);
    th_init = _acb_vec_init(3 * n * nb);
    easy_steps = flint_malloc(nb * sizeof(slong));
    acb_mat_init(new_tau, g, g);

    res = acb_theta_ql_setup(rts, t, &guard, easy_steps, zs, nb, tau, distances,
        nb_steps, all, prec);
    hp = prec + nb_steps * guard;
    acb_mat_scalar_mul_2exp_si(new_tau, tau, nb_steps);

    if (res && (split > 0))
    {
        /* Set list of zs for which we want theta to be computed. */
        /* We don't take advantage of the fact that some thetas differ by a real value */
        acb_ptr new_th, new_z;
        arb_ptr new_distances;
        slong new_nb, add;

        /* Count number of vectors */
        new_nb = 0;
        for (j = 0; j < nb; j++)
        {
            if (j == 0 && easy_steps[j] < nb_steps)
            {
                new_nb += 3;
            }
            else if (easy_steps[j] < nb_steps)
            {
                new_nb += 2;
            }
            else
            {
                new_nb += 1;
            }
        }

        new_th = _acb_vec_init(n * new_nb);
        new_z = _acb_vec_init(g * new_nb);
        new_distances = _arb_vec_init(n * new_nb);

        /* Set up input of ql_exact_lower_dim */
        new_nb = 0;
        for (j = 0; j < nb; j++)
        {
            if (j == 0 && easy_steps[j] < nb_steps)
            {
                _acb_vec_zero(new_z + new_nb * g, g);
                _acb_vec_set(new_z + (new_nb + 1) * g, t, g);
                _acb_vec_scalar_mul_2exp_si(new_z + (new_nb + 2) * g, t, g, 1);
                add = 3;
            }
            else if (easy_steps[j] < nb_steps)
            {
                _acb_vec_add(new_z + new_nb * g, zs + j * g, t, g, hp);
                _acb_vec_add(new_z + (new_nb + 1) * g, new_z + (new_nb + 1) * g, t, g, hp);
                add = 2;
            }
            else
            {
                _acb_vec_set(new_z + new_nb * g, zs + j * g, g);
                add = 1;
            }
            for (k = 0; k < add; k++)
            {
                _arb_vec_set(new_distances + (new_nb + k) * n, distances + j * n, n);
            }
            new_nb += add;
        }

        _acb_vec_scalar_mul_2exp_si(new_z, new_z, new_nb * g, nb_steps);
        _arb_vec_scalar_mul_2exp_si(new_distances, new_distances, new_nb * n, nb_steps);

        /* Recursive call */
        res = acb_theta_ql_exact_lower_dim(new_th, new_z, new_nb, new_tau, split,
            new_distances, 0, 1, hp);

        /* Set up th_init from computed data */
        new_nb = 0;
        for (j = 0; j < nb; j++)
        {
            if (j == 0 && easy_steps[j] < nb_steps)
            {
                _acb_vec_set(th_init, new_th, 3 * n);
                new_nb += 3;
            }
            else if (easy_steps[j] < nb_steps)
            {
                _acb_vec_set(th_init + 3 * j * n + n, new_th + new_nb * n, 2 * n);
                new_nb += 2;
            }
            else
            {
                _acb_vec_set(th_init + 3 * j * n, new_th + new_nb * n, n);
                new_nb += 1;
            }
        }

        _acb_vec_clear(new_th, new_nb * n);
        _acb_vec_clear(new_z, new_nb * g);
        _arb_vec_clear(new_distances, new_nb * n);
    }
    else if (res) /* split = 0; set up th_init efficiently using sum_a0_tilde */
    {
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * aux;
        acb_theta_ctx_z_t ctxt;
        acb_mat_t new_tau;
        acb_ptr new_z;
        arb_ptr d;

        acb_theta_ctx_tau_init(ctx_tau, g);
        aux = acb_theta_ctx_z_vec_init(3, g);
        new_z = _acb_vec_init(g);
        d = _arb_vec_init(n);
        if (easy_steps[0] < nb_steps)
        {
            acb_theta_ctx_z_init(ctxt, g);
        }

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
                if (j == 0)
                {
                    acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 3, ctx_tau, d, hp);
                }
                else
                {
                    acb_theta_sum_a0_tilde(th_init + 3 * n * j + n, aux + 1, 2, ctx_tau, d, hp);
                }
            }
        }

        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(aux, 3);
        _acb_vec_clear(new_z, g);
        _arb_vec_clear(d, n);
        if (easy_steps[0] < nb_steps)
        {
            acb_theta_ctx_z_clear(ctxt);
        }
    }

    if (res) /* th_init is set; now perform steps */
    {
        acb_theta_ql_steps(th, th_init, rts, nb, nb_steps, distances, easy_steps, g, hp);
    }

    _acb_vec_clear(rts, 3 * n * nb * nb_steps);
    _acb_vec_clear(t, g);
    _acb_vec_clear(th_init, 3 * n * nb);
    flint_free(easy_steps);
    acb_mat_clear(new_tau);
    return res;
}

int acb_theta_ql_exact(acb_ptr th, acb_srcptr zs, slong nb, const acb_mat_t tau,
    int all, int shifted_prec, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    acb_theta_ctx_tau_t ctx_tau;
    arb_ptr distances;
    slong split, nb_steps;
    slong lp = ACB_THETA_LOW_PREC;
    slong j;
    int res = 1;

    if (nb <= 0)
    {
        return 1;
    }

    acb_theta_ctx_tau_init(ctx_tau, g);
    distances = _arb_vec_init(n * nb);

    acb_theta_ctx_tau_set(ctx_tau, tau, lp);
    nb_steps = acb_theta_ql_nb_steps(&split, ctx_tau, lp);
    if (nb_steps > 0 || shifted_prec)
    {
        for (j = 0; j < nb; j++)
        {
            acb_theta_dist_a0(distances + j * n, zs + j * n, tau, lp);
        }
    }

    flint_printf("(ql_exact) using nb_steps = %wd, split = %wd\n", nb_steps, split);

    if (nb_steps == 0 && split == 0)
    {
        acb_theta_ql_exact_sum(th, zs, nb, tau, distances, all, shifted_prec, prec);
    }
    else if (nb_steps == 0 && split > 0)
    {
        res = acb_theta_ql_exact_lower_dim(th, zs, nb, tau, split, distances, all, shifted_prec, prec);
    }
    else
    {
        res = acb_theta_ql_exact_steps(th, zs, nb, tau, split, distances, nb_steps, all, prec);
    }

    acb_theta_ctx_tau_clear(ctx_tau);
    _arb_vec_clear(distances, n * nb);
    return res;
}
