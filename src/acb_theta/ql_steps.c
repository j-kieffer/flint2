/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "acb_theta.h"

void acb_theta_ql_step_1(acb_ptr res, acb_srcptr th0, acb_srcptr th, acb_srcptr rts,
    arb_srcptr d0, arb_srcptr d, slong g, slong prec)
{
    slong n = 1 << g;

    acb_theta_agm_mul_tight(res, th0, th, d0, d, g, prec);
    _acb_vec_scalar_mul_2exp_si(res, res, n, g);
    acb_theta_agm_sqrt(res, res, rts, n, prec);
}

void acb_theta_ql_steps(acb_ptr th, acb_ptr th_init, acb_srcptr rts,
    slong nb, slong nb_steps, arb_srcptr distances, const slong * easy_steps,
    slong g, slong prec)
{
    slong n = 1 << g;
    acb_ptr th_next, th_temp;
    arb_ptr d0, d;
    slong j, k, a;

    th_next = _acb_vec_init(3 * n * nb);
    d = _arb_vec_init(n);
    d0 = _arb_vec_init(n);

    for (k = nb_steps - 1; k >= 0; k--)
    {
        _arb_vec_scalar_mul_2exp_si(d0, distances, n, k);
        for (j = 0; j < nb; j++)
        {
            _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, k);
            if (easy_steps[j] > k)
            {
                /* theta(z, tau)^2 = sum of theta(0, 2tau) theta(2z, 2tau) */
                acb_theta_ql_step_1(th_next + 3 * j * n, th_init,
                    th_init + 3 * j * n, rts + k * (3 * n * nb_steps) + j * (3 * n),
                    d0, d, g, prec);
            }
            else
            {
                /* theta(z + t, tau)^2 = sum of theta(0, 2tau) theta(2z + 2t, 2tau) */
                acb_theta_ql_step_1(th_next + 3 * j * n + n, th_init,
                    th_init + 3 * j * n + n,
                    rts + k * (3 * n * nb_steps) + j * (3 * n) + n, d0, d, g, prec);
                /* theta(z + 2t, tau)^2 = sum of theta(0, 2tau) theta(2z + 4t, 2tau) */
                acb_theta_ql_step_1(th_next + 3 * j * n + 2 * n, th_init,
                    th_init + 3 * j * n + 2 * n,
                    rts + k * (3 * n * nb_steps) + j * (3 * n) + 2 * n, d0, d, g, prec);
                if ((easy_steps[j] == k) || (j == 0))
                {
                    /* theta(z, tau) theta(z + 2t, tau)
                       = sum of theta(2t, 2tau) theta(2z + 2t, 2tau) */
                    acb_theta_agm_mul_tight(th_next + 3 * j * n, th_init + n,
                        th_init + 3 * j * n + n, d0, d, g, prec);
                    _acb_vec_scalar_mul_2exp_si(th_next + 3 * j * n,
                        th_next + 3 * j * n, n, g);
                    for (a = 0; a < n; a++)
                    {
                        acb_div(&th_next[3 * j * n + a], &th_next[3 * j * n + a],
                            &th_next[3 * j * n + 2 * n + a], prec);
                    }
                }
            }
        }
        /* Swap th_next and th_init */
        th_temp = th_init;
        th_init = th_next;
        th_next = th_temp;
    }

    for (j = 0; j < nb; j++)
    {
        _acb_vec_set(th + j * n, th_init + 3 * j * n, n);
    }

    _acb_vec_clear(th_next, 3 * n * nb);
    _acb_vec_clear(th_temp, 3 * n * nb);
    _arb_vec_clear(d, n);
    _arb_vec_clear(d0, n);
}

    /* Do we need this f factor ?
    acb_siegel_yinv(Yinv, tau, prec);
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(w, Yinv, y, prec);
    arb_dot(f, NULL, 1, y, 1, w, 1, g, prec);
    arb_const_pi(c, prec);
    arb_mul(f, f, c, prec); */

