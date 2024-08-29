/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ql_exact, state)
{
    slong iter;

    /* Test: coincides with sum_a0_tilde */
    for (iter = 0; iter < 25 * flint_test_multiplier(); iter++)
    {
        slong g = 1 + n_randint(state, 3);
        slong n = 1 << g;
        slong prec = 100 + n_randint(state, 200);
        slong bits = n_randint(state, 4);
        slong nb = 1 + n_randint(state, 4);
        int all = 0;
        int shifted_prec = 1;
        acb_mat_t tau;
        acb_ptr zs, th, test;
        acb_theta_ctx_tau_t ctx_tau;
        acb_theta_ctx_z_struct * vec;
        arb_ptr distances;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        th = _acb_vec_init(n * nb);
        test = _acb_vec_init(n * nb);
        acb_theta_ctx_tau_init(ctx_tau, g);
        vec = acb_theta_ctx_z_vec_init(nb, g);
        distances = _arb_vec_init(nb * n);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec(zs + g, state, (nb - 1) * g, prec);

        /* Strip tau, z of errors */
        for (j = 0; j < g; j++)
        {
            for (k = j; k < g; k++)
            {
                acb_get_mid(acb_mat_entry(tau, j, k), acb_mat_entry(tau, j, k));
                acb_set(acb_mat_entry(tau, k, j), acb_mat_entry(tau, j, k));
            }
        }
        for (j = 0; j < nb * g; j++)
        {
            acb_get_mid(&zs[j], &zs[j]);
        }

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_ctx_z_set(&vec[j], zs + j * g, ctx_tau, prec);
            if (shifted_prec)
            {
                acb_theta_dist_a0(distances + j * n, zs + j * g, tau, prec);
            }
        }
        acb_theta_sum_a0_tilde(test, vec, nb, ctx_tau, distances, prec);

        flint_printf("\n\ng = %wd, prec = %wd, nb = %wd, all = %wd, shifted_prec = %wd\n",
            g, prec, nb, all, shifted_prec);
        acb_mat_printd(tau, 5);
        _acb_vec_printd(zs, nb * g, 5);
        flint_printf("sum_a0_tilde:\n");
        _acb_vec_printd(test, n * nb, 5);

        /* For now, all = 0 */
        res = acb_theta_ql_exact(th, zs, nb, tau, all, shifted_prec, prec);

        flint_printf("result of ql_exact: %wd, got theta:\n", res);
        _acb_vec_printd(th, n * nb, 5);

        if (res && !_acb_vec_overlaps(th, test, nb * n))
        {
            flint_printf("FAIL\n");
            flint_abort();
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(th, n * nb);
        _acb_vec_clear(test, n * nb);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_z_vec_clear(vec, nb);
        _arb_vec_clear(distances, nb * n);
    }

    TEST_FUNCTION_END(state);
}
