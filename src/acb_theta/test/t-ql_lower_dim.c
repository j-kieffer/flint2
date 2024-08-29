/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "acb_mat.h"
#include "acb_theta.h"

TEST_FUNCTION_START(acb_theta_ql_lower_dim, state)
{
    slong iter;

    /* Test: agrees with sum_a0_tilde */
    for (iter = 0; iter < 50 * flint_test_multiplier(); iter++)
    {
        slong g = 2 + n_randint(state, 2);
        slong n = 1 << g;
        slong s = 1 + n_randint(state, g - 1);
        slong n0 = 1 << s;
        slong nba = 1 << (g - s);
        ulong a = n_randint(state, nba);
        slong prec = 100 + n_randint(state, 100);
        slong bits = n_randint(state, 4);
        acb_mat_t tau, tau0;
        acb_ptr z;
        arb_ptr d, d0;
        acb_ptr z0s, cofactors;
        arf_t err;
        acb_theta_ctx_tau_t ctx_tau, ctx_tau0;
        acb_theta_ctx_z_t ctx_z, ctx_z0;
        acb_ptr th, th0s;
        acb_t test;
        slong nb, fullprec;
        slong j, a0;
        int res;

        acb_mat_init(tau, g, g);
        z = _acb_vec_init(g);
        d = _arb_vec_init(n);
        d0 = _arb_vec_init(n0);
        acb_theta_ctx_tau_init(ctx_tau, g);
        acb_theta_ctx_tau_init(ctx_tau0, s);
        acb_theta_ctx_z_init(ctx_z, g);
        acb_theta_ctx_z_init(ctx_z0, s);
        th = _acb_vec_init(n);
        arf_init(err);
        acb_init(test);

        acb_siegel_randtest_reduced(tau, state, prec, bits);
        acb_siegel_randtest_vec(z, state, g, prec);

        acb_theta_ctx_tau_set(ctx_tau, tau, prec);
        acb_theta_ctx_z_set(ctx_z, z, ctx_tau, prec);
        acb_theta_dist_a0(d, z, tau, ACB_THETA_LOW_PREC);
        acb_theta_sum_a0_tilde(th, ctx_z, 1, ctx_tau, d, prec);

        res = acb_theta_ql_lower_dim(&z0s, &cofactors, &nb, err,
            &fullprec, z, tau, d, s, a, prec);

        th0s = _acb_vec_init(nb * n0);
        acb_mat_window_init(tau0, tau, 0, 0, s, s);

        acb_theta_ctx_tau_set(ctx_tau0, tau0, prec);
        for (j = 0; j < nb; j++)
        {
            acb_theta_dist_a0(d0, z0s + j * s, tau0, ACB_THETA_LOW_PREC);
            acb_theta_ctx_z_set(ctx_z0, z0s + j * s, ctx_tau0, prec);
            acb_theta_sum_a0_tilde(th0s + j * n0, ctx_z0, 1, ctx_tau0, d0, prec);
            _acb_vec_scalar_mul(th0s + j * n0, th0s + j * n0, n0, &cofactors[j], prec);
        }

        if (res)
        {
            for (a0 = 0; a0 < n0; a0++)
            {
                acb_zero(test);
                for (j = 0; j < nb; j++)
                {
                    acb_add(test, test, &th0s[j * n0 + a0], prec);
                }
                acb_add_error_arf(test, err);

                if (!acb_overlaps(test, &th[(a0 << (g - s)) + a]))
                {
                    flint_printf("FAIL\n");
                    flint_printf("a0 = %wd, a = %wd\n", a0, a);
                    flint_printf("th: ");
                    acb_printd(&th[(a0 << (g - s)) + a], 5);
                    flint_printf("\ntest: ");
                    acb_printd(test, 5);
                    flint_printf("\n");
                    flint_abort();
                }
            }
        }

        acb_mat_clear(tau);
        acb_mat_window_clear(tau0);
        _acb_vec_clear(z, g);
        _arb_vec_clear(d, n);
        _arb_vec_clear(d0, n0);
        _acb_vec_clear(z0s, nb * s);
        _acb_vec_clear(cofactors, nb);
        arf_clear(err);
        acb_theta_ctx_tau_clear(ctx_tau);
        acb_theta_ctx_tau_clear(ctx_tau0);
        acb_theta_ctx_z_clear(ctx_z);
        acb_theta_ctx_z_clear(ctx_z0);
        _acb_vec_clear(th, n);
        _acb_vec_clear(th0s, nb * n0);
        acb_clear(test);
    }

    TEST_FUNCTION_END(state);
}
