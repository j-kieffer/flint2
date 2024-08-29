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

TEST_FUNCTION_START(acb_theta_ql_steps, state)
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
        slong nb_steps = n_randint(state, 6);
        acb_mat_t tau;
        acb_ptr zs, t, rts;
        acb_ptr th, th_init, test;
        arb_ptr distances;
        arb_t y;
        slong * easy_steps;
        slong guard, hp;
        slong j, k;
        int res;

        acb_mat_init(tau, g, g);
        zs = _acb_vec_init(nb * g);
        t = _acb_vec_init(g);
        rts = _acb_vec_init(nb * 3 * n * nb_steps);
        th = _acb_vec_init(n * nb);
        th_init = _acb_vec_init(3 * n * nb);
        test = _acb_vec_init(n * nb);
        distances = _arb_vec_init(nb * n);
        arb_init(y);
        easy_steps = flint_malloc(nb * sizeof(slong));

        /* Sample tau with reasonable imaginary part */
        res = 0;
        while(!res)
        {
            acb_siegel_randtest_reduced(tau, state, prec, bits);
            arb_sub_si(y, acb_imagref(acb_mat_entry(tau, g - 1, g - 1)), 200, prec);
            res = arb_is_negative(y);
        }
        acb_siegel_randtest_vec(zs + g, state, (nb - 1) * g, prec);

        /* Strip tau, z of error */
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

        /* Compute distances */
        for (j = 0; j < nb; j++)
        {
            acb_theta_dist_a0(distances + j * n, zs + j * g, tau, prec);
        }

        /* For now, all = 0 */
        res = acb_theta_ql_setup(rts, t, &guard, easy_steps, zs, nb, tau,
            distances, nb_steps, 0, prec);

        /* flint_printf("\n\ng = %wd, prec = %wd, nb = %wd, nb_steps = %wd\n", g, prec, nb, nb_steps);
           acb_mat_printd(tau, 5);
           _acb_vec_printd(zs, nb * g, 5); */
        hp = prec + guard * nb_steps;

        if (res)
        {
            /* Fill in th_init with correct values */
            /* Copied from ql_exact */
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
            acb_theta_ctx_tau_set(ctx, new_tau, hp);
            if (easy_steps[0] < nb_steps)
            {
                _acb_vec_scalar_mul_2exp_si(new_z, t, g, nb_steps);
                acb_theta_ctx_z_set(ctxt, new_z, ctx, hp);
            }

            for (j = 0; j < nb; j++)
            {
                _acb_vec_scalar_mul_2exp_si(new_z, zs + j * g, g, nb_steps);
                acb_theta_ctx_z_set(&aux[0], new_z, ctx, hp);
                _arb_vec_scalar_mul_2exp_si(d, distances + j * n, n, nb_steps);
                if (easy_steps[j] == nb_steps)
                {
                    acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 1, ctx, d, hp);
                }
                else
                {
                    acb_theta_ctx_z_add_real(&aux[1], &aux[0], ctxt, hp);
                    acb_theta_ctx_z_add_real(&aux[2], &aux[1], ctxt, hp);
                    acb_theta_sum_a0_tilde(th_init + 3 * n * j, aux, 3, ctx, d, hp);
                }
            }

            acb_theta_ql_steps(th, th_init, rts, nb, nb_steps, distances,
                easy_steps, g, hp);

            acb_theta_ctx_tau_set(ctx, tau, prec);
            for (j = 0; j < nb; j++)
            {
                acb_theta_ctx_z_set(aux, zs + j * g, ctx, prec);
                acb_theta_sum_a0_tilde(test + j * n, aux, 1, ctx, distances + j * n, prec);
            }

            /* flint_printf("After %wd steps:\n", nb_steps);
               _acb_vec_printd(th, nb * n, 5);
               flint_printf("Result of sum_a0_tilde:\n");
               _acb_vec_printd(test, nb * n, 5); */

            if (!_acb_vec_overlaps(th, test, nb * n))
            {
                flint_printf("FAIL\n");
                flint_abort();
            }

            /* Clear */
            acb_theta_ctx_tau_clear(ctx);
            acb_theta_ctx_z_vec_clear(aux, 3);
            acb_mat_clear(new_tau);
            _acb_vec_clear(new_z, g);
            _arb_vec_clear(d, n);
            if (easy_steps[0] < nb_steps)
            {
                acb_theta_ctx_z_clear(ctxt);
            }
        }

        acb_mat_clear(tau);
        _acb_vec_clear(zs, nb * g);
        _acb_vec_clear(t, g);
        _acb_vec_clear(rts, nb * 3 * n * nb_steps);
        _acb_vec_clear(th, n * nb);
        _acb_vec_clear(th_init, 3 * n * nb);
        _acb_vec_clear(test, n * nb);
        _arb_vec_clear(distances, nb * n);
        arb_clear(y);
        flint_free(easy_steps);
    }

    TEST_FUNCTION_END(state);
}
