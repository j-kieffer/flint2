/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_modular.h"
#include "acb_theta.h"

static void
acb_theta_ctx_copy_tau(acb_theta_ctx_t new_ctx, const acb_theta_ctx_t ctx)
{
    slong g = acb_theta_ctx_g(ctx);
    FLINT_ASSERT(acb_theta_ctx_g(new_ctx) == g);

    acb_mat_set(acb_theta_ctx_tau(new_ctx), acb_theta_ctx_tau(ctx));
    arb_mat_set(acb_theta_ctx_y(new_ctx), acb_theta_ctx_y(ctx));
    arb_mat_set(acb_theta_ctx_yinv(new_ctx), acb_theta_ctx_yinv(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau_div_4(new_ctx), acb_theta_ctx_exp_tau_div_4(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau_div_2(new_ctx), acb_theta_ctx_exp_tau_div_2(ctx));
    acb_mat_set(acb_theta_ctx_exp_tau(new_ctx), acb_theta_ctx_exp_tau(ctx));
    if (g > 1)
    {
        arb_mat_set(acb_theta_ctx_cho(new_ctx), acb_theta_ctx_cho(ctx));
        arb_mat_set(acb_theta_ctx_choinv(new_ctx), acb_theta_ctx_choinv(ctx));
        acb_mat_set(acb_theta_ctx_exp_tau_inv(new_ctx), acb_theta_ctx_exp_tau(new_ctx));
    }
}

static void
acb_theta_ctx_shift_z(acb_theta_ctx_t new_ctx, const acb_theta_ctx_t ctx,
    slong z_start, slong z_end, ulong a, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong nb = z_end - z_start + 1;
    acb_t c, cinv, csqr, csqrinv;
    arb_t abs;
    slong j, k;

    acb_init(c);
    acb_init(cinv);
    acb_init(csqr);
    acb_init(csqrinv);
    arb_init(abs);

    /* Replace exp_zs by analogs for z + tau a/2 */
    _acb_vec_set(acb_theta_ctx_exp_zs(new_ctx), acb_theta_ctx_exp_zs(ctx) + z_start * g, nb * g);
    _acb_vec_set(acb_theta_ctx_exp_zs_inv(new_ctx), acb_theta_ctx_exp_zs_inv(ctx) + z_start * g, nb * g);
    _acb_vec_set(acb_theta_ctx_exp_2zs(new_ctx), acb_theta_ctx_exp_2zs(ctx) + z_start * g, nb * g);
    _acb_vec_set(acb_theta_ctx_exp_2zs_inv(new_ctx), acb_theta_ctx_exp_2zs_inv(ctx) + z_start * g, nb * g);
    for (j = 0; j < g; j++)
    {
        if (!((a >> j) & 1))
        {
            continue; /* do nothing if aj = 0 */
        }
        acb_one(c);
        for (k = 0; k < g; k++)
        {
            if (!((a >> k) & 1))
            {
                continue;
            }
            if (k < j)
            {
                acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), k, j), prec);
            }
            else if (k == j)
            {
                acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_2(ctx), k, k), prec);
            }
            else
            {
                acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
            }
        }
        acb_inv(cinv, c, prec);
        acb_sqr(csqr, c, prec);
        acb_sqr(csqrinv, cinv, prec);
        for (k = 0; k < nb; k++)
        {
            acb_mul(&acb_theta_ctx_exp_zs(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_zs(new_ctx)[k * g + j], c, prec);
            acb_mul(&acb_theta_ctx_exp_zs_inv(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_zs_inv(new_ctx)[k * g + j], cinv, prec);
            acb_mul(&acb_theta_ctx_exp_2zs(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_2zs(new_ctx)[k * g + j], csqr, prec);
            acb_mul(&acb_theta_ctx_exp_2zs_inv(new_ctx)[k * g + j],
                &acb_theta_ctx_exp_2zs_inv(new_ctx)[k * g + j], csqrinv, prec);
        }
    }

    /* For each z, compute cofactor exp(pi i a^T z) */
    for (j = 0; j < nb; j++)
    {
        acb_one(c);
        for (k = 0; k < g; k++)
        {
            if (!((a >> k) & 1))
            {
                continue;
            }
            acb_mul(c, c, &acb_theta_ctx_exp_zs(ctx)[(z_start + j) * g + k], prec);
        }
        acb_set(&acb_theta_ctx_cs(new_ctx)[j], c);
    }

    /* Compute common cofactor exp(pi i/4 a^T tau a) */
    acb_one(c);
    for (j = 0; j < g; j++)
    {
        if (!((a >> j) & 1))
        {
            continue; /* do nothing if aj = 0 */
        }
        acb_mul(c, c, acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), j, k), prec);
    }
    _acb_vec_scalar_mul(acb_theta_ctx_cs(new_ctx), acb_theta_ctx_cs(new_ctx), nb, c, prec);

    /* Compute cs, us, as, vs */
    _arb_vec_set(acb_theta_ctx_as(new_ctx), acb_theta_ctx_as(ctx) + z_start * g, nb * g);
    for (j = 0; j < nb; j++)
    {
        acb_abs(abs, &acb_theta_ctx_cs(new_ctx)[j], prec);
        arb_mul(&acb_theta_ctx_us(new_ctx)[j], &acb_theta_ctx_us(ctx)[j], abs, prec);
        acb_mul(&acb_theta_ctx_cs(new_ctx)[j], &acb_theta_ctx_cs(new_ctx)[j], &acb_theta_ctx_cs(ctx)[j], prec);
    }
    if (g > 1)
    {
        arb_ptr v_shift;
        v_shift = _arb_vec_init(g);

        acb_theta_char_get_arb(v_shift, a, g);
        arb_mat_vector_mul_col(v_shift, acb_theta_ctx_cho(ctx), v_shift, prec);
        for (j = 0; j < nb; j++)
        {
            _arb_vec_add(acb_theta_ctx_vs(new_ctx) + j * g, v_shift,
                acb_theta_ctx_vs(new_ctx) + (z_start + j) * g, g, prec);
        }

        _arb_vec_clear(v_shift, g);
    }

    acb_clear(c);
    acb_clear(cinv);
    acb_clear(csqr);
    acb_clear(csqrinv);
    arb_clear(abs);
}

void
acb_theta_sum_a0(acb_ptr th, const acb_theta_ctx_t ctx, slong z_start,
    slong z_end, int z_is_real, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong n = 1 << g;
    slong nb = z_start - z_end + 1;
    acb_ptr res;
    acb_theta_ctx_t new_ctx;
    slong new_prec;
    slong a, j;

    if (g == 1)
    {
        res = _acb_vec_init(4);
        for (j = 0; j < nb; j++)
        {
            /* acb_modular_theta_sum takes shifted precisions into account */
            acb_modular_theta_sum(&res[0], &res[1], &res[2], &res[3],
                &acb_theta_ctx_exp_zs(ctx)[z_start + j], 0,
                acb_mat_entry(acb_theta_ctx_exp_tau(ctx), 0, 0), 1, prec);
            acb_mul(&th[2 * j], &res[2], &acb_theta_ctx_cs(ctx)[z_start + j], prec);
            acb_mul(&th[2 * j + 1], &res[1], &acb_theta_ctx_cs(ctx)[z_start + j], prec);
            acb_mul(&th[2 * j + 1], &th[2 * j + 1],
                acb_mat_entry(acb_theta_ctx_exp_tau_div_4(ctx), 0, 0), prec);
        }
        _acb_vec_clear(res, 4);
    }
    else
    {
        /* Update the context for each a, call sum_00 with the right precision */
        acb_theta_ctx_init(new_ctx, g, nb);
        res = _acb_vec_init(nb);

        acb_theta_ctx_copy_tau(new_ctx, ctx);
        for (a = 0; a < n; a++)
        {
            acb_theta_ctx_shift_z(new_ctx, ctx, z_start, z_end, a, prec);
            new_prec = prec + acb_theta_dist_addprec(
                (z_is_real ? &acb_theta_ctx_d0(ctx)[a] : &acb_theta_ctx_d(ctx)[a]));
            acb_theta_sum_00(res, new_ctx, new_prec);
            for (j = 0; j < nb; j++)
            {
                acb_set(&th[n * j + a], &res[j]);
            }
        }

        _acb_vec_clear(res, nb);
        acb_theta_ctx_clear(new_ctx);
    }
}
