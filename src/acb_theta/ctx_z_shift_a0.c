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
#include "acb_theta.h"

static int
acb_theta_aj_is_zero(ulong a, slong j, slong g)
{
    return !((a >> (g - 1 - j)) & 1);
}

void
acb_theta_ctx_z_shift_a0(acb_theta_ctx_z_t res, acb_t c, const acb_theta_ctx_z_t ctx,
    const acb_theta_ctx_tau_t ctx_tau, ulong a, slong prec)
{
    slong g = acb_theta_ctx_g(ctx_tau);
    arb_ptr v_shift;
    arb_t abs;
    slong j;

    FLINT_ASSERT(g > 1);
    v_shift = _arb_vec_init(g);
    arb_init(abs);

    /* Replace exp_z by analogs for z + tau a/2 */
    for (j = 0; j < g; j++)
    {
        acb_mul(&acb_theta_ctx_exp_z(res)[j], &acb_theta_ctx_exp_z(ctx)[j],
            &acb_theta_ctx_exp_tau_a_div_2(ctx_tau, a)[j], prec);
        acb_mul(&acb_theta_ctx_exp_2z(res)[j], &acb_theta_ctx_exp_2z(ctx)[j],
            &acb_theta_ctx_exp_tau_a(ctx_tau, a)[j], prec);
        acb_mul(&acb_theta_ctx_exp_z_inv(res)[j], &acb_theta_ctx_exp_z_inv(ctx)[j],
            &acb_theta_ctx_exp_tau_a_div_2_inv(ctx_tau, a)[j], prec);
        acb_mul(&acb_theta_ctx_exp_2z_inv(res)[j], &acb_theta_ctx_exp_2z_inv(ctx)[j],
            &acb_theta_ctx_exp_tau_a_inv(ctx_tau, a)[j], prec);
    }

    /* Compute cofactor exp(pi i a^T z), and multiply by common cofactor
       exp(pi i/4 a^T tau a) */
    acb_one(c);
    for (j = 0; j < g; j++)
    {
        if (acb_theta_aj_is_zero(a, j, g))
        {
            continue;
        }
        acb_mul(c, c, &acb_theta_ctx_exp_z(ctx)[j], prec);
    }
    acb_mul(c, c, acb_theta_ctx_exp_a_tau_a_div_4(ctx_tau, a), prec);

    /* Compute v; u and uinv must be multiplied by abs(c) */
    acb_abs(abs, c, prec);
    arb_div(acb_theta_ctx_u(res), acb_theta_ctx_u(ctx), abs, prec);
    arb_mul(acb_theta_ctx_uinv(res), acb_theta_ctx_uinv(ctx), abs, prec);

    acb_theta_char_get_arb(v_shift, a, g);
    arb_mat_vector_mul_col(v_shift, acb_theta_ctx_cho(ctx_tau), v_shift, prec);
    _arb_vec_add(acb_theta_ctx_v(res), v_shift, acb_theta_ctx_v(ctx), g, prec);

    _arb_vec_clear(v_shift, g);
    arb_clear(abs);
}
