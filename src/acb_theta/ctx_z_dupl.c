/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb.h"
#include "acb_theta.h"

void
acb_theta_ctx_z_dupl(acb_theta_ctx_z_t ctx, slong prec)
{
    slong g = ctx->g;
    acb_ptr temp;
    slong j, k;

    /* Compute exponentials, swapping vectors around if g > 1 */
    if (g == 1)
    {
        for (j = 0; j < g; j++)
        {
            acb_sqr(&acb_theta_ctx_exp_z(ctx)[k], &acb_theta_ctx_exp_z(ctx)[k], prec);
        }
    }
    else
    {
        temp = acb_theta_ctx_exp_z(ctx);
        acb_theta_ctx_exp_z(ctx) = acb_theta_ctx_exp_2z(ctx);
        acb_theta_ctx_exp_2z(ctx) = temp;
        temp = acb_theta_ctx_exp_z_inv(ctx);
        acb_theta_ctx_exp_z_inv(ctx) = acb_theta_ctx_exp_2z_inv(ctx);
        acb_theta_ctx_exp_2z_inv(ctx) = temp;
        for (j = 0; j < g; j++)
        {
            acb_sqr(&acb_theta_ctx_exp_2z(ctx)[k], &acb_theta_ctx_exp_z(ctx)[k], prec);
            acb_theta_ctx_sqr_inv(&acb_theta_ctx_exp_2z_inv(ctx)[k],
                &acb_theta_ctx_exp_z_inv(ctx)[k], &acb_theta_ctx_exp_2z(ctx)[k],
                acb_theta_ctx_is_real(ctx), prec);
        }
    }

    acb_sqr(acb_theta_ctx_c(ctx), acb_theta_ctx_c(ctx), prec);
    /* r does not change. */
    if (g > 1)
    {
        arb_t sqrt2;

        arb_init(sqrt2);
        arb_set_si(sqrt2, 2);
        arb_sqrt(sqrt2, sqrt2, prec);

        arb_sqr(acb_theta_ctx_u(ctx), acb_theta_ctx_u(ctx), prec);
        _arb_vec_scalar_mul(acb_theta_ctx_v(ctx), acb_theta_ctx_v(ctx), g, sqrt2, prec);

        arb_clear(sqrt2);
    }
}
