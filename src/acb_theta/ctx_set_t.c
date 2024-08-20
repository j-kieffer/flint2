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

void
acb_theta_ctx_set_t(acb_theta_ctx_t ctx, const acb_ptr t, slong prec)
{
    slong g = acb_theta_ctx_g(ctx);
    slong j;

    FLINT_ASSERT(_acb_vec_is_real(t, g));

    ctx->t_is_zero = _acb_vec_is_zero(t, g);
    if (ctx->t_is_zero)
    {
        return;
    }
    acb_theta_ctx_set_z(ctx, t, 1, prec);

    /* Propagate 2t using squares/conjs */
    for (j = 0; j < g; j++)
    {
        acb_sqr(&acb_theta_ctx_exp_zs(ctx)[2 * g + j],
            &acb_theta_ctx_exp_zs(ctx)[g + j], prec);
        acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[2 * g + j],
            &acb_theta_ctx_exp_zs_inv(ctx)[g + j]);
        /* Ignore cs and vs which are not used. */
    }

    /* Propagate z+t, z+2t using mults */
    if (!ctx->z_is_zero)
    {
        for (j = 0; j < g; j++)
        {
            acb_mul(&acb_theta_ctx_exp_zs(ctx)[4 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[3 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[g + j], prec);
            acb_mul(&acb_theta_ctx_exp_zs(ctx)[5 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[3 * g + j],
                &acb_theta_ctx_exp_zs(ctx)[2 * g + j], prec);
            if (ctx->z_is_real)
            {
                acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[4 * g + j],
                    &acb_theta_ctx_exp_zs(ctx)[4 * g + j]);
                acb_conj(&acb_theta_ctx_exp_zs_inv(ctx)[5 * g + j],
                    &acb_theta_ctx_exp_zs(ctx)[5 * g + j]);
            }
            else
            {
                acb_mul(&acb_theta_ctx_exp_zs_inv(ctx)[4 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[3 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[g + j], prec);
                acb_mul(&acb_theta_ctx_exp_zs_inv(ctx)[5 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[3 * g + j],
                    &acb_theta_ctx_exp_zs_inv(ctx)[2 * g + j], prec);
            }
        }
        /* Ignore cs and vs which are not used. */
    }
}
