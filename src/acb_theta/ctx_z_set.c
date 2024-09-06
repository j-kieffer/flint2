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
#include "arb_mat.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_ctx_z_set(acb_theta_ctx_z_t ctx, acb_srcptr z, const acb_theta_ctx_tau_t ctx_tau, slong prec)
{
    slong g = ctx_tau->g;
    arb_t u;
    arb_ptr y, t;
    acb_ptr s;
    slong k;
    int is_real;

    arb_init(u);
    y = _arb_vec_init(g);
    t = _arb_vec_init(g);
    s = _acb_vec_init(g);

    /* Set t = Yinv y */
    _acb_vec_get_imag(y, z, g);
    arb_mat_vector_mul_col(t, &ctx_tau->yinv, y, prec);

    /* u is exp(pi y^T Yinv y) and v is cho * t */
    arb_dot(u, NULL, 0, y, 1, t, 1, g, prec);
    arb_const_pi(acb_theta_ctx_u(ctx), prec);
    arb_mul(acb_theta_ctx_u(ctx), acb_theta_ctx_u(ctx), u, prec);
    arb_exp(acb_theta_ctx_u(ctx), acb_theta_ctx_u(ctx), prec);
    arb_inv(acb_theta_ctx_uinv(ctx), acb_theta_ctx_u(ctx), prec);
    if (g > 1)
    {
        arb_mat_vector_mul_col(ctx->v, &ctx_tau->cho, t, prec);
    }

    /* Set z_is_real, exp_z, exp_z_inv, exp_2z, exp_2z_inv */
    for (k = 0; k < g; k++)
    {
        acb_exp_pi_i(&acb_theta_ctx_exp_z(ctx)[k], &z[k], prec);
        if (g > 1)
        {
            is_real = acb_is_real(&z[k]);
            acb_sqr(&acb_theta_ctx_exp_2z(ctx)[k], &acb_theta_ctx_exp_z(ctx)[k], prec);
            acb_theta_ctx_exp_inv(&acb_theta_ctx_exp_z_inv(ctx)[k],
                &acb_theta_ctx_exp_z(ctx)[k], &z[k], is_real, prec);
            acb_theta_ctx_sqr_inv(&acb_theta_ctx_exp_2z_inv(ctx)[k],
                &acb_theta_ctx_exp_z_inv(ctx)[k], &acb_theta_ctx_exp_2z(ctx)[k],
                is_real, prec);
        }
    }
    acb_theta_ctx_is_real(ctx) = _acb_vec_is_real(z, g);

    arb_clear(u);
    _arb_vec_clear(y, g);
    _arb_vec_clear(t, g);
    _acb_vec_clear(s, g);
}
