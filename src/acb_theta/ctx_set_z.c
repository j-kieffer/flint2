/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

static void
acb_theta_naive_round(arb_ptr a, arb_srcptr v, slong g)
{
    slong j;
    fmpz_t m;

    fmpz_init(m);

    for (j = 0; j < g; j++)
    {
        if (arb_is_finite(&v[j])
            && arf_cmpabs_2exp_si(arb_midref(&v[j]), 1000000) <= 0)
        {
            arf_get_fmpz(m, arb_midref(&v[j]), ARF_RND_NEAR);
            arb_set_fmpz(&a[j], m);
        }
        else
        {
            arb_zero(&a[j]);
        }
    }

    fmpz_clear(m);
}

void
acb_theta_ctx_set_z(acb_theta_ctx_t ctx, arb_srcptr a, acb_t c,
	acb_srcptr z, slong j, slong prec)
{
	slong g = acb_theta_ctx_g(ctx);
	acb_t x;
	arb_t u;
	arb_ptr y;
	slong k;
	int res;

	FLINT_ASSERT(j >= 0 && j < acb_theta_ctx_nb(ctx));

	acb_init(x);
	arb_init(u);
	y = _arb_vec_init(g);

	/* Set v, a, c */
	/* todo: round... */

	/* Set exp_z, exp_z_inv */
	for (k = 0; k < g; k++)
	{
		acb_set_round(x, &z[k], prec);
		acb_mul_2exp_si(x, x, 1);
		acb_exp_pi_i(&acb_theta_ctx_exp_zs(ctx)[j * g + k], x, prec);
		if (acb_is_real(&z[k]))
		{
			acb_conj(x, x);
		}
		else
		{
			acb_inv(x, x, prec);
		}
		acb_set(&acb_theta_ctx_exp_zs_inv(ctx)[j * g + k], x);
	}

	/* Set c, v */
	if (g == 1)
	{
		arb_sqr(u, acb_imagref(z), prec);
		arb_neg(u, u);
		arb_const_pi(&acb_theta_ctx_exp_cs(ctx)[j], prec);
		arb_mul(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], u, prec);
		arb_exp(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], prec);
	}
	else
	{
		_acb_vec_get_imag(y, z, g);
		_arb_vec_neg(acb_theta_ctx_vs(ctx) + j * g, y, g);
		arb_mat_vector_mul_col(acb_theta_ctx_vs(ctx) + j * g, acb_theta_ctx_yinv(ctx),
			acb_theta_ctx_vs(ctx) + j * g, prec);
		arb_dot(u, y, 1, acb_theta_ctx_vs(ctx) + j * g, 1, g, prec);
		arb_const_pi(&acb_theta_ctx_exp_cs(ctx)[j], prec);
		arb_mul(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], u, prec);
		arb_exp(&acb_theta_ctx_exp_cs(ctx)[j], &acb_theta_ctx_exp_cs(ctx)[j], prec);
	}

	acb_clear(x);
	arb_clear(u);
	_arb_vec_clear(y, g);
}
