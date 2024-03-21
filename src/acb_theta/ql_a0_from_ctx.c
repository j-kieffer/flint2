/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

/* In this function, we assume that the input is exact and that the imaginary
   part is reasonable. */
int
acb_theta_ql_a0_from_ctx(acb_ptr th, const acb_theta_ql_ctx_t ctx, slong prec)
{
    slong g = acb_theta_ql_ctx_g(ctx);
    slong split, nb_steps, padding;
    int res;

    nb_steps = acb_theta_ql_a0_nb_steps_from_ctx(&split, ctx, prec);
    padding = nb_steps * (guard + g);
    res = acb_theta_ql_a0_steps_from_ctx(th, ctx, nb_steps, split, prec + padding, guard,
		&acb_theta_ql_a0_from_ctx);
    return res;
}
