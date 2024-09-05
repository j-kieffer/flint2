/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "arb_mat.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"

void
acb_theta_dist_a0(arb_ptr d, acb_srcptr z, const acb_mat_t tau, slong prec)
{
    slong g = acb_mat_nrows(tau);
    slong n = 1 << g;
    arb_mat_t yinv, cho;
    arb_ptr v, w;
    ulong a;

    arb_mat_init(yinv, g, g);
    arb_mat_init(cho, g, g);
    v = _arb_vec_init(g);
    w = _arb_vec_init(g);

    acb_siegel_cho_yinv(cho, yinv, tau, prec);

    _acb_vec_get_imag(v, z, g);
    arb_mat_vector_mul_col(v, yinv, v, prec);

    for (a = 0; a < n; a++)
    {
        acb_theta_char_get_arb(w, a, g);
        _arb_vec_add(w, v, w, g, prec);
        arb_mat_vector_mul_col(w, cho, w, prec);
        acb_theta_dist_lat(&d[a], w, cho, prec);
    }

    arb_mat_clear(yinv);
    arb_mat_clear(cho);
    _arb_vec_clear(v, g);
    _arb_vec_clear(w, g);
}
