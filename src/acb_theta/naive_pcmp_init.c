/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_naive_pcmp_init(acb_theta_naive_pcmp_t S, acb_srcptr zs,
	const acb_mat_t tau, slong nb, slong prec)
{
	S->nbz = nb;
}
