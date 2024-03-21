/*
    Copyright (C) 2023 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

void
acb_theta_all_notransform(acb_ptr th, acb_srcptr zs, slong nb,
	const acb_mat_t tau, int sqr, slong prec);
{
	int res;

	use_ql = (...);
	if (use_ql)
	{
		res = acb_theta_ql_all_notransform(th, zs, nb, tau, sqr, prec);
	}
	if (!res || !use_ql)
	{
		acb_theta_naive_all(th, zs, nb, tau, sqr, prec);
	}
}
