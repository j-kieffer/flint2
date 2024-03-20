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
acb_theta_naive_exp_z(acb_ptr exp_z, acb_scrptr z, slong g, slong prec)
{
	slong j;

	for (j = 0; j < g; j++)
	{
		acb_set_round(&exp_z[j], &z[j], prec);
		acb_mul_2exp_si(&exp_z[j], &exp_z[j], 1);
		acb_exp_pi_i(&exp_z[j], &exp_z[j], prec);
	}
}
