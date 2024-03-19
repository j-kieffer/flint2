/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb_theta.h"

slong
acb_theta_ql_ctx_nbz(const acb_theta_ql_ctx_t ctx)
{
	slong nbt = (ctx->t_is_zero ? 1 : 3);
	slong nbz = (ctx->z_is_zero ? 1 : 2);
	return nbt * nbz;
}
