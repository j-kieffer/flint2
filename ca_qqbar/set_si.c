/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "ca_qqbar.h"

void
ca_qqbar_set_si(ca_qqbar_t res, slong x)
{
    fmpz_t t;
    fmpz_init_set_si(t, x);
    ca_qqbar_set_fmpz(res, t);
    fmpz_clear(t);
}

