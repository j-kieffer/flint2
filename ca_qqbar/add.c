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
ca_qqbar_add(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y)
{
    if (ca_qqbar_is_zero(x))
    {
        ca_qqbar_set(res, y);
    }
    else if (ca_qqbar_is_zero(y))
    {
        ca_qqbar_set(res, x);
    }
    else
    {
        ca_qqbar_binary_op(res, x, y, 0);
    }
}

