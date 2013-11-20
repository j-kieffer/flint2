/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

******************************************************************************/


#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_divrem_basecase)(TEMPLATE(T, struct) * Q, TEMPLATE(T, struct) * R,
                         const TEMPLATE(T, struct) * A, long lenA, const TEMPLATE(T, struct) * B,
                         long lenB, const TEMPLATE(T, t) invB, const TEMPLATE(T, ctx_t) ctx)
{
    long iQ, iR;

    if (R != A)
        _TEMPLATE(T, poly_set)(R, A, lenA, ctx);

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (TEMPLATE(T, is_zero)(R + iR, ctx))
            TEMPLATE(T, zero)(Q + iQ, ctx);
        else
        {
            TEMPLATE(T, mul)(Q + iQ, R + iR, invB, ctx);

            _TEMPLATE(T, TEMPLATE(poly_scalar_submul, T))(R + iQ, B, lenB, Q + iQ, ctx);
        }
    }
}

void
TEMPLATE(T, poly_divrem_basecase)(TEMPLATE(T, poly_t) Q, TEMPLATE(T, poly_t) R,
                        const TEMPLATE(T, poly_t) A, const TEMPLATE(T, poly_t) B,
                        const TEMPLATE(T, ctx_t) ctx)
{
    const long lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1;
    TEMPLATE(T, struct) *q, *r;
    TEMPLATE(T, t) invB;

    if (lenA < lenB)
    {
        TEMPLATE(T, poly_set)(R, A, ctx);
        TEMPLATE(T, poly_zero)(Q, ctx);
        return;
    }

    TEMPLATE(T, init)(invB, ctx);
    TEMPLATE(T, inv)(invB, TEMPLATE(T, poly_lead)(B, ctx), ctx);

    if (Q == A || Q == B)
    {
        q = _TEMPLATE(T, vec_init)(lenQ, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length)(Q, lenQ, ctx);
        q = Q->coeffs;
    }
    if (R == B)
    {
        r = _TEMPLATE(T, vec_init)(lenA, ctx);
    }
    else
    {
        TEMPLATE(T, poly_fit_length)(R, lenA, ctx);
        r = R->coeffs;
    }

    _TEMPLATE(T, poly_divrem_basecase)(q, r, A->coeffs, lenA,
                             B->coeffs, lenB, invB, ctx);

    if (Q == A || Q == B)
    {
        _TEMPLATE(T, vec_clear)(Q->coeffs, Q->alloc, ctx);
        Q->coeffs = q;
        Q->alloc = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _TEMPLATE(T, poly_set_length)(Q, lenQ, ctx);
    }
    if (R == B)
    {
        _TEMPLATE(T, vec_clear)(R->coeffs, R->alloc, ctx);
        R->coeffs = r;
        R->alloc = lenA;
        R->length = lenA;
    }
    _TEMPLATE(T, poly_set_length)(R, lenB - 1, ctx);
    _TEMPLATE(T, poly_normalise)(R, ctx);

    TEMPLATE(T, clear)(invB, ctx);
}


#endif
