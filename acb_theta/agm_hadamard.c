
#include "acb_theta.h"

void
acb_theta_agm_hadamard(acb_ptr r, acb_srcptr a, slong g, slong prec)
{
    acb_ptr v;
    slong half;

    if (g == 0)
        acb_set(&r[0], &a[0]);
    else
    {
        half = 1 << (g - 1);
        v = _acb_vec_init(1 << g);

        acb_theta_agm_hadamard(v, a, g - 1, prec);
        acb_theta_agm_hadamard(v + half, a + half, g - 1, prec);
        _acb_vec_add(r, v, v + half, half, prec);
        _acb_vec_sub(r + half, v, v + half, half, prec);

        _acb_vec_clear(v, 1 << g);
    }
}
