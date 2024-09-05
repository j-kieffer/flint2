/*
    Copyright (C) 2023 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Include functions *********************************************************/

#include "t-00.c"
#include "t-agm_distances.c"
#include "t-agm_hadamard.c"
#include "t-agm_mul.c"
#include "t-agm_mul_tight.c"
#include "t-agm_sqrt.c"
#include "t-all.c"
#include "t-all_notransform.c"
#include "t-char_dot.c"
#include "t-char_is_even.c"
#include "t-char_is_goepel.c"
#include "t-char_is_syzygous.c"
#include "t-ctx_exp_inv.c"
#include "t-ctx_sqr_inv.c"
#include "t-ctx_tau_set.c"
#include "t-ctx_tau_dupl.c"
#include "t-ctx_z_set.c"
#include "t-ctx_z_add_real.c"
#include "t-ctx_z_copy.c"
#include "t-ctx_z_dupl.c"
#include "t-ctx_z_shift_a0.c"
#include "t-eld_border.c"
#include "t-eld_points.c"
#include "t-g2_character.c"
#include "t-g2_chi10.c"
#include "t-g2_chi12.c"
#include "t-g2_chi35.c"
#include "t-g2_chi3_6.c"
#include "t-g2_chi5.c"
#include "t-g2_covariants.c"
#include "t-g2_covariants_lead.c"
#include "t-g2_detk_symj.c"
#include "t-g2_psi4.c"
#include "t-g2_psi6.c"
#include "t-g2_sextic.c"
#include "t-g2_sextic_chi5.c"
#include "t-g2_transvectant.c"
#include "t-g2_transvectant_lead.c"
#include "t-jet_00.c"
#include "t-jet_all.c"
#include "t-jet_all_notransform.c"
#include "t-jet_one_notransform.c"
#include "t-jet_compose.c"
#include "t-jet_error_bounds.c"
#include "t-jet_mul.c"
#include "t-jet_naive_radius.c"
#include "t-jet_ql_bounds.c"
#include "t-jet_ql_radius.c"
#include "t-jet_tuples.c"
#include "t-naive_radius.c"
#include "t-naive_term.c"
#include "t-one_notransform.c"
#include "t-ql_setup.c"
#include "t-ql_steps.c"
#include "t-ql_lower_dim.c"
#include "t-ql_exact.c"
#include "t-siegel_cocycle.c"
#include "t-siegel_is_reduced.c"
#include "t-siegel_reduce.c"
#include "t-siegel_transform.c"
#include "t-sp2gz_decompose.c"
#include "t-sp2gz_inv.c"
#include "t-sp2gz_is_correct.c"
#include "t-sp2gz_set_blocks.c"
#include "t-sum_00.c"
#include "t-sum_a0_tilde.c"
#include "t-sum_all_tilde.c"
#include "t-sum_jet_00.c"
#include "t-sum_jet_all.c"
#include "t-transform_char.c"
#include "t-transform_kappa.c"
#include "t-transform_proj.c"
#include "t-transform_sqrtdet.c"

/* Array of test functions ***************************************************/

test_struct tests[] =
{
    TEST_FUNCTION(acb_theta_00),
    TEST_FUNCTION(acb_theta_agm_distances),
    /* TEST_FUNCTION(acb_theta_agm_hadamard),*/
    TEST_FUNCTION(acb_theta_agm_mul),
    /* TEST_FUNCTION(acb_theta_agm_mul_tight),
    TEST_FUNCTION(acb_theta_agm_sqrt), */
    TEST_FUNCTION(acb_theta_all),
    TEST_FUNCTION(acb_theta_all_notransform),
    TEST_FUNCTION(acb_theta_char_dot),
    /* TEST_FUNCTION(acb_theta_char_is_even),
    TEST_FUNCTION(acb_theta_char_is_goepel),
    TEST_FUNCTION(acb_theta_char_is_syzygous),*/
    TEST_FUNCTION(acb_theta_ctx_exp_inv),
    TEST_FUNCTION(acb_theta_ctx_sqr_inv),
    TEST_FUNCTION(acb_theta_ctx_tau_set),
    TEST_FUNCTION(acb_theta_ctx_tau_dupl),
    TEST_FUNCTION(acb_theta_ctx_z_set),
    TEST_FUNCTION(acb_theta_ctx_z_add_real),
    TEST_FUNCTION(acb_theta_ctx_z_copy),
    TEST_FUNCTION(acb_theta_ctx_z_dupl),
    TEST_FUNCTION(acb_theta_ctx_z_shift_a0),
    TEST_FUNCTION(acb_theta_eld_border),
    TEST_FUNCTION(acb_theta_eld_points),
    TEST_FUNCTION(acb_theta_g2_character),
    TEST_FUNCTION(acb_theta_g2_chi10),
    TEST_FUNCTION(acb_theta_g2_chi12),
    TEST_FUNCTION(acb_theta_g2_chi35),
    TEST_FUNCTION(acb_theta_g2_chi3_6),
    TEST_FUNCTION(acb_theta_g2_chi5),
    TEST_FUNCTION(acb_theta_g2_covariants),
    TEST_FUNCTION(acb_theta_g2_covariants_lead),
    TEST_FUNCTION(acb_theta_g2_detk_symj),
    TEST_FUNCTION(acb_theta_g2_psi4),
    TEST_FUNCTION(acb_theta_g2_psi6),
    TEST_FUNCTION(acb_theta_g2_sextic),
    TEST_FUNCTION(acb_theta_g2_sextic_chi5),
    TEST_FUNCTION(acb_theta_g2_transvectant),
    TEST_FUNCTION(acb_theta_g2_transvectant_lead),
    TEST_FUNCTION(acb_theta_jet_00),
    TEST_FUNCTION(acb_theta_jet_all),
    TEST_FUNCTION(acb_theta_jet_all_notransform),
    TEST_FUNCTION(acb_theta_jet_one_notransform),
    /* TEST_FUNCTION(acb_theta_jet_compose), */
    TEST_FUNCTION(acb_theta_jet_error_bounds),
    /* TEST_FUNCTION(acb_theta_jet_mul), */
    TEST_FUNCTION(acb_theta_jet_naive_radius),
    TEST_FUNCTION(acb_theta_jet_ql_bounds),
    /* TEST_FUNCTION(acb_theta_jet_ql_radius),
    TEST_FUNCTION(acb_theta_jet_tuples), */
    TEST_FUNCTION(acb_theta_naive_radius),
    /* TEST_FUNCTION(acb_theta_naive_term),*/
    TEST_FUNCTION(acb_theta_one_notransform),
    TEST_FUNCTION(acb_theta_ql_setup),
    TEST_FUNCTION(acb_theta_ql_steps),
    TEST_FUNCTION(acb_theta_ql_lower_dim),
    TEST_FUNCTION(acb_theta_ql_exact),
    /* TEST_FUNCTION(acb_theta_siegel_cocycle),
    TEST_FUNCTION(acb_theta_siegel_is_reduced),
    TEST_FUNCTION(acb_theta_siegel_reduce),
    TEST_FUNCTION(acb_theta_siegel_transform),
    TEST_FUNCTION(acb_theta_sp2gz_decompose),
    TEST_FUNCTION(acb_theta_sp2gz_inv),
    TEST_FUNCTION(acb_theta_sp2gz_is_correct),
    TEST_FUNCTION(acb_theta_sp2gz_set_blocks), */
    TEST_FUNCTION(acb_theta_sum_00),
    TEST_FUNCTION(acb_theta_sum_a0_tilde),
    TEST_FUNCTION(acb_theta_sum_all_tilde),
    TEST_FUNCTION(acb_theta_sum_jet_00),
    TEST_FUNCTION(acb_theta_sum_jet_all),
    /* TEST_FUNCTION(acb_theta_transform_char),
    TEST_FUNCTION(acb_theta_transform_kappa),
    TEST_FUNCTION(acb_theta_transform_proj),
    TEST_FUNCTION(acb_theta_transform_sqrtdet)*/
};

/* main function *************************************************************/

TEST_MAIN(tests)
