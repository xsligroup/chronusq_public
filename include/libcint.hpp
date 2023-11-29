/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#pragma once


// Every declaration in libcint should have the C ABI
extern "C" {

#include <cint.h>

// TODO: Why are we using the old interface instead of the new one?
//       Ask Xiaosong

#define CQ_CINT_INT1E_WRAPPER(name) \
  size_t c##name(double*, int*, int*, int, int*, int, double*);

#define CQ_CINT_INT2E_WRAPPER(name) \
  size_t c##name(double*, const FINT*, const FINT*, const FINT, const FINT*, const FINT, const double*, const CINTOpt*);

#define CQ_INT_WRAPPER(name) \
  size_t name(double*, FINT*, FINT*, FINT*, FINT, FINT*, FINT, double*, CINTOpt*, double*);

#define CQ_CINT_OPT_WRAPPER(name) \
  void c##name(CINTOpt**, const int*, const int, const int*, const int, const double*);

//
//  Old libcint declarations here
//

/* Plain ERI (ij|kl) */
CQ_INT_WRAPPER(int2e_cart);
CQ_INT_WRAPPER(int2e_sph);

/* <i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_ovlp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ovlp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ovlp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ovlp_spinor);
CQ_INT_WRAPPER(int1e_ovlp_cart);
CQ_INT_WRAPPER(int1e_ovlp_sph);

/* <i|NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_nuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_nuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_nuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_nuc_spinor);
CQ_INT_WRAPPER(int1e_nuc_cart);
CQ_INT_WRAPPER(int1e_nuc_sph);

/* <i|OVLP |P DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_kin_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_kin_cart);
CQ_CINT_INT1E_WRAPPER(int1e_kin_sph);
CQ_CINT_INT1E_WRAPPER(int1e_kin_spinor);
CQ_INT_WRAPPER(int1e_kin_cart);
CQ_INT_WRAPPER(int1e_kin_sph);

/* <i|NABLA-RINV |CROSS P j> */
CQ_CINT_OPT_WRAPPER(int1e_ia01p_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ia01p_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ia01p_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ia01p_spinor);

/* <i|OVLP |R CROSS P j> */
CQ_CINT_OPT_WRAPPER(int1e_giao_irjxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_giao_irjxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_giao_irjxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_giao_irjxp_spinor);

/* <i|OVLP |RC CROSS P j> */
CQ_CINT_OPT_WRAPPER(int1e_cg_irxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_cg_irxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_cg_irxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_cg_irxp_spinor);

/* <i|NABLA-RINV |R j> */
CQ_CINT_OPT_WRAPPER(int1e_giao_a11part_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_giao_a11part_cart);
CQ_CINT_INT1E_WRAPPER(int1e_giao_a11part_sph);
CQ_CINT_INT1E_WRAPPER(int1e_giao_a11part_spinor);

/* <i|NABLA-RINV |RC j> */
CQ_CINT_OPT_WRAPPER(int1e_cg_a11part_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_cg_a11part_cart);
CQ_CINT_INT1E_WRAPPER(int1e_cg_a11part_sph);
CQ_CINT_INT1E_WRAPPER(int1e_cg_a11part_spinor);

/* <G i|NABLA-RINV CROSS P |j> */
CQ_CINT_OPT_WRAPPER(int1e_a01gp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_a01gp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_a01gp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_a01gp_spinor);

/* <G i|OVLP |P DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_igkin_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_igkin_cart);
CQ_CINT_INT1E_WRAPPER(int1e_igkin_sph);
CQ_CINT_INT1E_WRAPPER(int1e_igkin_spinor);

/* <G i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_igovlp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_igovlp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_igovlp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_igovlp_spinor);

/* <G i|NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_ignuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ignuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ignuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ignuc_spinor);

/* <P* i|NUC DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_pnucp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_pnucp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_pnucp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_pnucp_spinor);
CQ_INT_WRAPPER(int1e_pnucp_cart);
CQ_INT_WRAPPER(int1e_pnucp_sph);

/* <i|ZC |j> */
CQ_CINT_OPT_WRAPPER(int1e_z_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_z_cart);
CQ_CINT_INT1E_WRAPPER(int1e_z_sph);
CQ_CINT_INT1E_WRAPPER(int1e_z_spinor);

/* <i|ZC ZC |j> */
CQ_CINT_OPT_WRAPPER(int1e_zz_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_zz_cart);
CQ_CINT_INT1E_WRAPPER(int1e_zz_sph);
CQ_CINT_INT1E_WRAPPER(int1e_zz_spinor);

/* <i|RC |j> */
CQ_CINT_OPT_WRAPPER(int1e_r_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_r_cart);
CQ_CINT_INT1E_WRAPPER(int1e_r_sph);
CQ_CINT_INT1E_WRAPPER(int1e_r_spinor);
CQ_INT_WRAPPER(int1e_r_sph);

/* <i|RC DOT RC |j> */
CQ_CINT_OPT_WRAPPER(int1e_r2_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_r2_cart);
CQ_CINT_INT1E_WRAPPER(int1e_r2_sph);
CQ_CINT_INT1E_WRAPPER(int1e_r2_spinor);

/* <i|RC RC |j> */
CQ_CINT_OPT_WRAPPER(int1e_rr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_rr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_rr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_rr_spinor);

/* <i|RC RC RC |j> */
CQ_CINT_OPT_WRAPPER(int1e_rrr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_rrr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_rrr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_rrr_spinor);

/* <i|RC RC RC RC |j> */
CQ_CINT_OPT_WRAPPER(int1e_rrrr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_rrrr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_rrrr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_rrrr_spinor);

/* <i|Z |j> */
CQ_CINT_OPT_WRAPPER(int1e_z_origj_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_z_origj_cart);
CQ_CINT_INT1E_WRAPPER(int1e_z_origj_sph);
CQ_CINT_INT1E_WRAPPER(int1e_z_origj_spinor);

/* <i|Z Z |j> */
CQ_CINT_OPT_WRAPPER(int1e_zz_origj_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_zz_origj_cart);
CQ_CINT_INT1E_WRAPPER(int1e_zz_origj_sph);
CQ_CINT_INT1E_WRAPPER(int1e_zz_origj_spinor);

/* <i|R |j> */
CQ_CINT_OPT_WRAPPER(int1e_r_origj_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_r_origj_cart);
CQ_CINT_INT1E_WRAPPER(int1e_r_origj_sph);
CQ_CINT_INT1E_WRAPPER(int1e_r_origj_spinor);

/* <i|R R |j> */
CQ_CINT_OPT_WRAPPER(int1e_rr_origj_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_rr_origj_cart);
CQ_CINT_INT1E_WRAPPER(int1e_rr_origj_sph);
CQ_CINT_INT1E_WRAPPER(int1e_rr_origj_spinor);

/* <i|R DOT R |j> */
CQ_CINT_OPT_WRAPPER(int1e_r2_origj_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_r2_origj_cart);
CQ_CINT_INT1E_WRAPPER(int1e_r2_origj_sph);
CQ_CINT_INT1E_WRAPPER(int1e_r2_origj_spinor);

/* <i|OVLP |R DOT R R DOT R j> */
CQ_CINT_OPT_WRAPPER(int1e_r4_origj_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_r4_origj_cart);
CQ_CINT_INT1E_WRAPPER(int1e_r4_origj_sph);
CQ_CINT_INT1E_WRAPPER(int1e_r4_origj_spinor);

/* <P DOT P i|OVLP |P DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_p4_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_p4_cart);
CQ_CINT_INT1E_WRAPPER(int1e_p4_sph);
CQ_CINT_INT1E_WRAPPER(int1e_p4_spinor);

/* <P* i|RINV DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_prinvp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_prinvp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_prinvp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_prinvp_spinor);

/* <P* i|RINV CROSS P |j> */
CQ_CINT_OPT_WRAPPER(int1e_prinvxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_prinvxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_prinvxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_prinvxp_spinor);

/* <P* i|NUC CROSS P |j> */
CQ_CINT_OPT_WRAPPER(int1e_pnucxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_pnucxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_pnucxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_pnucxp_spinor);
CQ_INT_WRAPPER(int1e_pnucxp_cart);
CQ_INT_WRAPPER(int1e_pnucxp_sph);

/* <i|RC NABLA |j> */
CQ_CINT_OPT_WRAPPER(int1e_irp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_irp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_irp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_irp_spinor);

/* <i|RC RC NABLA |j> */
CQ_CINT_OPT_WRAPPER(int1e_irrp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_irrp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_irrp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_irrp_spinor);

/* <i|RC NABLA RC |j> */
CQ_CINT_OPT_WRAPPER(int1e_irpr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_irpr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_irpr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_irpr_spinor);

/* <i|G G |j> */
CQ_CINT_OPT_WRAPPER(int1e_ggovlp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ggovlp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ggovlp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ggovlp_spinor);

/* <i|G G NABLA DOT NABLA |j> */
CQ_CINT_OPT_WRAPPER(int1e_ggkin_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ggkin_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ggkin_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ggkin_spinor);

/* <i|G G NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_ggnuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ggnuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ggnuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ggnuc_spinor);

/* <i|G R CROSS P |j> */
CQ_CINT_OPT_WRAPPER(int1e_grjxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_grjxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_grjxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_grjxp_spinor);

/* <i|RINV |j> */
CQ_CINT_OPT_WRAPPER(int1e_rinv_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_rinv_cart);
CQ_CINT_INT1E_WRAPPER(int1e_rinv_sph);
CQ_CINT_INT1E_WRAPPER(int1e_rinv_spinor);

/* <i|NABLA-RINV |j> */
CQ_CINT_OPT_WRAPPER(int1e_drinv_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_drinv_cart);
CQ_CINT_INT1E_WRAPPER(int1e_drinv_sph);
CQ_CINT_INT1E_WRAPPER(int1e_drinv_spinor);

/* (G i j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ig1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ig1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ig1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ig1_spinor);

/* (G G i j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_gg1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gg1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gg1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gg1_spinor);

/* (G i j|R12 |G k l) */
CQ_CINT_OPT_WRAPPER(int2e_g1g2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_g1g2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_g1g2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_g1g2_spinor);

/* (P* i CROSS P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_p1vxp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_p1vxp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_p1vxp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_p1vxp1_spinor);

/* (i RC j|NABLA-R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ip1v_rc1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip1v_rc1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip1v_rc1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip1v_rc1_spinor);

/* (i R j|NABLA-R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ip1v_r1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip1v_r1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip1v_r1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip1v_r1_spinor);

/* (G i j|NABLA-R12 CROSS P |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipvg1_xp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipvg1_xp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipvg1_xp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipvg1_xp1_spinor);

/* (i j|NABLA-R12 CROSS P |G k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipvg2_xp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipvg2_xp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipvg2_xp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipvg2_xp1_spinor);

/* <i|NUC |RC CROSS P j> */
CQ_CINT_OPT_WRAPPER(int1e_inuc_rcxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_inuc_rcxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_inuc_rcxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_inuc_rcxp_spinor);

/* <i|NUC |R CROSS P j> */
CQ_CINT_OPT_WRAPPER(int1e_inuc_rxp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_inuc_rxp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_inuc_rxp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_inuc_rxp_spinor);

/* <i|OVLP |SIGMA j> */
CQ_CINT_OPT_WRAPPER(int1e_sigma_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_sigma_cart);
CQ_CINT_INT1E_WRAPPER(int1e_sigma_sph);
CQ_CINT_INT1E_WRAPPER(int1e_sigma_spinor);

/* <SIGMA DOT P i|OVLP |SIGMA SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_spsigmasp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spsigmasp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spsigmasp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spsigmasp_spinor);

/* <SIGMA DOT R i|OVLP |SIGMA DOT R j> */
CQ_CINT_OPT_WRAPPER(int1e_srsr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_srsr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_srsr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_srsr_spinor);

/* <SIGMA DOT R i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_sr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_sr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_sr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_sr_spinor);

/* <SIGMA DOT R i|OVLP |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_srsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_srsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_srsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_srsp_spinor);

/* <SIGMA DOT P i|OVLP |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_spsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spsp_spinor);

/* <SIGMA DOT P i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_sp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_sp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_sp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_sp_spinor);

/* <SIGMA DOT P i|NUC |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_spnucsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spnucsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spnucsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spnucsp_spinor);

/* <SIGMA DOT P i|RINV |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_sprinvsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_sprinvsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_sprinvsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_sprinvsp_spinor);

/* <SIGMA DOT R i|NUC |SIGMA DOT R j> */
CQ_CINT_OPT_WRAPPER(int1e_srnucsr_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_srnucsr_cart);
CQ_CINT_INT1E_WRAPPER(int1e_srnucsr_sph);
CQ_CINT_INT1E_WRAPPER(int1e_srnucsr_spinor);

/* <SIGMA DOT P i|RC |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_sprsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_sprsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_sprsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_sprsp_spinor);
CQ_INT_WRAPPER(int1e_sprsp_sph);

/* <G i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_govlp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_govlp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_govlp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_govlp_spinor);

/* <G i|NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_gnuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_gnuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_gnuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_gnuc_spinor);

/* <SIGMA CROSS RC i|SIGMA CROSS NABLA-RINV |j> */
CQ_CINT_OPT_WRAPPER(int1e_cg_sa10sa01_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10sa01_cart);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10sa01_sph);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10sa01_spinor);

/* <RC CROSS SIGMA i|OVLP |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_cg_sa10sp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10sp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10sp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10sp_spinor);

/* <RC CROSS SIGMA i|NUC |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_cg_sa10nucsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10nucsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10nucsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_cg_sa10nucsp_spinor);

/* <SIGMA CROSS R i|SIGMA CROSS NABLA-RINV |j> */
CQ_CINT_OPT_WRAPPER(int1e_giao_sa10sa01_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10sa01_cart);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10sa01_sph);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10sa01_spinor);

/* <R CROSS SIGMA i|OVLP |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_giao_sa10sp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10sp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10sp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10sp_spinor);

/* <R CROSS SIGMA i|NUC |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_giao_sa10nucsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10nucsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10nucsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_giao_sa10nucsp_spinor);

/* <i|NABLA-RINV CROSS SIGMA |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_sa01sp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_sa01sp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_sa01sp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_sa01sp_spinor);

/* <G SIGMA DOT P i|OVLP |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_spgsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spgsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spgsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spgsp_spinor);

/* <G SIGMA DOT P i|NUC |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_spgnucsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spgnucsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spgnucsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spgnucsp_spinor);

/* <G SIGMA DOT P i|NABLA-RINV CROSS SIGMA |j> */
CQ_CINT_OPT_WRAPPER(int1e_spgsa01_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spgsa01_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spgsa01_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spgsa01_spinor);

/* (SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_spsp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spsp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spsp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spsp1_spinor);

/* (SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_spsp1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spsp1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spsp1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spsp1spsp2_spinor);

/* (SIGMA DOT R i SIGMA DOT R j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_srsr1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_srsr1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_srsr1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_srsr1_spinor);

/* (SIGMA DOT R i SIGMA DOT R j|R12 |SIGMA DOT R k SIGMA DOT R l) */
CQ_CINT_OPT_WRAPPER(int2e_srsr1srsr2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_srsr1srsr2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_srsr1srsr2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_srsr1srsr2_spinor);

/* (RC CROSS SIGMA i SIGMA DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_cg_sa10sp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_cg_sa10sp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_cg_sa10sp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_cg_sa10sp1_spinor);

/* (RC CROSS SIGMA i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_cg_sa10sp1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_cg_sa10sp1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_cg_sa10sp1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_cg_sa10sp1spsp2_spinor);

/* (R CROSS SIGMA i SIGMA DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_giao_sa10sp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_giao_sa10sp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_giao_sa10sp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_giao_sa10sp1_spinor);

/* (R CROSS SIGMA i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_giao_sa10sp1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_giao_sa10sp1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_giao_sa10sp1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_giao_sa10sp1spsp2_spinor);

/* (G i j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_g1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_g1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_g1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_g1_spinor);

/* (G SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_spgsp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spgsp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spgsp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spgsp1_spinor);

/* (G i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_g1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_g1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_g1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_g1spsp2_spinor);

/* (G SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_spgsp1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spgsp1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spgsp1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spgsp1spsp2_spinor);

/* (P* i DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_pp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_pp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_pp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_pp1_spinor);

/* (i j|R12 |P* k DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_pp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_pp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_pp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_pp2_spinor);

/* (P* i DOT P j|R12 |P* k DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_pp1pp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_pp1pp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_pp1pp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_pp1pp2_spinor);

/* <SIGMA DOT P i|OVLP |SIGMA DOT P SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_spspsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spspsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spspsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spspsp_spinor);

/* <SIGMA DOT P i|NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_spnuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_spnuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_spnuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_spnuc_spinor);

/* (SIGMA DOT P i j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_spv1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spv1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spv1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spv1_spinor);

/* (i SIGMA DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_vsp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1_spinor);

/* (i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spsp2_spinor);

/* (SIGMA DOT P i j|R12 |SIGMA DOT P k l) */
CQ_CINT_OPT_WRAPPER(int2e_spv1spv2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spv1spv2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spv1spv2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spv1spv2_spinor);

/* (i SIGMA DOT P j|R12 |SIGMA DOT P k l) */
CQ_CINT_OPT_WRAPPER(int2e_vsp1spv2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1spv2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1spv2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1spv2_spinor);

/* (SIGMA DOT P i j|R12 |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_spv1vsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spv1vsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spv1vsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spv1vsp2_spinor);

/* (i SIGMA DOT P j|R12 |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_vsp1vsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1vsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1vsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1vsp2_spinor);

/* (SIGMA DOT P i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_spv1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_spv1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_spv1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_spv1spsp2_spinor);

/* (i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_vsp1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_vsp1spsp2_spinor);

/* <NABLA i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipovlp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipovlp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipovlp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipovlp_spinor);

/* <i|OVLP |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ovlpip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ovlpip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ovlpip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ovlpip_spinor);

/* <NABLA i|OVLP |P DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_ipkin_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipkin_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipkin_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipkin_spinor);

/* <i|OVLP |P DOT P NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_kinip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_kinip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_kinip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_kinip_spinor);

/* <NABLA i|NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipnuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipnuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipnuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipnuc_spinor);

/* <NABLA i|RINV |j> */
CQ_CINT_OPT_WRAPPER(int1e_iprinv_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_iprinv_cart);
CQ_CINT_INT1E_WRAPPER(int1e_iprinv_sph);
CQ_CINT_INT1E_WRAPPER(int1e_iprinv_spinor);

/* <NABLA SIGMA DOT P i|NUC |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_ipspnucsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipspnucsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipspnucsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipspnucsp_spinor);

/* <NABLA SIGMA DOT P i|RINV |SIGMA DOT P j> */
CQ_CINT_OPT_WRAPPER(int1e_ipsprinvsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipsprinvsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipsprinvsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipsprinvsp_spinor);

/* <P* NABLA i|NUC DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ippnucp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ippnucp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ippnucp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ippnucp_spinor);

/* <P* NABLA i|RINV DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipprinvp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipprinvp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipprinvp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipprinvp_spinor);

/* (NABLA i j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip1_spinor);

/* (i j|R12 |NABLA k l) */
CQ_CINT_OPT_WRAPPER(int2e_ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip2_spinor);

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipspsp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipspsp1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipspsp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipspsp1_spinor);

/* (NABLA i j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_ip1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip1spsp2_spinor);

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |SIGMA DOT P k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_ipspsp1spsp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipspsp1spsp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipspsp1spsp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipspsp1spsp2_spinor);

/* (NABLA SIGMA DOT R i SIGMA DOT R j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipsrsr1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipsrsr1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipsrsr1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipsrsr1_spinor);

/* (NABLA i j|R12 |SIGMA DOT R k SIGMA DOT R l) */
CQ_CINT_OPT_WRAPPER(int2e_ip1srsr2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip1srsr2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip1srsr2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip1srsr2_spinor);

/* (NABLA SIGMA DOT R i SIGMA DOT R j|R12 |SIGMA DOT R k SIGMA DOT R l) */
CQ_CINT_OPT_WRAPPER(int2e_ipsrsr1srsr2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipsrsr1srsr2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipsrsr1srsr2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipsrsr1srsr2_spinor);

/* (i SIGMA DOT P j|GAUNT |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_ssp1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ssp1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ssp1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ssp1ssp2_spinor);

/* (i SIGMA DOT P j|GAUNT |SIGMA DOT P k l) */
CQ_CINT_OPT_WRAPPER(int2e_ssp1sps2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ssp1sps2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ssp1sps2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ssp1sps2_spinor);

/* (SIGMA DOT P i j|GAUNT |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_sps1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_sps1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_sps1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_sps1ssp2_spinor);

/* (SIGMA DOT P i j|GAUNT |SIGMA DOT P k l) */
CQ_CINT_OPT_WRAPPER(int2e_sps1sps2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_sps1sps2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_sps1sps2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_sps1sps2_spinor);

/* (RC CROSS SIGMA i j|GAUNT |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_cg_ssa10ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_cg_ssa10ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_cg_ssa10ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_cg_ssa10ssp2_spinor);

/* (R CROSS SIGMA i j|GAUNT |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_giao_ssa10ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_giao_ssa10ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_giao_ssa10ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_giao_ssa10ssp2_spinor);

/* (G i SIGMA DOT P j|GAUNT |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_gssp1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gssp1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gssp1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gssp1ssp2_spinor);

/* (i R0 SIGMA DOT P j|BREIT-R1 |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r1_ssp1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_ssp1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_ssp1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_ssp1ssp2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r1_ssp1ssp2_cart);
CQ_INT_WRAPPER(int2e_gauge_r1_ssp1ssp2_sph);

/* (i R0 SIGMA DOT P j|BREIT-R1 |SIGMA DOT P k l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r1_ssp1sps2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_ssp1sps2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_ssp1sps2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_ssp1sps2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r1_ssp1sps2_sph);
CQ_INT_WRAPPER(int2e_gauge_r1_ssp1sps2_cart);


/* (SIGMA DOT P i R0 j|BREIT-R1 |k SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r1_sps1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_sps1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_sps1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_sps1ssp2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r1_sps1ssp2_cart);
CQ_INT_WRAPPER(int2e_gauge_r1_sps1ssp2_sph);

/* (SIGMA DOT P i R0 j|BREIT-R1 |SIGMA DOT P k l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r1_sps1sps2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_sps1sps2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_sps1sps2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r1_sps1sps2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r1_sps1sps2_cart);
CQ_INT_WRAPPER(int2e_gauge_r1_sps1sps2_sph);

/* (i SIGMA DOT P j|BREIT-R2 |k R0 SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r2_ssp1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_ssp1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_ssp1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_ssp1ssp2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r2_ssp1ssp2_cart);
CQ_INT_WRAPPER(int2e_gauge_r2_ssp1ssp2_sph);

/* (i SIGMA DOT P j|BREIT-R2 |SIGMA DOT P k R0 l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r2_ssp1sps2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_ssp1sps2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_ssp1sps2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_ssp1sps2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r2_ssp1sps2_cart);
CQ_INT_WRAPPER(int2e_gauge_r2_ssp1sps2_sph);

/* (SIGMA DOT P i j|BREIT-R2 |k R0 SIGMA DOT P l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r2_sps1ssp2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_sps1ssp2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_sps1ssp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_sps1ssp2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r2_sps1ssp2_cart);
CQ_INT_WRAPPER(int2e_gauge_r2_sps1ssp2_sph);

/* (SIGMA DOT P i j|BREIT-R2 |SIGMA DOT P k R0 l) */
CQ_CINT_OPT_WRAPPER(int2e_gauge_r2_sps1sps2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_sps1sps2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_sps1sps2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_gauge_r2_sps1sps2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r2_sps1sps2_cart);
CQ_INT_WRAPPER(int2e_gauge_r2_sps1sps2_sph);

/* <NABLA NABLA i|OVLP |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipipovlp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipipovlp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipipovlp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipipovlp_spinor);

/* <NABLA i|OVLP |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ipovlpip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipovlpip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipovlpip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipovlpip_spinor);

/* <NABLA NABLA i|P DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipipkin_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipipkin_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipipkin_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipipkin_spinor);

/* <NABLA i|P DOT P |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ipkinip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipkinip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipkinip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipkinip_spinor);

/* <NABLA NABLA i|NUC |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipipnuc_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipipnuc_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipipnuc_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipipnuc_spinor);

/* <NABLA i|NUC |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ipnucip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipnucip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipnucip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipnucip_spinor);

/* <NABLA NABLA i|RINV |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipiprinv_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipiprinv_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipiprinv_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipiprinv_spinor);

/* <NABLA i|RINV |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_iprinvip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_iprinvip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_iprinvip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_iprinvip_spinor);

/* (NABLA NABLA i j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipip1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipip1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipip1_spinor);

/* (NABLA i NABLA j|R12 |k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipvip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipvip1_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipvip1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_pp1_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipvip1_spinor);
CQ_INT_WRAPPER(int2e_ipvip1_cart);
CQ_INT_WRAPPER(int2e_ipvip1_sph);
CQ_INT_WRAPPER(int2e_pp1_sph);

/* (NABLA i j|R12 |NABLA k l) */
CQ_CINT_OPT_WRAPPER(int2e_ip1ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ip1ip2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ip1ip2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ip1ip2_spinor);
CQ_INT_WRAPPER(int2e_ip1ip2_cart);
CQ_INT_WRAPPER(int2e_ip1ip2_sph);

/* <P* NABLA NABLA i|NUC DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipippnucp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipippnucp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipippnucp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipippnucp_spinor);

/* <P* NABLA i|NUC DOT P |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ippnucpip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ippnucpip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ippnucpip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ippnucpip_spinor);

/* <P* NABLA NABLA i|RINV DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipipprinvp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipipprinvp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipipprinvp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipipprinvp_spinor);

/* <P* NABLA i|RINV DOT P |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ipprinvpip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipprinvpip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipprinvpip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipprinvpip_spinor);

/* <NABLA NABLA SIGMA DOT P i|NUC SIGMA DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipipspnucsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipipspnucsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipipspnucsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipipspnucsp_spinor);

/* <NABLA SIGMA DOT P i|NUC SIGMA DOT P |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ipspnucspip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipspnucspip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipspnucspip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipspnucspip_spinor);

/* <NABLA NABLA SIGMA DOT P i|RINV SIGMA DOT P |j> */
CQ_CINT_OPT_WRAPPER(int1e_ipipsprinvsp_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipipsprinvsp_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipipsprinvsp_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipipsprinvsp_spinor);

/* <NABLA SIGMA DOT P i|RINV SIGMA DOT P |NABLA j> */
CQ_CINT_OPT_WRAPPER(int1e_ipsprinvspip_optimizer);
CQ_CINT_INT1E_WRAPPER(int1e_ipsprinvspip_cart);
CQ_CINT_INT1E_WRAPPER(int1e_ipsprinvspip_sph);
CQ_CINT_INT1E_WRAPPER(int1e_ipsprinvspip_spinor);

/* (NABLA NABLA i j|R12 |NABLA NABLA k l) */
CQ_CINT_OPT_WRAPPER(int2e_ipip1ipip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipip1ipip2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipip1ipip2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipip1ipip2_spinor);

/* (NABLA i NABLA j|R12 |NABLA k NABLA l) */
CQ_CINT_OPT_WRAPPER(int2e_ipvip1ipvip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2e_ipvip1ipvip2_cart);
CQ_CINT_INT2E_WRAPPER(int2e_ipvip1ipvip2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_pp1pp2_sph);
CQ_CINT_INT2E_WRAPPER(int2e_ipvip1ipvip2_spinor);
CQ_INT_WRAPPER(int2e_ipvip1ipvip2_cart);
CQ_INT_WRAPPER(int2e_ipvip1ipvip2_sph);
CQ_INT_WRAPPER(int2e_pp1pp2_sph);

/* (NABLA i j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip1_spinor);

/* (i j|R12 |NABLA k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip2_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip2_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip2_spinor);

/* (P* i DOT P j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_pvp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_pvp1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_pvp1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_pvp1_spinor);

/* (P* i CROSS P j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_pvxp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_pvxp1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_pvxp1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_pvxp1_spinor);

/* (NABLA i |R12 |j) */
CQ_CINT_OPT_WRAPPER(int2c2e_ip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip1_cart);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip1_sph);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip1_spinor);

/* (i |R12 |NABLA j) */
CQ_CINT_OPT_WRAPPER(int2c2e_ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip2_cart);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip2_sph);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip2_spinor);

/* (G i j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ig1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ig1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ig1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ig1_spinor);

/* (SIGMA DOT P i SIGMA DOT P j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_spsp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_spsp1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_spsp1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_spsp1_spinor);

/* (NABLA SIGMA DOT P i SIGMA DOT P j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ipspsp1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipspsp1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipspsp1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipspsp1_spinor);

/* (SIGMA DOT P i SIGMA DOT P j|R12 |NABLA k) */
CQ_CINT_OPT_WRAPPER(int3c2e_spsp1ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_spsp1ip2_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_spsp1ip2_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_spsp1ip2_spinor);

/* (NABLA NABLA i j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ipip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipip1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipip1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipip1_spinor);

/* (i j|R12 |NABLA NABLA k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ipip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipip2_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipip2_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipip2_spinor);

/* (NABLA i NABLA j|R12 |k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ipvip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipvip1_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipvip1_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ipvip1_spinor);

/* (NABLA i j|R12 |NABLA k) */
CQ_CINT_OPT_WRAPPER(int3c2e_ip1ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip1ip2_cart);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip1ip2_sph);
CQ_CINT_INT2E_WRAPPER(int3c2e_ip1ip2_spinor);

/* (NABLA NABLA i |R12 |j) */
CQ_CINT_OPT_WRAPPER(int2c2e_ipip1_optimizer);
CQ_CINT_INT2E_WRAPPER(int2c2e_ipip1_cart);
CQ_CINT_INT2E_WRAPPER(int2c2e_ipip1_sph);
CQ_CINT_INT2E_WRAPPER(int2c2e_ipip1_spinor);

/* (NABLA i |R12 |NABLA j) */
CQ_CINT_OPT_WRAPPER(int2c2e_ip1ip2_optimizer);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip1ip2_cart);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip1ip2_sph);
CQ_CINT_INT2E_WRAPPER(int2c2e_ip1ip2_spinor);

/* 3-center 1-electron integral <(i) (j) (P DOT P k)> */
CQ_CINT_OPT_WRAPPER(int3c1e_p2_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c1e_p2_cart);
CQ_CINT_INT2E_WRAPPER(int3c1e_p2_sph);
CQ_CINT_INT2E_WRAPPER(int3c1e_p2_spinor);

/* 3-center 1-electron integral <(P i) (j) (k)> */
CQ_CINT_OPT_WRAPPER(int3c1e_iprinv_optimizer);
CQ_CINT_INT2E_WRAPPER(int3c1e_iprinv_cart);
CQ_CINT_INT2E_WRAPPER(int3c1e_iprinv_sph);
CQ_CINT_INT2E_WRAPPER(int3c1e_iprinv_spinor);

/* spin-free Gaunt */
CQ_INT_WRAPPER(int2e_gaunt_ps1ps2_cart);
CQ_INT_WRAPPER(int2e_gaunt_ps1ps2_sph);
CQ_INT_WRAPPER(int2e_gaunt_ps1ps2_spinor);

/* spin-free Gauge */
CQ_INT_WRAPPER(int2e_gauge_r1_sp1sp2_cart);
CQ_INT_WRAPPER(int2e_gauge_r1_sp1sp2_sph);
CQ_INT_WRAPPER(int2e_gauge_r1_sp1sp2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r2_sp1sp2_cart);
CQ_INT_WRAPPER(int2e_gauge_r2_sp1sp2_sph);
CQ_INT_WRAPPER(int2e_gauge_r2_sp1sp2_spinor);

CQ_INT_WRAPPER(int2e_gauge_r1_sp1ps2_cart);
CQ_INT_WRAPPER(int2e_gauge_r1_sp1ps2_sph);
CQ_INT_WRAPPER(int2e_gauge_r1_sp1ps2_spinor);
CQ_INT_WRAPPER(int2e_gauge_r2_sp1ps2_cart);
CQ_INT_WRAPPER(int2e_gauge_r2_sp1ps2_sph);
CQ_INT_WRAPPER(int2e_gauge_r2_sp1ps2_spinor);


}  // extern C
