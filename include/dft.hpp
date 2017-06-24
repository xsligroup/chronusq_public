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

#include <chronusq_sys.hpp>
#include <xc.h>

namespace ChronusQ {

  /**
   *  \brief Implements a DFT functional
   */
  class DFTFunctional {
    // (polarized) rho[2*np], vrho[2*np], sigma[3*np], vsigma[3*np], v2rho2[3*np], v2rhosigma[6*np], v2sigma2[6*np]
    //  in rho [0]    = UP   rho 
    //  in rho [1]    = DOWN rho
    //  in sigma[0]   = UP   - UP   Del N dot Del N 
    //  in sigma[1]   = UP   - DOWN Del N dot Del N 
    //  in sigma[2]   = DOWN - DOWN Del N dot Del N 
    //  out eps       = the energy per unit particle
    //  out vxc/vrho[0]   = UP   first part der of the energy per unit volume in terms of the dens 
    //  out vxc/vrho[1]   = DOWN first part der of the energy per unit volume in terms of the dens N
    //  out vsigma[0]  = UP   UP    first part der of the energy per unit volume in terms of sigma 
    //  out vsigma[1]  = UP   DOWN  first part der of the energy per unit volume in terms of sigma 
    //  out vsigma[2]  = DOWN DOWN  first part der of the energy per unit volume in terms of sigma 
    //  out fxc/v2rho2[0]  = UP UP      second part der of the energy per unit volume in terms of the dens 
    //  out fxc/v2rho2[1]  = UP UP      second part der of the energy per unit volume in terms of the dens 
    //  out fxc/v2rho2[2]  = UP UP      second part der of the energy per unit volume in terms of the dens 
    //  out v2rhosigma[0]  = UP   - UP UP     second part der of the energy per unit volume in terms of the dens and sigma 
    //  out v2rhosigma[1]  = UP   - UP DOWN   second part der of the energy per unit volume in terms of the dens and sigma 
    //  out v2rhosigma[2]  = UP   - DOWN DOWN second part der of the energy per unit volume in terms of the dens and sigma 
    //  out v2rhosigma[3]  = DOWN - UP UP     second part der of the energy per unit volume in terms of the dens and sigma 
    //  out v2rhosigma[4]  = DOWN - UP DOWN   second part der of the energy per unit volume in terms of the dens and sigma 
    //  out v2rhosigma[5]  = DOWN - DOWN DOWN second part der of the energy per unit volume in terms of the dens and sigma 
    //  out v2sigma2[0]    = UP UP     - UP UP     second part der of the energy per unit volume in terms of sigma sigma
    //  out v2sigma2[1]    = UP UP     - UP DOWN   second part der of the energy per unit volume in terms of sigma sigma
    //  out v2sigma2[2]    = UP UP     - DOWN DOWN second part der of the energy per unit volume in terms of sigma sigma
    //  out v2sigma2[3]    = UP DONW   - UP DONW   second part der of the energy per unit volume in terms of sigma sigma
    //  out v2sigma2[4]    = UP DONW   - DOWN DOWN second part der of the energy per unit volume in terms of sigma sigma
    //  out v2sigma2[5]    = DOWN DOWN - DOWN DOWN second part der of the energy per unit volume in terms of sigma sigma

  protected:
    bool isGGA_;

    bool isEPC_ = false;
    

    xc_func_type functional_; ///< LibXC functional definition

  public:

    double xHFX; ///< Scaling factor for HFX

    DFTFunctional(int FUNC_IDENT) { 
      xc_func_init(&functional_,FUNC_IDENT,XC_POLARIZED);
      xHFX = xc_hyb_exx_coef(&functional_);
      //std::cerr << "HYB " << xHFX << std::endl;
    };

    // Constructor specifically for EPC functional
    DFTFunctional(std::string epcName) { xHFX = 1.; };

    // TODO: Implement COPY / MOVE...
      
    void evalEXC_VXC(size_t N, double *rho, double *eps, double *vxc) {
      assert(not isGGA_);
      xc_lda_exc_vxc(&this->functional_,N,rho,eps,vxc);
    }

    void evalEXC_VXC(size_t N, double *rho, double *sigma, double *eps, double *vrho, 
      double * vsigma) {

      assert(isGGA_);
      xc_gga_exc_vxc(&this->functional_,N,rho,sigma,eps,vrho,vsigma);

    };

    virtual void evalEXC_VXC(size_t N, double *rho, double *aux_rho, double *eps, 
                             double *vxc, bool electron) {}; 

    virtual void evalEXC_VXC(size_t N, double *rho, double *aux_rho, double *sigma,
                             double *aux_sigma, double *cross_sigma, double *eps,
                             double *vrho, double *vsigma, double *vcsigma, 
                             bool electron) {};


    void evalEXC_VXC_FXC(size_t N, double *rho, double *eps, double *vxc, double *fxc) {

      evalEXC_VXC(N,rho,eps,vxc);
      xc_lda_fxc(&this->functional_,N,rho,fxc);

    }

    void evalEXC_VXC_FXC(size_t N, double *rho, double *sigma, double *eps, double *vrho, 
      double * vsigma, double *v2rho2, double *v2rhosigma, double *v2sigma2) {

      evalEXC_VXC(N,rho,sigma,eps,vrho,vsigma);
      xc_gga_fxc(&this->functional_,N,rho,sigma,v2rho2,v2rhosigma,v2sigma2);

    };


/*
    virtual void evalFXC(size_t N, double *rho, double *fxc) = 0;
    virtual void evalFXC(size_t N, double *rho, double *sigma, double *v2rho2, double *v2rhosigma, double *v2sigma2) = 0;
*/

    bool isGGA() { return this->isGGA_; }

    bool isEPC() { return this->isEPC_; }

  }; // class DFTFunctional



  /**
   *  \brief Base class for all LDA functionals
   */ 
  class LDA : public DFTFunctional {

  public:
    LDA(int FUNC_IDENT) : DFTFunctional(FUNC_IDENT) { this->isGGA_ = false;};

  }; // class LDA



  /**
   *  \brief Base class for all GGA functionals
   */ 
  class GGA : public DFTFunctional {

  public:
    GGA(int FUNC_IDENT) : DFTFunctional(FUNC_IDENT) { this->isGGA_ = true;};
    
  }; // class GGA

};

#include <dft/libxc_wrappers.hpp>
#include <dft/epc.hpp>
