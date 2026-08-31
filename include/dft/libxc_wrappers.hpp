/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

#include <dft.hpp>

namespace ChronusQ {

  class SlaterExchange : public LDA {

  public:
  
    SlaterExchange() : LDA(xc_functional_get_number("LDA_X")) { }

  }; // class SlaterExchage


  class VWNIII : public LDA {

  public:

    //VWNIII() : LDA(XC_LDA_C_VWN_RPA) { }
    VWNIII() : LDA(xc_functional_get_number("LDA_C_VWN_3")) { }


  }; // class VWNIII


  class VWNV : public LDA {

  public:

    VWNV() : LDA(xc_functional_get_number("LDA_C_VWN_RPA")) { }

  }; // class VWNV



  class VWNV_G : public LDA {

  public:

    VWNV_G() : LDA(xc_functional_get_number("LDA_C_VWN")) { }

  }; // class VWNV_GAUSSIAN

  class BEightyEight : public GGA {

  public:
  
    BEightyEight() : GGA(xc_functional_get_number("GGA_X_B88")) { }

  }; // class B88

  class PBEX : public GGA {

  public:
  
    PBEX() : GGA(xc_functional_get_number("GGA_X_PBE")) { }

  }; // class PBE exchange



  class LYP : public GGA {

  public:
  
    LYP() : GGA(xc_functional_get_number("GGA_C_LYP")) { }

  }; // class B88

  class PBEC : public GGA {

  public:
  
    PBEC() : GGA(xc_functional_get_number("GGA_C_PBE")) { }

  }; // class PBE correlation

  class B3LYP : public GGA {

  public:
  
    B3LYP() : GGA(xc_functional_get_number("HYB_GGA_XC_B3LYP")) { }

  }; // class B3LYP hybrid

  class B3PW91 : public GGA {

  public:
  
    B3PW91() : GGA(xc_functional_get_number("HYB_GGA_XC_B3PW91")) { }

  }; // class B3LYP hybrid

  class PBE0 : public GGA {

  public:
  
    PBE0() : GGA(xc_functional_get_number("HYB_GGA_XC_PBEH")) { }

  }; // class PBE0 hybrid

  class BHANDH : public GGA {

  public:
  
    BHANDH() : GGA(xc_functional_get_number("HYB_GGA_XC_BHANDH")) { }

  }; // class BHANDH hybrid

  class BHANDHLYP : public GGA {

  public:
  
    BHANDHLYP() : GGA(xc_functional_get_number("HYB_GGA_XC_BHANDHLYP")) { }

  }; // class BHANDHLYP hybrid

} // namespace ChronusQ
