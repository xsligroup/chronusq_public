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
#include <cubegen/density.hpp>
#include <cubegen/orbital.hpp>
#include <cubegen/impl.hpp>

namespace ChronusQ {

    template
    void CubeGen::evalDenCube(std::string filePrefix, std::shared_ptr<cqmatrix::PauliSpinorMatrices<double>> oPDM_, double particleCharge, bool skipoutput);
    template
    void CubeGen::evalDenCube(std::string filePrefix, std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> oPDM_, double particleCharge, bool skipoutput);
    template
    void CubeGen::evalDenCompCube(double* oPDM_, double particleCharge);
    template
    void CubeGen::evalDenCompCube(dcomplex* oPDM_, double particleCharge);
    template 
    void CubeGen::evalOrbCube(std::string filePrefix, double* MOs, size_t LDMO, std::vector<size_t> whichMO, std::function<double(double)> Op);
    template 
    void CubeGen::evalOrbCube(std::string filePrefix, dcomplex* MOs, size_t LDMO, std::vector<size_t> whichMO, std::function<double(dcomplex)> Op);
    template
    void CubeGen::evalOrbCompCube(double * MO,size_t LDMO, size_t MOIndex,std::function<double(double)>Op);
    template
    void CubeGen::evalOrbCompCube(dcomplex * MO,size_t LDMO, size_t MOIndex,std::function<double(dcomplex)>Op);


}
