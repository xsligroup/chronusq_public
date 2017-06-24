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
#include <realtime/realtimesingleslater/fields.hpp>

// This file instantiates the various field envelop definitions
namespace ChronusQ {

  /**
   *  \brief Definition of a step function field envelope
   *  
   *  \f[
   *    F(t) = \begin{cases}
   *              1 &  \mathrm{tOn} <= t <= \mathrm{tOff} \\
   *              0 &  \mathrm{else}
   *           \end{cases}
   *  \f]
   *
   *  Also handle the case where t is close to the time
   *  boundaries due to numerical precision.
   *
   *  \param [in] t Time point
   *  \returns      F(t)
   *
   */ 
  double StepField::getAmp(double t) {
    const double epsilonDeltaT = 1e-10;
    if(
        ((t >= this->tOn - epsilonDeltaT) or (t >= this->tOn + epsilonDeltaT))
        and 
        ((t <= this->tOff - epsilonDeltaT) or (t <= this->tOff + epsilonDeltaT))
      ) {
        return 1.;
    } else return 0.;
  };

  /**
   *  \brief Definition of a linear function field envelope
   *  
   *  \f[
   *    F(t) = \begin{cases}
   *              t - tOn / (tOff - tOn) &  \mathrm{tOn} <= t <= \mathrm{tOff} \\
   *              0 &  \mathrm{else}
   *           \end{cases}
   *  \f]
   *
   *  Also handle the case where t is close to the time
   *  boundaries due to numerical precision.
   *
   *  \param [in] t Time point
   *  \returns      F(t)
   *
   */ 
  double LinRampField::getAmp(double t) {
    const double epsilonDeltaT = 1e-10;
    if(
        ((t >= this->tOn - epsilonDeltaT) or (t >= this->tOn + epsilonDeltaT))
        and 
        ((t <= this->tOff - epsilonDeltaT) or (t <= this->tOff + epsilonDeltaT))
      ) {
        return (t - this->tOn) / (this->tOff - this->tOn);
    } else {
        return 0.;
    }
  };

  /**
   *  \brief Definition of a Gaussian function field envelope
   *  
   *  \f[
   *    F(t) = \begin{cases}
   *              exp(-\alpha (t - tOn)^2) &  \mathrm{tOn} <= t <= \mathrm{tOff} \\
   *              0 &  \mathrm{else}
   *           \end{cases}
   *  \f]
   *  
   *  Note the negative sign is included in getAmp so $\alpha$ should be positive (unless you want an exponentially growing field).
   *
   *  Also handle the case where t is close to the time
   *  boundaries due to numerical precision.
   *
   *  \param [in] t Time point
   *  \returns      F(t)
   *
   */ 
  double GaussianField::getAmp(double t) {
    const double epsilonDeltaT = 1e-10;
    if(
        ((t >= this->tOn - epsilonDeltaT) or (t >= this->tOn + epsilonDeltaT))
        and 
        ((t <= this->tOff - epsilonDeltaT) or (t <= this->tOff + epsilonDeltaT))
      ) {
        return std::exp( - this->alpha * std::pow(t - this->tOn, 2.));
    } else {
        return 0.;
    }
  };

/**
   *  \brief Definition of a Plane Wave function field envelope
   *  
   *  \f[
   *    F(t) = \begin{cases}
   *              sin/cos(\omega (t - tOn)) &  \mathrm{tOn} <= t <= \mathrm{tOff} \\
   *              0 &  \mathrm{else}
   *           \end{cases}
   *  \f]
   *  
   *  Also handle the case where t is close to the time
   *  boundaries due to numerical precision.
   *
   *  \param [in] t Time point
   *  \returns      F(t)
   *
   */ 
  double PlaneWaveField::getAmp(double t) {
    const double epsilonDeltaT = 1e-10;
    if(
        ((t >= this->tOn - epsilonDeltaT) or (t >= this->tOn + epsilonDeltaT))
        and 
        ((t <= this->tOff - epsilonDeltaT) or (t <= this->tOff + epsilonDeltaT))
      ) {
        if (this->doCos)
            return std::cos( this->omega * (t - this->tOn));
        else
            return std::sin( this->omega * (t - this->tOn));
    } else {
        return 0.;
    }
  };


}; // namespace ChronusQ
