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
#include <singleslater/impl.hpp>
#include <corehbuilder/impl.hpp>
#include <fockbuilder/impl.hpp>
#include <mointstransformer/impl.hpp>
#include <singleslater/neoss/impl.hpp>

namespace ChronusQ {

  template class SingleSlater<double,double>;
  template class SingleSlater<dcomplex,double>;
  template class SingleSlater<dcomplex,dcomplex>;

  // Handle MatsT on scr file
  // TODO: implement IntsT=dcomplex
  template void SingleSlater<double,double>::getScr1PDM<double>(SafeFile &);
  template void SingleSlater<dcomplex,double>::getScr1PDM<double>(SafeFile &);
  template void SingleSlater<dcomplex,double>::getScr1PDM<dcomplex>(SafeFile &);

  // Explicit template instantiation for RT functions that require dcomplex matrix types
  template void SingleSlater<dcomplex, double>::addTauToFock<dcomplex>();
  template void SingleSlater<dcomplex, dcomplex>::addTauToFock<dcomplex>();
  template void SingleSlater<dcomplex, double>::RK4Propagation<dcomplex>(bool, double, bool, EMPerturbation&, EMPerturbation&);
  template void SingleSlater<dcomplex, dcomplex>::RK4Propagation<dcomplex>(bool, double, bool, EMPerturbation&, EMPerturbation&);
  template void SingleSlater<dcomplex, double>::unitaryPropagation<dcomplex>(bool, double, bool, EMPerturbation&);
  template void SingleSlater<dcomplex, dcomplex>::unitaryPropagation<dcomplex>(bool, double, bool, EMPerturbation&);
  template cqmatrix::PauliSpinorMatrices<dcomplex> SingleSlater<dcomplex, double>::getTimeDerDen<dcomplex>(bool);
  template cqmatrix::PauliSpinorMatrices<dcomplex> SingleSlater<dcomplex, dcomplex>::getTimeDerDen<dcomplex>(bool);

  // Instantiate copy constructors
  template SingleSlater<dcomplex,double>::SingleSlater(const SingleSlater<double,double> &, int);
  template SingleSlater<dcomplex,double>::SingleSlater(const SingleSlater<dcomplex,double> &, int);
  template SingleSlater<dcomplex,dcomplex>::SingleSlater(const SingleSlater<dcomplex,dcomplex> &, int);

  // Instantiate move ctors
  template SingleSlater<dcomplex,double>::SingleSlater( SingleSlater<double,double> &&, int);

  template class HartreeFock<double,double>;
  template class HartreeFock<dcomplex,double>;
  template class HartreeFock<dcomplex,dcomplex>;
  // Instantiate copy constructors
  template HartreeFock<dcomplex,double>::HartreeFock(const HartreeFock<double,double> &, int);
  // Instantiate copy ructors
  template HartreeFock<dcomplex,double>::HartreeFock( HartreeFock<double,double> &&, int);

  template class KohnSham<double,double>;
  template class KohnSham<dcomplex,double>;
  template class KohnSham<dcomplex,dcomplex>;
  // Instantiate copy constructors
  template KohnSham<dcomplex,double>::KohnSham(const KohnSham<double,double> &, int);
  // Instantiate copy ructors
  template KohnSham<dcomplex,double>::KohnSham( KohnSham<double,double> &&, int);

  template void KohnSham<double,double>::formFXC(MPI_Comm,std::vector<TwoBodyContraction<double>> &, EMPerturbation&);
  template void KohnSham<double,double>::formFXC(MPI_Comm,std::vector<TwoBodyContraction<dcomplex>> &, EMPerturbation&);

  template void KohnSham<dcomplex,double>::formFXC(MPI_Comm,std::vector<TwoBodyContraction<dcomplex>> &, EMPerturbation&);

  template class NEOSS<double,double>;
  template class NEOSS<dcomplex,double>;
  template class NEOSS<dcomplex,dcomplex>;
  template NEOSS<dcomplex,double>::NEOSS(const NEOSS<double, double>&, int);
  template NEOSS<dcomplex,double>::NEOSS(const NEOSS<dcomplex, double>&, int);
  template NEOSS<dcomplex,dcomplex>::NEOSS(const NEOSS<dcomplex, dcomplex>&, int);
  
  template NEOSS<dcomplex,double>::NEOSS(NEOSS<double, double>&&, int);

  template std::shared_ptr<SingleSlater<double,double>>
  SingleSlater<double,double>::convert(const std::shared_ptr<SingleSlater<double,double>>& ss);
  template std::shared_ptr<SingleSlater<dcomplex,double>>
  SingleSlater<double,double>::convert(const std::shared_ptr<SingleSlater<double,double>>& ss);
  template std::shared_ptr<SingleSlater<dcomplex,double>>
  SingleSlater<dcomplex,double>::convert(const std::shared_ptr<SingleSlater<dcomplex,double>>& ss);
  template std::shared_ptr<SingleSlater<dcomplex,dcomplex>>
  SingleSlater<dcomplex,dcomplex>::convert(const std::shared_ptr<SingleSlater<dcomplex,dcomplex>>& ss);


  }; // namespace ChronusQ
