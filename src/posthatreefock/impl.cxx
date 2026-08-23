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
#include <posthartreefock/impl.hpp>

namespace ChronusQ {

template class PostHartreeFock<double,double>;
template class PostHartreeFock<dcomplex,double>;
template class PostHartreeFock<dcomplex,dcomplex>;

// Instantiate copy constructors
template PostHartreeFock<dcomplex,double>::PostHartreeFock(const PostHartreeFock<double,double> &, int);
template PostHartreeFock<dcomplex,dcomplex>::PostHartreeFock(const PostHartreeFock<dcomplex,dcomplex> &, int);

// Instantiate move ctors
template PostHartreeFock<dcomplex,double>::PostHartreeFock(PostHartreeFock<double,double> &&, int);

template <>
void MOIntsTransformer<dcomplex,dcomplex>::directTransformTPIBatch(EMPerturbation & pert,
  dcomplex* MOTPI, const std::vector<std::pair<size_t,size_t>> & off_sizes) const {
  CErr("Complex integral is an invalid option",std::cout);
}

template
void MOIntsTransformer<double,double>::directTransformTPIBatch(EMPerturbation & pert,
  double* MOTPI, const std::vector<std::pair<size_t,size_t>> & off_sizes) const;

template
void MOIntsTransformer<dcomplex,double>::directTransformTPIBatch(EMPerturbation & pert,
  dcomplex* MOTPI, const std::vector<std::pair<size_t,size_t>> & off_sizes) const;



} // namespace ChronusQ
