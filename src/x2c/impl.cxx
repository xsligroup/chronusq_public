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
#include <corehbuilder/x2c/impl.hpp>

namespace ChronusQ {

  template class X2C<double,double>;
  template class X2C<dcomplex,double>;
  template class X2C<dcomplex,dcomplex>;

  // Instantiate copy constructors
  template X2C<dcomplex,double>::X2C(const X2C<double,double> &, int);
  template X2C<dcomplex,dcomplex>::X2C(const X2C<dcomplex,dcomplex> &, int);

}; // namespace ChronusQ
