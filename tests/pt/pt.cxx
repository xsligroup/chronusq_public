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

#include "pt.hpp"


TEST( MRPT, AL_631G_X2C) {
  CQPTTEST("pt/serial/al_6-31G_x2c",
    "al_6-31G_x2c.bin.ref", false);
}

TEST( MRPT, H2O_631G_X2C) {
  CQPTTEST("pt/serial/water_631g_x2c_readci",
    "water_631g_x2c_readci.bin.ref", false, "water_631g_x2c_ref.bin.ref");
}

#ifdef _CQ_DO_PARTESTS

TEST( MRPT, PAR_AL_631G_X2C) {
  CQPTTEST("pt/parallel/al_6-31G_x2c",
           "al_6-31G_x2c.bin.ref", false);
}

TEST( MRPT, PAR_H2O_631G_X2C) {
  CQPTTEST("pt/parallel/water_631g_x2c_readci",
    "water_631g_x2c_readci.bin.ref", false, "water_631g_x2c_ref.bin.ref");
}
#endif

