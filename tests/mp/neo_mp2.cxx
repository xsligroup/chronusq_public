/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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

#include "mp2.hpp"

// MP2 tests, serial job
TEST(NEO_MP2,h2o_ccpvdz_pb4d_incore)
{
  CQMPTEST("mp/serial/neo_mp2/water_ccpvdz_pb4d_incore","water_ccpvdz_pb4d.bin.ref");
};

TEST(NEO_MP2,h2o_ccpvdz_pb4d_direct)
{
  CQMPTEST("mp/serial/neo_mp2/water_ccpvdz_pb4d_direct","water_ccpvdz_pb4d.bin.ref");
};

TEST(NEO_MP2,hcn_ccpvdz_pb4d)
{
  CQMPTEST("mp/serial/neo_mp2/hcn_ccpvdz_pb4d","hcn_ccpvdz_pb4d.bin.ref");
};

TEST(NEO_MP2,hcn_ccpvdz_pb4d_altref)
{
  CQMPTEST("mp/serial/neo_mp2/hcn_ccpvdz_pb4d_altref","hcn_ccpvdz_pb4d_altref.bin.ref");
};

// MP2 tests, parallel job
TEST(NEO_MP2,par_h2o_ccpvdz_pb4d_incore)
{
  CQMPTEST("mp/parallel/neo_mp2/water_ccpvdz_pb4d_incore","water_ccpvdz_pb4d.bin.ref");
};

TEST(NEO_MP2,par_h2o_ccpvdz_pb4d_direct)
{
  CQMPTEST("mp/parallel/neo_mp2/water_ccpvdz_pb4d_direct","water_ccpvdz_pb4d.bin.ref");
};

TEST(NEO_MP2,par_hcn_ccpvdz_pb4d)
{
  CQMPTEST("mp/parallel/neo_mp2/hcn_ccpvdz_pb4d","hcn_ccpvdz_pb4d.bin.ref");
};

TEST(NEO_MP2,par_hcn_ccpvdz_pb4d_altref)
{
  CQMPTEST("mp/parallel/neo_mp2/hcn_ccpvdz_pb4d_altref","hcn_ccpvdz_pb4d_altref.bin.ref");
};

