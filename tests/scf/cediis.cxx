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

#include "scf.hpp"
#include "cediis/bfgs-test.hpp"
#include "cediis/goldenSection-test.hpp"
#include "cediis/interpolate-test.hpp"


// Test of BFGS algorithm using a simple surface opt
TEST( CEDIIS, BFGS_Solver_Test ) {

  {
    CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,256e6,2048);
    BFGSTest bfgsTest;
    ASSERT_TRUE( bfgsTest.converged ) << "BFGS Solver Test did not converge";
    EXPECT_NEAR(bfgsTest.x[0],  1., 1E-8) << "BFGS Solver Test Failed for first variable";
    EXPECT_NEAR(bfgsTest.x[1],  0., 1E-8) << "BFGS Solver Test Failed for second variable";
    EXPECT_NEAR(bfgsTest.gx[0], 0., 1E-8) << "BFGS Solver Test Failed for gradient[0]. ie. gradient is not zero";
    EXPECT_NEAR(bfgsTest.gx[1], 0., 1E-8) << "BFGS Solver Test Failed for gradient[1]. ie. gradient is not zero";
    EXPECT_NEAR(bfgsTest.yMin, -1., 1E-8) << "BFGS Solver Test Failed for function value";
  }
 
};

// Test of GoldenSectionSearch Algorithm using 1D parabola
TEST( CEDIIS, GoldenSection_Solver_Test ) {

  {
    CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,256e6,2048);
    GoldenSectionTest goldenSectionTest;
    EXPECT_NEAR(goldenSectionTest.x[0], 0., 1E-8) << "Golden Section Search Solver Test Failed for variable";
    EXPECT_NEAR(goldenSectionTest.yMin, 0., 1E-8) << "Golden Section Search Solver Test Failed for ordinate";
  }
 
};

// Test of Interpolate optimization using a simple 2D case
TEST( CEDIIS, Interpolate_2D_Test ) {

  {
    CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,256e6,2048);
    InterpolateTest interpolateTest;
    ASSERT_TRUE( interpolateTest.converged ) << "Interpolate Test did not converge";
    EXPECT_NEAR(interpolateTest.coeffs[0], 0.75, 1E-6) << "Interpolate Test Failed for first coefficient";
    EXPECT_NEAR(interpolateTest.coeffs[1], 0.25, 1E-6) << "Interpolate Test Failed for second coefficient";
    EXPECT_NEAR(interpolateTest.gOpt[0], 0., 1E-6) << "Interpolate Test Failed for gradient[0]. i.e. gradient is not zero";
    EXPECT_NEAR(interpolateTest.gOpt[1], 0., 1E-6) << "Interpolate Test Failed for gradient[1]. i.e. gradient is not zero";
  }
 
};

// Acetaldehyde/CEDIIS/RHF/Sapporo-DZP-no test for convergence
TEST( CEDIIS, Acetaldehyde_CEDIIS_RHF ) {

  CQSCFTEST( "scf/serial/cediis/acetaldehyde-cediis-rhf", "acetaldehyde-cediis-rhf.bin.ref" );
 
};

// Potassium/CEDIIS/ROHF/6-31G(d) test for convergence
TEST( CEDIIS, Potassium_CEDIIS_ROHF ) {

  CQSCFTEST( "scf/serial/cediis/k_6-31Gd_cediis_ro", "k_6-31Gd_cediis_ro.bin.ref" );
 
};

// Acetaldehyde/CEDIIS/UHF/Sapporo-DZP-no test for convergence
TEST( CEDIIS, Acetaldehyde_CEDIIS_UHF ) {

  CQSCFTEST( "scf/serial/cediis/acetaldehyde-cediis-uhf", "acetaldehyde-cediis-uhf.bin.ref" );
 
};

// Acetaldehyde/CEDIIS/GHF/Sapporo-DZP-no test for convergence
TEST( CEDIIS, Acetaldehyde_CEDIIS_GHF ) {

  CQSCFTEST( "scf/serial/cediis/acetaldehyde-cediis-ghf", "acetaldehyde-cediis-ghf.bin.ref" );
 
};

// Silver/CEDIIS/X2C/Sapporo-DZP test for convergence
TEST( CEDIIS, Ag_CEDIIS_X2C ) {

  CQSCFTEST( "scf/serial/cediis/Ag-cediis-x2c", "Ag-cediis-x2c.bin.ref" ,1e-8,
    false, false, false, false, false, true );
 
};

// Silver/CEDIIS/4C/Sapporo-DZP test for convergence
TEST( CEDIIS, Ag_CEDIIS_4C ) {

  CQSCFTEST( "scf/serial/cediis/Ag-cediis-dcg", "Ag-cediis-dcg.bin.ref",1e-8,
    false, false, false, false, false, true );
 
};

// Acetaldehyde/EDIIS/RHF/Sapporo-DZP-no test for convergence
TEST( CEDIIS, Acetaldehyde_EDIIS_RHF ) {

  CQSCFTEST( "scf/serial/cediis/acetaldehyde-ediis-rhf", "acetaldehyde-ediis-rhf.bin.ref" ,1e-8,
    false, false, false, false, false, true );
 
};

// Potassium/EDIIS/ROHF/6-31G(d) test for convergence
TEST( CEDIIS, Potassium_EDIIS_ROHF ) {

  CQSCFTEST( "scf/serial/cediis/k_6-31Gd_ediis_ro", "k_6-31Gd_ediis_ro.bin.ref" ,1e-8,
    false, false, false, false, false, true );
 
};

// Acetaldehyde/EDIIS/UHF/Sapporo-DZP-no test for convergence
TEST( CEDIIS, Acetaldehyde_EDIIS_UHF ) {

  CQSCFTEST( "scf/serial/cediis/acetaldehyde-ediis-uhf", "acetaldehyde-ediis-uhf.bin.ref" ,1e-8,
    false, false, false, false, false, true );
 
};

// Acetaldehyde/EDIIS/GHF/Sapporo-DZP-no test for convergence
TEST( CEDIIS, Acetaldehyde_EDIIS_GHF ) {

  CQSCFTEST( "scf/serial/cediis/acetaldehyde-ediis-ghf", "acetaldehyde-ediis-ghf.bin.ref" ,1e-8,
    false, false, false, false, false, true );
 
};


#ifdef _CQ_DO_PARTESTS


#endif



