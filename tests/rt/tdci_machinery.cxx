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

#include <ut.hpp>

#include <cerr.hpp>
#include <memmanager.hpp>
#include <util/matout.hpp>
#include <realtime.hpp>

using namespace ChronusQ;

void tdci_machinery_TEST_real() {

  std::cout << "**** Test TDCI Infrastructure ****" << std::endl;
  
  const double TDCI_TEST_EPSILON = 1e-12;
  size_t mem     = 256e6; // Default 256 MB allocation
  size_t blkSize = 2048;  // Default 2KB block size

  typedef double MatsT;

  CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,mem,blkSize);
  
  size_t vec_size = 5;
  std::shared_ptr<SolverVectors<MatsT>> a,b,c;
  a = std::make_shared<RawVectors<MatsT>>(MPI_COMM_WORLD, vec_size,  1);
  b = std::make_shared<RawVectors<MatsT>>(MPI_COMM_WORLD, vec_size,  1);
  c = std::make_shared<RawVectors<MatsT>>(MPI_COMM_WORLD, vec_size,  1);

  // Test vecs
  //a 0 1 2 3 4
  //b 10 11 12 13 14
  
  auto print_vec = [](auto title, auto vec, auto size) {
      std::cout << "Vec " << title << "        : " << std::endl;
      for (int i = 0; i < size; i++){
         std::cout << vec->get(i,0) << "  ";
      }
      std::cout << std::endl;
  };

  // checking size
  std::cout << "Vec A size   : " << a->size() << std::endl;
  std::cout << "Vec A length : " << a->length() << std::endl;
  print_vec("A", a, vec_size);
  std::cout << "Vec B size   : " << b->size() << std::endl;
  std::cout << "Vec B length : " << b->length() << std::endl;
  print_vec("B", b, vec_size);

  ASSERT_EQ( a->size(), 1 );
  ASSERT_EQ( a->length(), vec_size );
  ASSERT_EQ( b->size(), 1 );
  ASSERT_EQ( b->length(), vec_size );

  // checking fill and copy

  RTMS::fill(a,MatsT(99));
  RTMS::copy(a,b);

  print_vec("A", a, vec_size);
  print_vec("B", b, vec_size);

  for (int i = 0; i < vec_size; i++){
    ASSERT_EQ( a->get(i,0), MatsT(99) );
    ASSERT_EQ( b->get(i,0), MatsT(99) );
    a->set(i,0, MatsT(i));
    b->set(i,0, MatsT(i+10));
  }

  print_vec("A", a, vec_size);
  print_vec("B", b, vec_size);
 
  // check add
  RTMS::fill(c,MatsT(0));
  RTMS::add(a,c, MatsT(0.5));
  RTMS::add(b,c, MatsT(0.5));
  print_vec("A+B", c, vec_size);
  ASSERT_EQ( c->get(0,0), MatsT(5) );
  ASSERT_EQ( c->get(1,0), MatsT(6) );
  ASSERT_EQ( c->get(2,0), MatsT(7) );
  ASSERT_EQ( c->get(3,0), MatsT(8) );
  ASSERT_EQ( c->get(4,0), MatsT(9) );
  
  // check dot
  MatsT result1= MatsT(0);
  MatsT result2= MatsT(0);
  RTMS::dot(a,b, result1);
  RTMS::dot(b,a, result2);
  ASSERT_EQ( result1, MatsT(130) );
  ASSERT_EQ( result2, MatsT(130) );
  std::cout << "a dot b      : " << result1 << std::endl;
  std::cout << "b dot a      : " << result2 << std::endl;

  // check scal
  RTMS::scal(c, MatsT(6));
  print_vec("6*0.5*(A+B)", c, vec_size);
  ASSERT_EQ( c->get(0,0), MatsT(30) );
  ASSERT_EQ( c->get(1,0), MatsT(36) );
  ASSERT_EQ( c->get(2,0), MatsT(42) );
  ASSERT_EQ( c->get(3,0), MatsT(48) );
  ASSERT_EQ( c->get(4,0), MatsT(54) );
  // check normalize
  MatsT result;
  RTMS::normalize(c, result);
  // c / np.sqrt(c@c) array([0.31311215, 0.37573457, 0.438357  , 0.50097943, 0.56360186])
  print_vec("normalized 3*(A+B)", c, vec_size);
  ASSERT_NEAR( c->get(0,0), MatsT(0.3131121455425747), TDCI_TEST_EPSILON );
  ASSERT_NEAR( c->get(1,0), MatsT(0.37573457465108967), TDCI_TEST_EPSILON );
  ASSERT_NEAR( c->get(2,0), MatsT(0.4383570037596046), TDCI_TEST_EPSILON );
  ASSERT_NEAR( c->get(3,0), MatsT(0.5009794328681195), TDCI_TEST_EPSILON );
  ASSERT_NEAR( c->get(4,0), MatsT(0.5636018619766345), TDCI_TEST_EPSILON );
  
} // Test End

void tdci_machinery_TEST_complex() {

  std::cout << "**** Test TDCI Infrastructure ****" << std::endl;
  
  const double TDCI_TEST_EPSILON = 1e-12;
  size_t mem     = 256e6; // Default 256 MB allocation
  size_t blkSize = 2048;  // Default 2KB block size

  typedef std::complex<double> MatsT;

  CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,mem,blkSize);
  
  size_t vec_size = 5;
  std::shared_ptr<SolverVectors<MatsT>> a,b,c;
  a = std::make_shared<RawVectors<MatsT>>(MPI_COMM_WORLD, vec_size,  1);
  b = std::make_shared<RawVectors<MatsT>>(MPI_COMM_WORLD, vec_size,  1);
  c = std::make_shared<RawVectors<MatsT>>(MPI_COMM_WORLD, vec_size,  1);

  // Test vecs
  //a 0 1 2 3 4
  //b 10 11 12 13 14
  
  auto print_vec = [](auto title, auto vec, auto size) {
      std::cout << "Vec " << title << "        : " << std::endl;
      for (int i = 0; i < size; i++){
         std::cout << vec->get(i,0) << "  ";
      }
      std::cout << std::endl;
  };
  
  auto near_complex = [TDCI_TEST_EPSILON](auto a, auto b) {
     ASSERT_NEAR(std::real(a), std::real(b), TDCI_TEST_EPSILON);
     ASSERT_NEAR(std::imag(a), std::imag(b), TDCI_TEST_EPSILON);
  };

  // checking size
  std::cout << "Vec A size   : " << a->size() << std::endl;
  std::cout << "Vec A length : " << a->length() << std::endl;
  std::cout << "Vec B size   : " << b->size() << std::endl;
  std::cout << "Vec B length : " << b->length() << std::endl;
  print_vec("B", b, vec_size);

  ASSERT_EQ( a->size(), 1 );
  ASSERT_EQ( a->length(), vec_size );
  ASSERT_EQ( b->size(), 1 );
  ASSERT_EQ( b->length(), vec_size );

  // checking fill and copy

  RTMS::fill(a,MatsT(99,0.5));
  RTMS::copy(a,b);

  print_vec("A", a, vec_size);
  print_vec("B", b, vec_size);

  for (int i = 0; i < vec_size; i++){
    near_complex( a->get(i,0), MatsT(99,0.5) );
    near_complex( b->get(i,0), MatsT(99,0.5) );
    a->set(i,0, MatsT(i,0.5));
    b->set(i,0, MatsT(i+10,0.5));
  }

  print_vec("A", a, vec_size);
  print_vec("B", b, vec_size);
 
  // check add
  RTMS::fill(c,MatsT(0,0));
  RTMS::add(a,c,MatsT(0.5,0));
  RTMS::add(b,c,MatsT(0.5,0));
  print_vec("A+B", c, vec_size);
  near_complex( c->get(0,0), MatsT(5,0.5) );
  near_complex( c->get(1,0), MatsT(6,0.5) );
  near_complex( c->get(2,0), MatsT(7,0.5) );
  near_complex( c->get(3,0), MatsT(8,0.5) );
  near_complex( c->get(4,0), MatsT(9,0.5) );
  
  // check dot
  MatsT result1= MatsT(0,0);
  MatsT result2= MatsT(0,0);
  RTMS::dot(a,b,result1);
  RTMS::dot(b,a,result2);
  near_complex( result1, MatsT(131.25,-25) );
  near_complex( result2, MatsT(131.25,25) );
  std::cout << "a dot b      : " << result1 << std::endl;
  std::cout << "b dot a      : " << result2 << std::endl;

  // check scal
  RTMS::scal(c, MatsT(6,0));
  print_vec("6*0.5*(A+B)", c, vec_size);
  near_complex( c->get(0,0), MatsT(30,3) );
  near_complex( c->get(1,0), MatsT(36,3) );
  near_complex( c->get(2,0), MatsT(42,3) );
  near_complex( c->get(3,0), MatsT(48,3) );
  near_complex( c->get(4,0), MatsT(54,3) );
  // check normalize
  MatsT result;
  RTMS::normalize(c,result);
  // c / np.sqrt(c@c) array([0.31311215, 0.37573457, 0.438357  , 0.50097943, 0.56360186])
  print_vec("normalized 3*(A+B)", c, vec_size);
  near_complex( c->get(0,0), MatsT(0.31234752377721214,0.031234752377721213) );
  near_complex( c->get(1,0), MatsT(0.37481702853265453,0.031234752377721213) );
  near_complex( c->get(2,0), MatsT(0.43728653328809697,0.031234752377721213) );
  near_complex( c->get(3,0), MatsT(0.4997560380435394, 0.031234752377721213) );
  near_complex( c->get(4,0), MatsT(0.5622255427989818, 0.031234752377721213) );
  
} // Test End


TEST(TDCI_MACHINERY, tdci_machinery_1) {
  tdci_machinery_TEST_real();
  tdci_machinery_TEST_complex();
}
