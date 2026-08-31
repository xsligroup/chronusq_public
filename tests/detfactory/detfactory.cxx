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
#include <detfactory.hpp>
#include <tensor.hpp>

using namespace ChronusQ;

void DetFactory_TEST() {

  std::cout << "**** Test DASCI Infrastructure ****" << std::endl;
  
  size_t mem     = 256e6; // Default 256 MB allocation
  size_t blkSize = 2048;  // Default 2KB block size

  CQMemManager::get().initialize(CQMemBackendType::PREALLOCATED,mem,blkSize);
  
#if 1
  std::cout << "/*                           " << std::endl; 
  std::cout << " * TEST SECTION 1: BASIC BITS OPERATIONS" << std::endl; 
  std::cout << " */                          " << std::endl; 
  uint8_t  * a1 = CQMemManager::get().malloc<uint8_t>(20);
  uint16_t * a2 = CQMemManager::get().malloc<uint16_t>(20);
  uint32_t * a3 = CQMemManager::get().malloc<uint32_t>(20);
  uint64_t * a4 = CQMemManager::get().malloc<uint64_t>(20);

  for (auto i = 0; i < 20; i++) {
    a1[i] = i;
    a2[i] = i;
    a3[i] = i;
    a4[i] = i;
  }
  
  std::cout << "**** Conversion from int to binary strings ****" << std::endl;
  
  std::cout << "uint8ToString :  " << +a1[15] << " to " << determinantToString(a1[15], 5) << std::endl;
  std::cout << "uint16ToString : " << a2[15]  << " to " << determinantToString(a2[15], 15) << std::endl;
  std::cout << "uint32ToString : " << a3[15]  << " to " << determinantToString(a3[15], 25) << std::endl;
  std::cout << "uint64ToString : " << a4[15]  << " to " << determinantToString(a4[15], 35) << std::endl;
  
  EXPECT_TRUE(determinantToString(a1[15],  5) == "11110");
  EXPECT_TRUE(determinantToString(a2[15], 15) == "111100000000000");
  EXPECT_TRUE(determinantToString(a3[15], 25) == "1111000000000000000000000");
  EXPECT_TRUE(determinantToString(a4[15], 35) == "11110000000000000000000000000000000");

  std::vector<uint8_t>  d1 = {a1[15], a1[5]};
  std::vector<uint16_t> d2 = {a2[15], a2[5]};
  std::vector<uint32_t> d3 = {a3[15], a3[5]};
  std::vector<uint64_t> d4 = {a4[15], a4[5]};
  std::vector<size_t> nO1 = {5, 5};
  std::vector<size_t> nO2 = {15, 15};
  std::vector<size_t> nO3 = {25, 25};
  std::vector<size_t> nO4 = {35, 35};
  
  std::cout << "uint8ToString :  (" << +a1[15] << "," <<  +a1[5] << ") to " 
            << determinantToString(d1, nO1) << std::endl;
  std::cout << "uint16ToString : (" << a2[15]  << "," <<  a2[5] << ") to " 
            << determinantToString(d2, nO2) << std::endl;
  std::cout << "uint32ToString : (" << a3[15]  << "," <<  a3[5] << ") to " 
            << determinantToString(d3, nO3) << std::endl;
  std::cout << "uint64ToString : (" << a4[15]  << "," <<  a4[5] << ") to " 
            << determinantToString(d4, nO4) << std::endl;
  
  EXPECT_TRUE(determinantToString(d1, nO1) == "11110 10100 ");
  EXPECT_TRUE(determinantToString(d2, nO2) == "111100000000000 101000000000000 ");
  EXPECT_TRUE(determinantToString(d3, nO3) == "1111000000000000000000000 1010000000000000000000000 ");
  EXPECT_TRUE(determinantToString(d4, nO4) == "11110000000000000000000000000000000 10100000000000000000000000000000000 ");

  CQMemManager::get().free(a1, a2, a3, a4);

  std::cout << "**** Conversion from binary strings to int ****" << std::endl;
  
  std::string s1 = "11010";

  std::cout << "stringToUInt8_t  : " << s1 << " to " << +stringToDeterminant<uint8_t >(s1) << std::endl;
  std::cout << "stringToUInt16_t : " << s1 << " to " << stringToDeterminant<uint16_t>(s1) << std::endl;
  std::cout << "stringToUInt32_t : " << s1 << " to " << stringToDeterminant<uint32_t>(s1) << std::endl;
  std::cout << "stringToUInt64_t : " << s1 << " to " << stringToDeterminant<uint64_t>(s1) << std::endl;
  
  EXPECT_TRUE(stringToDeterminant<uint8_t>(s1) == 11);
  EXPECT_TRUE(stringToDeterminant<uint16_t>(s1) == 11);
  EXPECT_TRUE(stringToDeterminant<uint32_t>(s1) == 11);
  EXPECT_TRUE(stringToDeterminant<uint64_t>(s1) == 11);

  std::string s2 = s1 + "10011";
  d1 = stringToDeterminant<uint8_t>(s2, nO1); 
  d2 = stringToDeterminant<uint16_t>(s2, nO1); 
  d3 = stringToDeterminant<uint32_t>(s2, nO1); 
  d4 = stringToDeterminant<uint64_t>(s2, nO1); 

  std::cout << "stringToUInt8_t  : " << s2 << " to (" << +d1[0] << "," << +d1[1] << ")" << std::endl;
  std::cout << "stringToUInt16_t : " << s2 << " to (" << d2[0] << "," << d2[1] << ")" << std::endl;
  std::cout << "stringToUInt32_t : " << s2 << " to (" << d3[0] << "," << d3[1] << ")" << std::endl;
  std::cout << "stringToUInt64_t : " << s2 << " to (" << d4[0] << "," << d4[1] << ")" << std::endl;
  
  EXPECT_TRUE(d1[0] == 11);
  EXPECT_TRUE(d1[1] == 25);
  EXPECT_TRUE(d2[0] == 11);
  EXPECT_TRUE(d2[1] == 25);
  EXPECT_TRUE(d3[0] == 11);
  EXPECT_TRUE(d3[1] == 25);
  EXPECT_TRUE(d4[0] == 11);
  EXPECT_TRUE(d4[1] == 25);

  auto occpos1 = occupiedInfo(d1[0], 5, 3);
  auto occpos2 = occupiedInfo(d2[0], 5, 3);
  auto occpos3 = occupiedInfo(d3[0], 5, 3);
  auto occpos4 = occupiedInfo(d4[0], 5, 3);
  auto virpos1 = virtualInfo(d1[0], 5, 3);
  auto virpos2 = virtualInfo(d2[0], 5, 3);
  auto virpos3 = virtualInfo(d3[0], 5, 3);
  auto virpos4 = virtualInfo(d4[0], 5, 3);
  
  std::cout << " string in uint8_t : " << s1 << std::endl 
            << "   occupied orbital: " << occpos1[0] << ", " << occpos1[1] << ", " << occpos1[2] << std::endl
            << "   virtual orbital:  " << virpos1[0] << ", " << virpos1[1] << std::endl; 
  std::cout << " string in uint16_t : " << s1 << std::endl 
            << "   occupied orbital: " << occpos2[0] << ", " << occpos2[1] << ", " << occpos2[2] << std::endl
            << "   virtual orbital:  " << virpos2[0] << ", " << virpos2[1] << std::endl; 
  std::cout << " string in uint32_t : " << s1 << std::endl 
            << "   occupied orbital: " << occpos3[0] << ", " << occpos3[1] << ", " << occpos3[2] << std::endl
            << "   virtual orbital:  " << virpos3[0] << ", " << virpos3[1] << std::endl; 
  std::cout << " string in uint64_t : " << s1 << std::endl 
            << "   occupied orbital: " << occpos4[0] << ", " << occpos4[1] << ", " << occpos4[2] << std::endl
            << "   virtual orbital:  " << virpos4[0] << ", " << virpos4[1] << std::endl; 
  
  EXPECT_TRUE(occpos1[0] == 0);
  EXPECT_TRUE(occpos1[1] == 1);
  EXPECT_TRUE(occpos1[2] == 3);
  EXPECT_TRUE(virpos1[0] == 2);
  EXPECT_TRUE(virpos1[1] == 4);

  EXPECT_TRUE(occpos2[0] == 0);
  EXPECT_TRUE(occpos2[1] == 1);
  EXPECT_TRUE(occpos2[2] == 3);
  EXPECT_TRUE(virpos2[0] == 2);
  EXPECT_TRUE(virpos2[1] == 4);

  EXPECT_TRUE(occpos3[0] == 0);
  EXPECT_TRUE(occpos3[1] == 1);
  EXPECT_TRUE(occpos3[2] == 3);
  EXPECT_TRUE(virpos3[0] == 2);
  EXPECT_TRUE(virpos3[1] == 4);

  EXPECT_TRUE(occpos4[0] == 0);
  EXPECT_TRUE(occpos4[1] == 1);
  EXPECT_TRUE(occpos4[2] == 3);
  EXPECT_TRUE(virpos4[0] == 2);
  EXPECT_TRUE(virpos4[1] == 4);

  std::cout << "* Bit Mask Test for uint8_t " << std::endl;
  std::cout << "  - Type A, bitAt" << std::endl;
  for (auto i = 0ul; i < 8; i++) {
    std::cout << "    " << i << ": " << determinantToString(BITMASKS<uint8_t>::bitAt[i], 8) << std::endl;
    std::string ref = "00000000";
    ref[i] = '1';
    EXPECT_TRUE(determinantToString(BITMASKS<uint8_t>::bitAt[i], 8) == ref);
  }
  std::cout << "  - Type B, bitsPrevTo" << std::endl;
  for (auto i = 0ul; i < 9; i++) {
    std::cout << "    " << i << ": " << determinantToString(BITMASKS<uint8_t>::bitsPrevTo[i], 8) << std::endl;
    std::string ref = "00000000";
    for (auto j = 0ul; j < i; j++){
      ref[j] = '1';
    }
    EXPECT_TRUE(determinantToString(BITMASKS<uint8_t>::bitsPrevTo[i], 8) == ref);
  }
  std::cout << "  - Type C, bitsAfter" << std::endl;
  for (auto i = 0ul; i < 9; i++) {
    std::cout << "    " << i << ": " << determinantToString(BITMASKS<uint8_t>::bitsAfter[i], 8) << std::endl;
    std::string ref = "00000000";
    for (auto j = i; j < ref.size(); j++){
      ref[j] = '1';
    }
    EXPECT_TRUE(determinantToString(BITMASKS<uint8_t>::bitsAfter[i], 8) == ref);
  }

  std::cout << "* Test for bit manipulation" << std::endl;
  uint8_t dtest = stringToDeterminant<uint8_t>("00011010");
  std::cout << "  - Test for flip bit of string 00011010" << std::endl;
  for (auto i = 0ul; i < 8; i++) {
    std::cout << "    flip bit " << i << ": " << determinantToString(flipBit(dtest, i), 8) << std::endl;
    std::string ref = "00011010";
    ref[i] = (ref[i] == '1') ? '0' : '1';
    EXPECT_TRUE(determinantToString(flipBit(dtest, i), 8) == ref);
  }
  
  std::cout << "  - Test for check bit of string 00011010" << std::endl;
  for (auto i = 0ul; i < 8; i++) {
    std::cout << "    check bit " << i << "? " << checkBit(dtest, i) << std::endl;
    std::string ref = "00011010";
    EXPECT_TRUE(checkBit(dtest,i) == (ref[i] == '1'));
  }
  
  std::cout << "  - Test for negativeSign1eExcitation of string 00011010" << std::endl;
  std::cout << "    Exctiation 3 -> 0, negativeSign ? " << +negativeSign1eExcitation(dtest, 0, 3) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,0,3) == 0);
  std::cout << "    Exctiation 3 -> 7, negativeSign ? " << +negativeSign1eExcitation(dtest, 7, 3) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,7,3) == 0);
  std::cout << "    Exctiation 3 -> 5, negativeSign ? " << +negativeSign1eExcitation(dtest, 5, 3) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,5,3) == 1);
  std::cout << "    Exctiation 4 -> 1, negativeSign ? " << +negativeSign1eExcitation(dtest, 1, 4) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,1,4) == 1);
  std::cout << "    Exctiation 4 -> 4, negativeSign ? " << +negativeSign1eExcitation(dtest, 4, 4) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,4,4) == 0);
  std::cout << "    Exctiation 4 -> 5, negativeSign ? " << +negativeSign1eExcitation(dtest, 4, 5) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,4,5) == 0);
  std::cout << "    Exctiation 4 -> 7, negativeSign ? " << +negativeSign1eExcitation(dtest, 4, 7) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,4,7) == 1);
  
  std::cout << "  - Test for negativeSign1eExcitation of string 00011010 to 11010101" << std::endl;
  uint8_t dtest2 = stringToDeterminant<uint8_t>("11010101");
  std::cout << "    Exctiation 3 -> 2, negativeSign ? " << +negativeSign1eExcitation(dtest, dtest2, 3, 2) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,dtest2,3,2) == 0);
  std::cout << "    Exctiation 3 -> 4, negativeSign ? " << +negativeSign1eExcitation(dtest, dtest2, 3, 4) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,dtest2,3,4) == 1);
  std::cout << "    Exctiation 3 -> 6, negativeSign ? " << +negativeSign1eExcitation(dtest, dtest2, 3, 6) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,dtest2,3,6) == 0);
  std::cout << "    Exctiation 4 -> 2, negativeSign ? " << +negativeSign1eExcitation(dtest, dtest2, 4, 2) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,dtest2,4,2) == 1);
  std::cout << "    Exctiation 4 -> 4, negativeSign ? " << +negativeSign1eExcitation(dtest, dtest2, 4, 4) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,dtest2,4,4) == 0);
  std::cout << "    Exctiation 4 -> 6, negativeSign ? " << +negativeSign1eExcitation(dtest, dtest2, 4, 6) << std::endl;
  EXPECT_TRUE(negativeSign1eExcitation(dtest,dtest2,4,6) == 1);
 
  constexpr static uint16_t comb_10_4 = combination<uint16_t, 10,4>();
  constexpr static uint16_t comb_8_5  = combination<uint16_t, 8, 5>();
  std::cout << "combination(10, 4) = " << comb_10_4 << std::endl;
  EXPECT_TRUE(comb_10_4 == 210);
  std::cout << "combination(8, 5) = " << comb_8_5 << std::endl;
  EXPECT_TRUE(comb_8_5 == 56);
  constexpr static uint16_t bw_5_4  = bitWeight<uint16_t, 5, 4>();
  constexpr static uint16_t bw_5_5  = bitWeight<uint16_t, 5, 5>();
  constexpr static uint16_t bw_5_6  = bitWeight<uint16_t, 5, 6>();
  constexpr static uint16_t bw_5_7  = bitWeight<uint16_t, 5, 7>();
  constexpr static uint16_t bw_5_8  = bitWeight<uint16_t, 5, 8>();
  std::cout << "bitWeight(5, 4) = "  << bw_5_4 << std::endl;
  EXPECT_TRUE(bw_5_4 == 0);
  std::cout << "bitWeight(5, 5) = "  << bw_5_5 << std::endl;
  EXPECT_TRUE(bw_5_5 == 0);
  std::cout << "bitWeight(5, 6) = "  << bw_5_6 << std::endl;
  EXPECT_TRUE(bw_5_6 == 1);
  std::cout << "bitWeight(5, 7) = "  << bw_5_7 << std::endl;
  EXPECT_TRUE(bw_5_7 == 7);
  std::cout << "bitWeight(5, 8) = "  << bw_5_8 << std::endl;
  EXPECT_TRUE(bw_5_8 == 28);
  
  std::cout << "bitWeight(5, 4) = "  << BITADDRESSING<uint16_t>::array[5][4] << std::endl;
  EXPECT_TRUE(BITADDRESSING<uint16_t>::array[5][4] == 0);
  std::cout << "bitWeight(5, 5) = "  << BITADDRESSING<uint16_t>::array[5][5] << std::endl;
  EXPECT_TRUE(BITADDRESSING<uint16_t>::array[5][5] == 0);
  std::cout << "bitWeight(5, 6) = "  << BITADDRESSING<uint16_t>::array[5][6] << std::endl;
  EXPECT_TRUE(BITADDRESSING<uint16_t>::array[5][6] == 1);
  std::cout << "bitWeight(5, 7) = "  << BITADDRESSING<uint16_t>::array[5][7] << std::endl;
  EXPECT_TRUE(BITADDRESSING<uint16_t>::array[5][7] == 7);
  std::cout << "bitWeight(5, 8) = "  << BITADDRESSING<uint16_t>::array[5][8] << std::endl;
  EXPECT_TRUE(BITADDRESSING<uint16_t>::array[5][8] == 28);
  
  std::cout << " * Test TensorLooper " << std::endl;
  auto TL_ptr = constructTensorLooper({4,3,2}, {1,0,0}, {0,1,0});
  auto const & TL = *TL_ptr;
  auto i = 0, iLim = 5;
  auto j = 0, jLim = 4;
  auto k = 0, kLim = 3;
  auto count = 0;
  
  for (TL_ptr->setIndex(); not TL_ptr->isEnd(); TL_ptr->increment()) {
    std::cout << " index = " << TL.index() 
              << " Address = " << TL.address() 
              << " auxAddress = " << TL.auxAddress()
              << " indices = (" << TL[0] << ", " << TL[1] << ", " << TL[2] << ")" 
              << std::endl;
    EXPECT_TRUE(TL.index() == count);
    EXPECT_TRUE(TL.address() == i);
    EXPECT_TRUE(TL.auxAddress() == j);
    EXPECT_TRUE(TL[0] == i);
    EXPECT_TRUE(TL[1] == j);
    EXPECT_TRUE(TL[2] == k);
    count++;
    i++;
    if ((i+1) % iLim == 0){
      i = 0;
      j++;
    }
    if ((j+1) % jLim == 0){
      j = 0;
      k++;
    }
  }
  std::cout << "/*                           " << std::endl; 
  std::cout << " * END OF TEST SECTION 1: BASIC BITS OPERATIONS" << std::endl; 
  std::cout << " */                          " << std::endl; 
#endif  

#if 1
  std::cout << "/*                           " << std::endl; 
  std::cout << " * TEST SECTION 2: DETERMINANTS MANAGEMENT AND GENERATION " << std::endl; 
  std::cout << " */                          " << std::endl; 
  
  size_t testNE = 2;
  size_t testNO = 4;
  count = 0;
  DeterminantGroup testG(testNE, testNO);
  std::string testCAS = "CAS(" + std::to_string(testNE) + ", " + std::to_string(testNO) + ")";
  std::vector<std::string> ref_strings = {"1100", "1010", "0110", "1001", "0101", "0011"};
  std::cout << "* Iterate through determiants in "<< testCAS << ": " << std::endl;
  auto testGGen = testG.template generator<uint8_t>();
  testGGen.visitDeterminants(0ul, testG.nDeterminants(), 
      [&](size_t addr, uint8_t det) {
        std::cout << "  - string " << addr << ": " << determinantToString(det, testNO) 
                  << ", addr = "   << +detToLexicographicAddr(det, testNE, testNO) 
                  << ", string = " << determinantToString(lexicographicAddrToBitString(uint8_t(addr), testNE, testNO), testNO)
                  << std::endl;
        EXPECT_TRUE(addr == count);
        EXPECT_TRUE(detToLexicographicAddr(det, testNE, testNO) == count);
        EXPECT_TRUE(determinantToString(det, testNO) == ref_strings[count]);
        EXPECT_TRUE(determinantToString(lexicographicAddrToBitString(uint8_t(addr), testNE, testNO), testNO) == ref_strings[count]);
        count++;
      }
  );

  std::cout << testG << std::endl;

  std::cout << "* Iterate through determiants in " << testCAS << "-" << testCAS << ": " << std::endl;
  std::vector<size_t> testNOs = {testNO, testNO};
  std::vector<size_t> testNEs = {testNE, testNE};
  count = 0;
  FullDeterminantCategory testCat(testNEs, testNOs);
  auto testCatGen = testCat.template generator<uint8_t>();
  auto testCatAddresser = testCatGen.addresser();
  testCatGen.visitDeterminants(0ul, testCat.nDeterminants(), 
      [&](size_t addr, std::vector<uint8_t>& dets) {
        std::cout << " - string " << addr << ": " <<  testCatAddresser.detsToString(dets)
                  << ", addr = " << testCatAddresser.detsToAddress(dets);
        testCatAddresser.addressToBitStrings(addr, dets);
        std::cout << ", string = " << testCatAddresser.detsToString(dets) 
                 << std::endl;
        int i = count / ref_strings.size();
        int j = count % ref_strings.size();

        std::string ref_det = ref_strings[j] + ref_strings[i];
        EXPECT_TRUE(addr == count);
        EXPECT_TRUE(testCatAddresser.detsToAddress(dets) == count);

        EXPECT_TRUE(testCatAddresser.detsToString(dets) == ref_det);
        count++;
      }
  );
  
  std::cout << testCat << std::endl;

  std::vector<ActiveSpaceParameters> actSpaces = {ActiveSpaceParameters({0, 8, 0, 8}),
                                                  ActiveSpaceParameters({8, 8, 0, 8}),
                                                  ActiveSpaceParameters({16, 8, 0, 8}) };
  
  std::cout << "Active Space:" << std::endl;
  std::cout << actSpaces << std::endl;

  CategoricalSpace detsSpace(actSpaces);
  
  std::vector<size_t> contTestNOs = {8,8,8};
  std::vector<size_t> contTestNEs = {8,4,0};
  
  auto refCat = std::make_shared<FullDeterminantCategory>(contTestNEs, contTestNOs);
  
  detsSpace.addCategory(refCat);
  detsSpace.expandCategoryWithInterGroupExcitation();
  detsSpace.output(std::cout, "Testing RAS Category Generation");
  std::cout << detsSpace << std::endl;

  detsSpace.expandCategoryWithInterGroupExcitation(2);
  EXPECT_TRUE(detsSpace.nDeterminants() == 70);
  
  // detsSpace.output(std::cout, "Testing RAS Category Generation");
  std::cout << "/*                           " << std::endl; 
  std::cout << " * END OF TEST SECTION 2: DETERMINANTS MANAGEMENT AND GENERATION " << std::endl; 
  std::cout << " */                          " << std::endl; 
#endif  
  
#if 1
  std::cout << "/*                           " << std::endl; 
  std::cout << " * TEST SECTION 3: EXCITATIONLIST AND ITERATION" << std::endl; 
  std::cout << " */                          " << std::endl; 

  IntraSpaceFullCD1eExList<uint8_t, uint8_t> test1eIntraSpaceExL(testG);
  std::cout << "computeExctiationList" << std::endl;
  test1eIntraSpaceExL.computeExcitationList();
  
  std::cout << "* Iterate through test intra excitation in " << testCAS << std::endl;
  auto test1eIntraSpaceExLGen = test1eIntraSpaceExL.generator();
  {
  size_t count = 0ul;
  const auto& KEx = test1eIntraSpaceExLGen->braExAddress();
  std::vector<int> ref_KEx = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5};

  std::vector<int> ref_p = {0, 1, 0, 0, 1, 1, 0, 2, 0, 0, 2, 2, 1, 2, 1, 1, 2, 2, 0, 3, 0, 0, 3, 3, 1, 3, 1, 1, 3, 3, 2, 3, 2, 2, 3, 3};



  std::vector<int> ref_q = {0, 1, 2, 3, 2, 3, 0, 2, 1, 3, 1, 3, 1, 2, 0, 3, 0, 3, 0, 3, 1, 2, 1, 2, 1, 3, 0, 2, 0, 2, 2, 3, 0, 1, 0, 1};

  std::vector<int> ref_LEx = {0, 0, 2, 4, 1, 3, 1, 1, 2, 5, 0, 3, 2, 2, 1, 5, 0, 4, 3, 3, 4, 5, 0, 1, 4, 4, 3, 5, 0, 2, 5, 5, 3, 4, 1, 2};

  std::vector<int> ref_sign = {0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1};

  test1eIntraSpaceExLGen->visitAllExcitations(
      [&] (size_t K, const auto& excitations) {
        for (const auto& [p, q, LEx, sign] : excitations) {
          std::cout << " - Excitation " << count
            << ": " << KEx    
            << " " << p     
            << " " << q     
            << " " << LEx     
            << " " << sign  
            << std::endl;
            EXPECT_TRUE(KEx == ref_KEx[count]);
            EXPECT_TRUE(p == ref_p[count]);
            EXPECT_TRUE(q == ref_q[count]);
            EXPECT_TRUE(LEx == ref_LEx[count]);
            EXPECT_TRUE(sign == ref_sign[count]);
            ++count;
        }
      }
  );
  }

  size_t testNE2 = 1;
  size_t testNO2 = 6;
  DeterminantGroup testG2(testNE2, testNO2);
  std::string testCAS2 = "CAS(" + std::to_string(testNE2) + ", " + std::to_string(testNO2) + ")";
  
  InterSpaceFullCD1eExList<uint8_t, uint8_t> test1eInterSpaceExL(testG, testG2, false);
  std::cout << "* Test on inter space excitation list, with NDetK = " 
            << test1eInterSpaceExL.braCategory()->nDeterminants() 
            << " NDetL = " 
            << test1eInterSpaceExL.ketCategory()->nDeterminants() 
            << std::endl; 
  test1eInterSpaceExL.computeExcitationList();
  
  std::cout << "* Iterate through test inter excitation in " << testCAS << "-" << testCAS2 << std::endl;
  auto test1eInterSpaceExLGen = test1eInterSpaceExL.generator();
  {
  size_t count = 0ul; 
  std::vector<int> ref_KEx = {0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,1 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,2 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,3 ,4 ,4 ,4 ,4 ,4 ,4 ,4 ,4 ,4 ,4 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,5 ,6 ,6 ,6 ,6 ,6 ,6 ,6 ,6 ,6 ,6 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,7 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,8 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,9 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,10 ,11 ,11 ,11 ,11 ,11 ,11 ,11 ,11 ,11 ,11 ,12 ,12 ,12 ,12 ,12 ,12 ,12 ,12 ,12 ,12 ,13 ,13 ,13 ,13 ,13 ,13 ,13 ,13 ,13 ,13 ,14 ,14 ,14 ,14 ,14 ,14 ,14 ,14 ,14 ,14 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,15 ,16 ,16 ,16 ,16 ,16 ,16 ,16 ,16 ,16 ,16 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,17 ,18 ,18 ,18 ,18 ,18 ,18 ,18 ,18 ,18 ,18 ,19 ,19 ,19 ,19 ,19 ,19 ,19 ,19 ,19 ,19 ,20 ,20 ,20 ,20 ,20 ,20 ,20 ,20 ,20 ,20 ,21 ,21 ,21 ,21 ,21 ,21 ,21 ,21 ,21 ,21 ,22 ,22 ,22 ,22 ,22 ,22 ,22 ,22 ,22 ,22 ,23 ,23 ,23 ,23 ,23 ,23 ,23 ,23 ,23 ,23 ,24 ,24 ,24 ,24 ,24 ,24 ,24 ,24 ,24 ,24 ,25 ,25 ,25 ,25 ,25 ,25 ,25 ,25 ,25 ,25 ,26 ,26 ,26 ,26 ,26 ,26 ,26 ,26 ,26 ,26 ,27 ,27 ,27 ,27 ,27 ,27 ,27 ,27 ,27 ,27 ,28 ,28 ,28 ,28 ,28 ,28 ,28 ,28 ,28 ,28 ,29 ,29 ,29 ,29 ,29 ,29 ,29 ,29 ,29 ,29 ,30 ,30 ,30 ,30 ,30 ,30 ,30 ,30 ,30 ,30 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,31 ,32 ,32 ,32 ,32 ,32 ,32 ,32 ,32 ,32 ,32 ,33 ,33 ,33 ,33 ,33 ,33 ,33 ,33 ,33 ,33 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,34 ,35 ,35 ,35 ,35 ,35 ,35 ,35 ,35 ,35 ,35 };

  std::vector<int> ref_p = {0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,2 ,2 ,2 ,2 ,2 ,1 ,1 ,1 ,1 ,1 ,2 ,2 ,2 ,2 ,2 ,0 ,0 ,0 ,0 ,0 ,3 ,3 ,3 ,3 ,3 ,1 ,1 ,1 ,1 ,1 ,3 ,3 ,3 ,3 ,3 ,2 ,2 ,2 ,2 ,2 ,3 ,3 ,3 ,3 ,3 ,0 ,0 ,0 ,0 ,0 ,1 ,1 ,1 ,1 ,1 ,0 ,0 ,0 ,0 ,0 ,2 ,2 ,2 ,2 ,2 ,1 ,1 ,1 ,1 ,1 ,2 ,2 ,2 ,2 ,2 ,0 ,0 ,0 ,0 ,0 ,3 ,3 ,3 ,3 ,3 , 1 , 1 , 1 , 1 , 1 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 3 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 3 , 3 , 3 , 3 , 3 , 1 , 1 , 1 , 1 , 1 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 3 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 3 , 3 , 3 , 3 , 3 , 1 , 1 , 1 , 1 , 1 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 3 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 3 , 3 , 3 , 3 , 3 , 1 , 1 , 1 , 1 , 1 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 3 , 0 , 0 , 0 , 0 , 0 , 1 , 1 , 1 , 1 , 1 , 0 , 0 , 0 , 0 , 0 , 2 , 2 , 2 , 2 , 2 , 1 , 1 , 1 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 0 , 0 , 0 , 0 , 0 , 3 , 3 , 3 , 3 , 3 , 1 , 1 , 1 , 1 , 1 , 3 , 3 , 3 , 3 , 3 , 2 , 2 , 2 , 2 , 2 , 3 , 3 , 3 , 3 , 3 };

  std::vector<int> ref_q = {1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,1 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,2 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,3 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,4 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,5 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 ,0 ,1 ,2 ,3 ,4 };

  std::vector<int> ref_LEx = {1 ,5 ,13 ,25 ,41 ,0 ,4 ,12 ,24 ,40 ,2 ,6 ,14 ,26 ,42 ,0 ,4 ,12 ,24 ,40 ,2 ,6 ,14 ,26 ,42 ,1 ,5 ,13 ,25 ,41 ,3 ,7 ,15 ,27 ,43 ,0 ,4 ,12 ,24 ,40 ,3 ,7 ,15 ,27 ,43 ,1 ,5 ,13 ,25 ,41 ,3 ,7 ,15 ,27 ,43 ,2 ,6 ,14 ,26 ,42 ,1 ,9 ,17 ,29 ,45 ,0 ,8 ,16 ,28 ,44 ,2 ,10 ,18 ,30 ,46 ,0 ,8 ,16 ,28 ,44 ,2 ,10 ,18 ,30 ,46 ,1 ,9 ,17 ,29 ,45 ,3 ,11 ,19 ,31 ,47 ,0 ,8 ,16 ,28 ,44 ,3 ,11 ,19 ,31 ,47 ,1 ,9 ,17 ,29 ,45 ,3 ,11 ,19 ,31 ,47 ,2 ,10 ,18 ,30 ,46 ,5 ,9 ,21 ,33 ,49 ,4 ,8 ,20 ,32 ,48 ,6 ,10 ,22 ,34 ,50 ,4 ,8 ,20 ,32 ,48 ,6 ,10 ,22 ,34 ,50 ,5 ,9 ,21 ,33 ,49 ,7 ,11 ,23 ,35 ,51 ,4 ,8 ,20 ,32 ,48 ,7 ,11 ,23 ,35 ,51 ,5 ,9 ,21 ,33 ,49 ,7 ,11 ,23 ,35 ,51 ,6 ,10 ,22 ,34 ,50 ,13 ,17 ,21 ,37 ,53 ,12 ,16 ,20 ,36 ,52 ,14 ,18 ,22 ,38 ,54 ,12 ,16 ,20 ,36 ,52 ,14 ,18 ,22 ,38 ,54 ,13 ,17 ,21 ,37 ,53 ,15 ,19 ,23 ,39 ,55 ,12 ,16 ,20 ,36 ,52 ,15 ,19 ,23 ,39 ,55 ,13 ,17 ,21 ,37 ,53 ,15 ,19 ,23 ,39 ,55 ,14 ,18 ,22 ,38 ,54 ,25 ,29 ,33 ,37 ,57 ,24 ,28 ,32 ,36 ,56 ,26 ,30 ,34 ,38 ,58 ,24 ,28 ,32 ,36 ,56 ,26 ,30 ,34 ,38 ,58 ,25 ,29 ,33 ,37 ,57 ,27 ,31 ,35 ,39 ,59 ,24 ,28 ,32 ,36 ,56 ,27 ,31 ,35 ,39 ,59 ,25 ,29 ,33 ,37 ,57 ,27 ,31 ,35 ,39 ,59 ,26 ,30 ,34 ,38 ,58 ,41 ,45 ,49 ,53 ,57 ,40 ,44 ,48 ,52 ,56 ,42 ,46 ,50 ,54 ,58 ,40 ,44 ,48 ,52 ,56 ,42 ,46 ,50 ,54 ,58 ,41 ,45 ,49 ,53 ,57 ,43 ,47 ,51 ,55 ,59 ,40 ,44 ,48 ,52 ,56 ,43 ,47 ,51 ,55 ,59 ,41 ,45 ,49 ,53 ,57 ,43 ,47 ,51 ,55 ,59 ,42 ,46 ,50 ,54 ,58 };

  std::vector<int> ref_sign = {0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0};


  const auto& KEx = test1eInterSpaceExLGen->braExAddress();
  test1eInterSpaceExLGen->visitAllExcitations(
      [&] (size_t K, const auto& excitations) {
        for (const auto& [p, q, LEx, sign] : excitations) {
          std::cout << " - Excitation " << count
            << ": " << KEx    
            << " " << p     
            << " " << q     
            << " " << LEx     
            << " " << sign  
            << std::endl;
            EXPECT_TRUE(KEx == ref_KEx[count]);
            EXPECT_TRUE(p == ref_p[count]);
            EXPECT_TRUE(q == ref_q[count]);
            EXPECT_TRUE(LEx == ref_LEx[count]);
            EXPECT_TRUE(sign == ref_sign[count]);
            ++count;
        }
      }
  );
  }


  InterSpaceFullCD1eExList<uint8_t, uint8_t> test1eInterSpaceExL2(testG2, testG, true);
  std::cout << "* Test on inter space excitation list, with NDetK = " 
            << test1eInterSpaceExL2.braCategory()->nDeterminants() 
            << " NDetL = " 
            << test1eInterSpaceExL2.ketCategory()->nDeterminants() 
            << std::endl; 
  test1eInterSpaceExL2.computeExcitationList();

  std::cout << "* Iterate through test inter excitation in " << testCAS2 << "-" << testCAS << std::endl;
  auto test1eInterSpaceExLGen2 = test1eInterSpaceExL2.generator();
  {
  size_t count = 0ul; 
  std::vector<int> ref_KEx = { 0, 0, 1, 1,  2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,  7,  8,  8,  9,  9,  10,  10,  11,  11,  12,  12,  13,  13,  14,  14,  15,  15,  16,  16,  17,  17,  18,  18,  19,  19,  20,  20,  21,  21,  22,  22,  23,  23,  24,  24,  25,  25,  26,  26,  27,  27,  28,  28,  29,  29,  30,  30,  31,  31,  32,  32,  33,  33,  34,  34,  35,  35};
  std::vector<int> ref_p = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5};
  std::vector<int> ref_q = {2,3,2,3,2,3,2,3,2,3,2,3,1,3,1,3,1,3,1,3,1,3,1,3,0,3,0,3,0,3,0,3,0,3,0,3,1,2,1,2,1,2,1,2,1,2,1,2,0,2,0,2,0,2,0,2,0,2,0,2,0,1,0,1,0,1,0,1,0,1,0,1};
  std::vector<int> ref_LEx = { 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 2, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 0, 3, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3};
  std::vector<int> ref_sign = { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0};


  const auto& KEx = test1eInterSpaceExLGen2->braExAddress();
  test1eInterSpaceExLGen2->visitAllExcitations(
      [&] (size_t K, const auto& excitations) {
        for (const auto& [p, q, LEx, sign] : excitations) {
          std::cout << " - Excitation " << count
            << ": " << KEx    
            << " " << p     
            << " " << q     
            << " " << LEx     
            << " " << sign  
            << std::endl;
            EXPECT_TRUE(KEx == ref_KEx[count]);
            EXPECT_TRUE(p == ref_p[count]);
            EXPECT_TRUE(q == ref_q[count]);
            EXPECT_TRUE(LEx == ref_LEx[count]);
            EXPECT_TRUE(sign == ref_sign[count]);
            ++count;
        }
      }
  );
  }

  std::cout << "/*                           " << std::endl; 
  std::cout << " * END OF TEST SECTION 3: EXCITATIONLIST AND ITERATION" << std::endl; 
  std::cout << " */                          " << std::endl; 
#endif

#if 1
  std::cout << "/*                           " << std::endl; 
  std::cout << " * TEST SECTION 4: INTERMEDIATES AND BUILD CONTRACTION" << std::endl; 
  std::cout << " */                          " << std::endl; 
  contTestNOs = {3,2,5,2,4,3};
  contTestNEs = {1,1,2,1,1,1};
  FullDeterminantCategory contCatK(contTestNEs, contTestNOs);
  DeterminantGroup contTestG1(2, 5);
  IntraSpaceFullCD1eExList<uint8_t, uint8_t> test1eContIntraSpaceExL(contTestG1);
  test1eContIntraSpaceExL.computeExcitationList();
  
  std::cout << " * " << test1eContIntraSpaceExL << std::endl;

  DeterminantGroup contTestG2(1, 4);
  InterSpaceFullCD1eExList<uint8_t, uint8_t> test1eContInterSpaceExL(contTestG1, contTestG2, false);
  contTestNEs[2]--;
  contTestNEs[4]++;
  FullDeterminantCategory contCatL(contTestNEs, contTestNOs);
  test1eContInterSpaceExL.computeExcitationList();
  
  DoubleFullCD1eExListGenerator double1eExListX(test1eContIntraSpaceExL, test1eContInterSpaceExL, {2, 2, 2, 4});
  {
    size_t count = 0;
    double1eExListX.visitAllExcitations({1, 2}, {1, 2}, 
        [&](const auto& JExIndex, 
           const auto& qpContractions, 
           const auto& rsContractions) {
          std::cout << "In Builder now, JExIndex = " << JExIndex << std::endl; 
          EXPECT_TRUE(JExIndex == count);
          ++count;
        }
    );
  } 
  
  std::cout << "/*                           " << std::endl; 
  std::cout << " * END OF TEST SECTION 4: INTERMEDIATES AND BUILD CONTRACTION" << std::endl; 
  std::cout << " */                          " << std::endl; 
#endif

#if 1
  std::cout << "/*                           " << std::endl; 
  std::cout << " * TEST SECTION 5: DETFACTORY" << std::endl; 
  std::cout << " */                          " << std::endl; 
  std::cout << "*---Test for DetFactory" << std::endl;
  MPI_Comm comm(MPI_COMM_WORLD);
  DeterminantFactory detsFac(comm, 12, actSpaces);
  
  auto dS = std::make_shared<CategoricalSpace>(detsSpace);
  detsFac.setKetCategoricalSpace(dS);
  detsFac.setBraCategoricalSpace(dS);
  detsFac.generateComputingGraph();
  detsFac.output(std::cout, "Test Det Factory");  
  
  std::cout << "*--- Test for interaction parsing "  << std::endl;
     
  for( const auto& term: detsFac.oneEExTerms() ) {
    std::cout << " - " << term << std::endl; 
    std::vector<std::string> oneETerms;
    std::vector<size_t> oneESpan;
    parseTermSpan(term, oneETerms, oneESpan, "+RI[]"); 
    
    std::cout << "    * " << oneETerms[0] << ", span: ";
    for (const auto & i: oneESpan) 
      std::cout << i << " ";
    std::cout << std::endl;
    for (auto i = 0ul; i < oneETerms.size(); ++i) {
      std::cout << "   * " << oneETerms[i] << std::endl;
    };  
  }

  for( const auto term: detsFac.twoEExTerms() ) {
    std::cout << " - " << term << std::endl; 
    std::vector<std::string> twoETerms;
    std::vector<size_t> twoESpan;
    parseTermSpan(term, twoETerms, twoESpan, "-"); 
    
    std::cout << "    * " <<  twoETerms[0] << ", span: ";
    for (const auto & i: twoESpan) 
      std::cout << i << " ";
    std::cout << std::endl;
    for (auto i = 0ul; i < twoETerms.size(); ++i) {
      std::cout << "   * " << twoETerms[i] << std::endl;
    };  
  }
  std::cout << "/*                           " << std::endl; 
  std::cout << " * END OF TEST SECTION 5: DETFACTORY" << std::endl; 
  std::cout << " */                          " << std::endl; 
#endif

} // Test End


TEST(DETFACTORY, DETFACTORY_TEST_1) {
  DetFactory_TEST();
}
