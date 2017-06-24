#include <grid/quadrature.hpp>

namespace ChronusQ {

  void Lebedev::generateQuadrature() {
    if( nPts == 6  ) {
      #include "lebedev_6.cxx"
    } else if( nPts == 14 ) {
      #include "lebedev_14.cxx"
    } else if( nPts == 26 ) {
      #include "lebedev_26.cxx"
    } else if( nPts == 38 ) {
      #include "lebedev_38.cxx"
    } else if( nPts == 50 ) {
      #include "lebedev_50.cxx"
    } else if( nPts == 74 ) {
      #include "lebedev_74.cxx"
    } else if( nPts == 86 ) {
      #include "lebedev_86.cxx"
    } else if( nPts == 110 ) {
      #include "lebedev_110.cxx"
    } else if( nPts == 146 ) {
      #include "lebedev_146.cxx"
    } else if( nPts == 170 ) {
      #include "lebedev_170.cxx"
    } else if( nPts == 194 ) {
      #include "lebedev_194.cxx"
    } else if( nPts == 230 ) {
      #include "lebedev_230.cxx"
    } else if( nPts == 266 ) {
      #include "lebedev_266.cxx"
    } else if( nPts == 302 ) {
      #include "lebedev_302.cxx"
    } else if( nPts == 590 ) {
      #include "lebedev_590.cxx"
    } else if( nPts == 974 ) {
      #include "lebedev_974.cxx"
    } else {
      std::cout << "Invalid Lebedev Grid Specification" << std::endl;
    }
  };


};
