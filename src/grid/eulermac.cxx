#include <grid/quadrature.hpp>

namespace ChronusQ {


/**
 *  Implemented as Eq 6. and 7 page 1001.
 *
 *  See
 *
 *  Christopher W. Murray, Nicholas C. Handy & Gregory J. Laming (1993):
 *    Quadrature schemes for integrals of density functional theory, Molecular Physics: 
 *    An International Journal at the Interface Between Chemistry and Physics, 78(4), 
 *    997-1014 
 *
 *  FIXME: get a full reference for this
 *  Johnson 169 Theor. Comp. Chem. Vol 2 (Jacobian included in the weights) 
 *
 *  Note that the Jacobian r^2 is not included (to added later on in the weights) and
 *  also the equeations are mapped from 0,inf respect to the paper according the following
 *  substitution 1 -> N+1, q -> i and w -> w * (N+1).
 *  The resulting equation including the Jacobian with m=2 are eq 24 and 25 in:
 *  Development, implementation and applications of efficient
 *  methodologies for density functional calculations
*/
  void EulerMac::generateQuadrature() {

    //std::cerr << "EulerMac 0" << std::endl;
    Quadrature::generateQuadrature();

    //std::cerr << "EulerMac 1" << std::endl;

    for(size_t i = 0; i < nPts; i++) {

      pts[i] = (i + 1.0) * (i + 1.0) / ((nPts - i) * (nPts - i));

      weights[i] = 2.0 * (i + 1.0) * (nPts + 1.0) / 
        ((nPts - i) * (nPts - i) * (nPts - i));

    }; // Loop over points

  };

}; // namespace ChronusQ
