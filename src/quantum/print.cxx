#include <quantum/base.hpp>
#include <util/matout.hpp>
#include <physcon.hpp>

namespace ChronusQ {

  void QuantumBase::printMultipoles(std::ostream &out) {

    out << "\nMultipole Information:\n" << bannerTop << "\n\n";

    out << std::fixed << std::setprecision(10);


    out << "Electric Dipole Moment                            (Debye)\n";
    out << "X= " << std::setw(20) << elecDipole[0] / EBohrPerDebye()
      << " Y= " << std::setw(20) << elecDipole[1] / EBohrPerDebye()
      << " Z= " << std::setw(20) << elecDipole[2] / EBohrPerDebye()
        << "\n\n";

    out << "Electric Quadrupole Moment                        (Debye-Å)\n";
    out << "XX= " << std::setw(20) << elecQuadrupole[0][0] * AngPerBohr() / EBohrPerDebye()
      << " XY= " << std::setw(20) << elecQuadrupole[0][1] * AngPerBohr() / EBohrPerDebye()
      << " XZ= " << std::setw(20) << elecQuadrupole[0][2] * AngPerBohr() / EBohrPerDebye()
        << "\n";
    out << "YX= " << std::setw(20) << elecQuadrupole[1][0] * AngPerBohr() / EBohrPerDebye()
        << " YY= " << std::setw(20) << elecQuadrupole[1][1] * AngPerBohr() / EBohrPerDebye()
        << " YZ= " << std::setw(20) << elecQuadrupole[1][2] * AngPerBohr() / EBohrPerDebye()
        << "\n";
    out << "ZX= " << std::setw(20) << elecQuadrupole[2][0] * AngPerBohr() / EBohrPerDebye()
        << " ZY= " << std::setw(20) << elecQuadrupole[2][1] * AngPerBohr() / EBohrPerDebye()
        << " ZZ= " << std::setw(20) << elecQuadrupole[2][2] * AngPerBohr() / EBohrPerDebye()
        << "\n\n";

    out << "\n\n";
    out << "Electric Octupole Moment                          (Debye-Å²)\n";
    for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j) {
        for (int k=0; k<3; ++k) {
          out << static_cast<char>('X'+i) << static_cast<char>('X'+j) << static_cast<char>('X'+k) << "= "
              << std::setw(20) << elecOctupole[i][j][k] * AngPerBohr() * AngPerBohr() / EBohrPerDebye();
        }
        out << "\n";
      }
    }

    out << bannerEnd << std::endl;

  }

  




}
