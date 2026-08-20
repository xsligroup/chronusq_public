#include <quantum/base.hpp>
#include <util/matout.hpp>
#include <physcon.hpp>

namespace ChronusQ {

  void QuantumBase::printMultipoles(std::ostream &out) {

    out << "\nMultipole Information:\n" << bannerTop << "\n\n";

    out << std::fixed << std::setprecision(10);

    // Use constants from physcon.hpp: AngPerBohr and EBohrPerDebye
    const double BOHR_TO_ANGSTROM = AngPerBohr();
    const double CHARGE_AU_TO_DEBYE = 1.0 / EBohrPerDebye();

    out << "Electric Dipole Moment                            (Debye)\n";
    out << "X= " << std::setw(20) << elecDipole[0] * BOHR_TO_ANGSTROM * CHARGE_AU_TO_DEBYE
      << " Y= " << std::setw(20) << elecDipole[1] * BOHR_TO_ANGSTROM * CHARGE_AU_TO_DEBYE
      << " Z= " << std::setw(20) << elecDipole[2] * BOHR_TO_ANGSTROM * CHARGE_AU_TO_DEBYE
        << "\n\n";

    out << "Electric Quadrupole Moment                        (Debye-Å)\n";
    out << "XX= " << std::setw(20) << elecQuadrupole[0][0] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
      << " XY= " << std::setw(20) << elecQuadrupole[0][1] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
      << " XZ= " << std::setw(20) << elecQuadrupole[0][2] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << "\n";
    out << "YX= " << std::setw(20) << elecQuadrupole[1][0] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << " YY= " << std::setw(20) << elecQuadrupole[1][1] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << " YZ= " << std::setw(20) << elecQuadrupole[1][2] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << "\n";
    out << "ZX= " << std::setw(20) << elecQuadrupole[2][0] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << " ZY= " << std::setw(20) << elecQuadrupole[2][1] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << " ZZ= " << std::setw(20) << elecQuadrupole[2][2] * pow(BOHR_TO_ANGSTROM,2) * CHARGE_AU_TO_DEBYE
        << "\n\n";

    out << "\n\n";
    out << "Electric Octupole Moment                          (Debye-Å²)\n";
    for(int i=0; i<3; ++i) {
      for(int j=0; j<3; ++j) {
        for (int k=0; k<3; ++k) {
          out << static_cast<char>('X'+i) << static_cast<char>('X'+j) << static_cast<char>('X'+k) << "= "
              << std::setw(20) << elecOctupole[i][j][k] * pow(BOHR_TO_ANGSTROM,3) * CHARGE_AU_TO_DEBYE;
        }
        out << "\n";
      }
    }

    out << bannerEnd << std::endl;

  }

  




}
