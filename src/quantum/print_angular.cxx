#include <quantum/base.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  void QuantumBase::printAngularProperties(std::ostream &out, bool withBanner) {

    if (withBanner)
      out << "\nAngular Momentum Information:\n" << bannerTop << "\n\n";

    out << std::fixed << std::setprecision(12);

    if (isX2CReference())
      out << "  NOTE: X2C picture change error" << std::endl;

    if (this->nC == 4) {
        out << "  NOTE: Angular momentum expectation values are not computed for 4C." << std::endl;
      return;
    }

    // Orbital angular momentum expectation values
    out << "  <Lx> = " << std::setw(10) << LExpect[0] << std::endl;
    out << "  <Ly> = " << std::setw(10) << LExpect[1] << std::endl;
    out << "  <Lz> = " << std::setw(10) << LExpect[2] << std::endl;
    out << "  <L^2> = " << std::setw(10) << LSq << std::endl;
    out << "  Orbital quantum number = " << std::setw(10) << LQuantNum << std::endl;
    out << std::endl;

    out << "  <SL> = " << std::setw(9) << SL << std::endl;
    out << "  <LS> = " << std::setw(9) << LS << std::endl;
    out << std::endl;

    // Total angular momentum expectation values
    out << "  <Jx> = " << std::setw(10) << JExpect[0] << std::endl;
    out << "  <Jy> = " << std::setw(10) << JExpect[1] << std::endl;
    out << "  <Jz> = " << std::setw(10) << JExpect[2] << std::endl;
    out << "  <J^2> = " << std::setw(10) << JSq << std::endl;
    out << "  Total angular momentum quantum number = " << std::setw(10) << JQuantNum << std::endl;
    out << std::endl;

    if (withBanner)
      out << "\n" << bannerEnd << std::endl;
  }

}
