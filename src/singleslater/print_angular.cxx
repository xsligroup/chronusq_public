#include <singleslater/print.hpp>

namespace ChronusQ {

    template<typename MatsT, typename IntsT>
    void SingleSlater<MatsT, IntsT>::printAngularProperties(std::ostream& out, bool withBanner) {
        QuantumBase::printAngularProperties(out, withBanner);

        // If the reference itself is X2C, note that there is a picture change error.
        if (this->nC == 2 && this->aoints_ != nullptr && this->aoints_->options_.x2cType != X2C_TYPE::OFF) {
            out << "  NOTE: X2C picture change error" << std::endl;
        }

        if (this->nC == 4) {
            out << "  NOTE: Angular momentum expectation values are not computed for 4C." << std::endl;
            return;
        }

        if (this->aoints_ == nullptr || this->aoints_->angmom == nullptr) return;

        auto res = computeAngularOrbitalRows(true);

        // out << "Angular property timing (s):\n";
        // out << "  Operator prep:      " << res.operatorBuildTime << '\n';
        // out << "  Expectation eval:   " << res.expectationTime << '\n';
        // out << "  Total:              " << res.totalTime << "\n\n";

        out << "Angular Momentum Expectation Values Per Orbital:\n" << bannerTop << "\n\n";
        out << std::fixed << std::setprecision(5);
        out << std::setw(8)  << std::left  << "MO"
            << std::setw(10) << std::right << "Occ"
            << std::setw(14) << std::right << "Energy"
            << std::setw(12) << std::right << "<s_z>"
            << std::setw(12) << std::right << "<l_z>"
            << std::setw(12) << std::right << "<j_z>"
            << std::setw(12) << std::right << "<s^2>"
            << std::setw(12) << std::right << "<l^2>"
            << std::setw(12) << std::right << "<j^2>" << std::endl;
        out << bannerMid << std::endl;

        const auto &rows = res.rows;
        for (size_t iRow = 0; iRow < rows.size(); ++iRow) {
            const auto &row = rows[iRow];
            out << std::setw(8) << std::left << std::to_string(iRow + 1);
            out << std::setw(10) << std::right << row.occupation;
            out << std::setw(14) << std::right << row.energy;
            out << std::setw(12) << std::right << row.sz;
            out << std::setw(12) << std::right << row.lz;
            out << std::setw(12) << std::right << row.jz;
            out << std::setw(12) << std::right << row.s2;
            out << std::setw(12) << std::right << row.l2;
            out << std::setw(12) << std::right << row.j2 << std::endl;
        }

        out << std::endl << bannerEnd << std::endl << std::endl;
    }

    template void SingleSlater<double, double>::printAngularProperties(std::ostream& out, bool withBanner);
    template void SingleSlater<dcomplex, double>::printAngularProperties(std::ostream& out, bool withBanner);
    template void SingleSlater<dcomplex, dcomplex>::printAngularProperties(std::ostream& out, bool withBanner);

}
