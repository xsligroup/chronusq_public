 /*
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
#pragma once

#include <mcscf.hpp>
#include <cibuilder/neo.hpp>
#include <cibuilder/neo/impl.hpp>
#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>
#include <mointstransformer/moranges.hpp>

namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void NEOMCWaveFunction<MatsT,IntsT>::computeOneRDM(size_t i)
    {
        this->ciBuilder->computeOneRDM(*this,this->CIVecs[i],this->ewfn_->oneRDM[i]);
        NEOCIBuilder->computePOneRDM(*this,this->CIVecs[i],this->pwfn_->oneRDM[i]);
    }

    template <typename MatsT, typename IntsT>
    void NEOMCWaveFunction<MatsT,IntsT>::computeOneRDM()
    {
        for(size_t i = 0; i < this->NStates; i++)
        {
            computeOneRDM(i);
        }
    }

    template <typename MatsT, typename IntsT>
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> NEOMCWaveFunction<MatsT,IntsT>::getOnePDM()
    {
        // Get the orbital offsets for the state
        size_t nInact = this->ewfn_->MOPartition.nInact;
        size_t nCorrO = this->ewfn_->MOPartition.nCorrO;
        size_t nAO = this->ewfn_->ref_->mo[0].nRows();
        MatsT* MO = this->ewfn_->ref_->mo[0].pointer() + nAO * nInact;

        std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> PDMs;
        PDMs.reserve(this->NStates);

        cqmatrix::Matrix<MatsT> SCR(nAO);
        for(size_t i = 0; i < this->NStates; i++)
        {
            cqmatrix::Matrix<MatsT> PDM(nAO);
            MatsT * rdm = this->ewfn_->oneRDM[i].pointer();
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nAO,nCorrO,nCorrO,MatsT(1.0),MO,nAO,rdm,nCorrO,0.0,SCR.pointer(),nAO);
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,nAO,nAO,nCorrO,MatsT(1.0),SCR.pointer(),nAO,MO,nAO,0.0,PDM.pointer(),nAO);
            PDMs.emplace_back(std::make_shared<cqmatrix::Matrix<MatsT>>(PDM));
        }
        return PDMs;
    }

    template <typename MatsT, typename IntsT>
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> NEOMCWaveFunction<MatsT,IntsT>::getPOnePDM()
    {
        // Get the orbital offsets for the state
        size_t nInact = this->pwfn_->MOPartition.nInact;
        size_t nCorrO = this->pwfn_->MOPartition.nCorrO;
        size_t nAO = this->pwfn_->ref_->mo[0].nRows();
        MatsT * MO = this->pwfn_->ref_->mo[0].pointer() + nAO * nInact;

        std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> PDMs;
        PDMs.reserve(this->NStates);

        cqmatrix::Matrix<MatsT> SCR(nAO);
        for(size_t i = 0; i < this->NStates; i++)
        {
            cqmatrix::Matrix<MatsT> PDM(nAO);
            MatsT* rdm = this->pwfn_->oneRDM[i].pointer();
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nAO,nCorrO,nCorrO,MatsT(1.0),MO,nAO,rdm,nCorrO,MatsT(0.0),SCR.pointer(),nAO);
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,nAO,nAO,nCorrO,MatsT(1.0),SCR.pointer(),nAO,MO,nAO,MatsT(0.0),PDM.pointer(),nAO);
            PDMs.emplace_back(std::make_shared<cqmatrix::Matrix<MatsT>>(PDM));
        }
        return PDMs;
    }

}; // namespace ChronusQ
