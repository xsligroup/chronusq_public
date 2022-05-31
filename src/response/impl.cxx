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

#include <response/impl.hpp>

namespace ChronusQ {

  template class ResponseTBase<double>;
  template class ResponseTBase<dcomplex>;

  template class PolarizationPropagator<HartreeFock<double,double>>;
  template class PolarizationPropagator<HartreeFock<dcomplex,double>>;
  template class PolarizationPropagator<HartreeFock<dcomplex,dcomplex>>;
  template class PolarizationPropagator<KohnSham<double,double>>;
  template class PolarizationPropagator<KohnSham<dcomplex,double>>;
  template class PolarizationPropagator<KohnSham<dcomplex,dcomplex>>;
  template class PolarizationPropagator<NEOSS<double,double>>;
  template class PolarizationPropagator<NEOSS<dcomplex,double>>;
  template class PolarizationPropagator<NEOSS<dcomplex,dcomplex>>;

  template class ParticleParticlePropagator<HartreeFock<double,double>>;
  template class ParticleParticlePropagator<HartreeFock<dcomplex,double>>;
  template class ParticleParticlePropagator<HartreeFock<dcomplex,dcomplex>>;
  template class ParticleParticlePropagator<KohnSham<double,double>>;
  template class ParticleParticlePropagator<KohnSham<dcomplex,double>>;
  template class ParticleParticlePropagator<KohnSham<dcomplex,dcomplex>>;

  template void ResponseTBase<double>::
    runFullFDR(FDResponseResults<double,double> &);
  template void ResponseTBase<double>::
    runFullFDR(FDResponseResults<double,dcomplex> &);
  template void ResponseTBase<dcomplex>::
    runFullFDR(FDResponseResults<dcomplex,dcomplex> &);

  template void ResponseTBase<double>::
    runIterFDR(FDResponseResults<double,double> &,
               std::function< void(size_t,double,SolverVectors<double>&,SolverVectors<double>&) > & );
  template void ResponseTBase<double>::
    runIterFDR(FDResponseResults<double,dcomplex> &,
               std::function< void(size_t,dcomplex,SolverVectors<dcomplex>&,SolverVectors<dcomplex>&) > & );
  template void ResponseTBase<dcomplex>::
    runIterFDR(FDResponseResults<dcomplex,dcomplex> &,
               std::function< void(size_t,dcomplex,SolverVectors<dcomplex>&,SolverVectors<dcomplex>&) > & );

  template void ResponseTBase<double>::
    fdrObservables(FDResponseResults<double,double> &);
  template void ResponseTBase<double>::
    fdrObservables(FDResponseResults<double,dcomplex> &);
  template void ResponseTBase<dcomplex>::
    fdrObservables(FDResponseResults<dcomplex,dcomplex> &);



  template void ResponseTBase<double>::
    printRF(FDResponseResults<double,double>&,std::ostream&);
  template void ResponseTBase<double>::
    printRF(FDResponseResults<double,dcomplex>&,std::ostream&);
  template void ResponseTBase<dcomplex>::
    printRF(FDResponseResults<dcomplex,dcomplex>&,std::ostream&);




  template void PolarizationPropagator<SingleSlater<double,double>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<double>>,
    SINGLESLATER_POLAR_COPT
  );
  template void PolarizationPropagator<SingleSlater<double,double>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT
  );
  template void PolarizationPropagator<SingleSlater<dcomplex,double>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT
  );
  template void PolarizationPropagator<SingleSlater<dcomplex,dcomplex>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT
  );

  template void PolarizationPropagator<HartreeFock<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<double>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<HartreeFock<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<HartreeFock<dcomplex,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<HartreeFock<dcomplex,dcomplex>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<NEOSS<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<double>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<NEOSS<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );

  template void PolarizationPropagator<NEOSS<dcomplex,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<NEOSS<dcomplex,dcomplex>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );

  template void ParticleParticlePropagator<SingleSlater<double,double>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<double>>
  );
  template void ParticleParticlePropagator<SingleSlater<double,double>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );
  template void ParticleParticlePropagator<SingleSlater<dcomplex,double>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );
  template void ParticleParticlePropagator<SingleSlater<dcomplex,dcomplex>>::formLinearTrans_incore_impl(
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );


  template void ParticleParticlePropagator<HartreeFock<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<double>>
  );
  template void ParticleParticlePropagator<HartreeFock<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );
  template void ParticleParticlePropagator<HartreeFock<dcomplex,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );
  template void ParticleParticlePropagator<HartreeFock<dcomplex,dcomplex>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );



  template void PolarizationPropagator<KohnSham<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<double>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<KohnSham<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<KohnSham<dcomplex,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );
  template void PolarizationPropagator<KohnSham<dcomplex,dcomplex>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>,
    SINGLESLATER_POLAR_COPT,
    bool
  );

  template void ParticleParticlePropagator<KohnSham<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<double>>
  );
  template void ParticleParticlePropagator<KohnSham<double,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );
  template void ParticleParticlePropagator<KohnSham<dcomplex,double>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );
  template void ParticleParticlePropagator<KohnSham<dcomplex,dcomplex>>::formLinearTrans_direct_impl(
    MPI_Comm,
    std::vector<RESPONSE_CONTRACTION<dcomplex>>
  );


  template void PolarizationPropagator<SingleSlater<double,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,double*,double*);
  template void PolarizationPropagator<SingleSlater<double,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,dcomplex*,dcomplex*);
  template void PolarizationPropagator<SingleSlater<dcomplex,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,dcomplex*,dcomplex*);
  template void PolarizationPropagator<SingleSlater<dcomplex,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,double*,double*);
  template void PolarizationPropagator<SingleSlater<dcomplex,dcomplex>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,dcomplex*,dcomplex*);
  template void PolarizationPropagator<SingleSlater<dcomplex,dcomplex>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,double*,double*);

  template void ParticleParticlePropagator<SingleSlater<double,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,double*,double*);
  template void ParticleParticlePropagator<SingleSlater<double,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,dcomplex*,dcomplex*);
  template void ParticleParticlePropagator<SingleSlater<dcomplex,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,dcomplex*,dcomplex*);
  template void ParticleParticlePropagator<SingleSlater<dcomplex,double>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,double*,double*);
  template void ParticleParticlePropagator<SingleSlater<dcomplex,dcomplex>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,dcomplex*,dcomplex*);
  template void ParticleParticlePropagator<SingleSlater<dcomplex,dcomplex>>::printResMO_impl(
      std::ostream&,size_t,double*,std::vector<std::pair<std::string,double*>>,double*,double*);

  template void PolarizationPropagator<SingleSlater<double,double>>::preConditioner(size_t nVec,
      double shift, SolverVectors<double> &V, SolverVectors<double> &AV);
  template void PolarizationPropagator<SingleSlater<double,double>>::preConditioner(size_t nVec,
      dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV);
  template void PolarizationPropagator<SingleSlater<dcomplex,double>>::preConditioner(size_t nVec,
      dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV);
  template void PolarizationPropagator<SingleSlater<dcomplex,dcomplex>>::preConditioner(size_t nVec,
      dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV);

  template void ParticleParticlePropagator<SingleSlater<double,double>>::preConditioner(size_t nVec,
      double shift, SolverVectors<double> &V, SolverVectors<double> &AV);
  template void ParticleParticlePropagator<SingleSlater<double,double>>::preConditioner(size_t nVec,
      dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV);
  template void ParticleParticlePropagator<SingleSlater<dcomplex,double>>::preConditioner(size_t nVec,
      dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV);
  template void ParticleParticlePropagator<SingleSlater<dcomplex,dcomplex>>::preConditioner(size_t nVec,
      dcomplex shift, SolverVectors<dcomplex> &V, SolverVectors<dcomplex> &AV);
          
}; // namespace ChronusQ

