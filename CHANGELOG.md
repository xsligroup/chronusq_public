                     CHRONUSQ RELEASE HISTORY


  This file simply provides a brief overview of the release history of the
  ChronusQ program as released by the Li Research Group at the University
  of Washington. For the most up-to-date record of functionality, etc,
  please refer to the [ChronusQ Wiki](https://urania.chem.washington.edu/chronusq/chronusq_public/wikis/home).


  FORMAT: YYYY-MM-DD

  - 2022-03-07 0.6.0 (BETA)
    - Added Functionality:
      - Add many varieties of configuration interaction
        - Complete active space (CASCI and CASSCF) for all Hamiltonians
        - Restricted active space (RASCI) for 2 component Hamiltonians
      - Add nuclear electronic orbital Hartree-Fock and DFT (SCF & Real Time)
      - Add quasi-Newton SCF optimization methods (BFGS/SR1)
      - Add efficient Cholesky decomposition based RI approximation algorithms
      - Allow classical nuclei to carry variable and fractional charges
      - Add option to swap molecular orbitals in initial guess
      - Add overall program timer and some summary sections
      - Improve molecular orbital printing
      - Allow guesses from references of different types
      - Add time dependent population analysis in real time methods
    - Internal Refactoring:
      - Minimum C++ standard now supported is C++ 17
      - Support Libcint as a possible integral module
      - Move linear algebra interface to BLAS++/LAPACK++
    - Bug Fixes:
      - Remove QZ convergence failure
      - Corrected spelling of "Schwarz"
      - Add missing template instantiation in GCC 10 compilation
      - Fix segfault when reading basis in input file
  <br>

  - 2020-12-08 0.5.0 (BETA)
    - Added Functionality:
      - Add Relativistic (X2C) Coupled Cluster Singles and Doubles (MPI only)
      - Add Linear Response for two component DFT (X2C/GKS)
    - Internal Refactoring:
      - Change CI pipelines to build and deploy Docker images to Docker Hub
  <br>

  - 2020-10-20 0.4.0 (BETA)
    - Added Functionality:
      - Add Resolution of identity (RI) approximation for all methods
      - Add Restricted Open shell Hartree Fock (ROHF) wavefunction optimization
      - Accept Gaussian formatted checkpoint files as an initial guess
      - Add SKIP option to SCF.ALG for post-SCF methods
      - Add Peterson correlation consistent relativistic bases
      - Add examples of input and output to documentation
    - Internal Refactoring:
      - Change AOIntegrals class to Integrals to support RI and AO/MO bases
      - Rework documentation through wiki
      - Add transpose to GEMM and MatAdd
      - Improve CMake discovery of Libint
      - Update README and author list
    - Bugfixes:
      - Fail on bad input to RESPONSE.TYPE input section
      - Error on empty MOLECULE.GEOM input section
      - Enforce multiplicity consistent with number of electrons
      - Ensure input file is not the same as output file
      - Error on unspecified BASIS input section
  <br>

  - 2020-07-17 0.3.3 (BETA)
    - Support multiplicities up to element Es for SAD guess
    - Fix incorrect multiplicities for elements S, P, and Mn in SAD guess
    - Restructure SingleSlater class to have CoreHBuilder and FockBuilder objects
    - Fix incorrect reference data for X2C RT unit tests
    - Fix uncontracted basis by removing duplicate primitives
    - Allow arbitrary indentation now allowed in input file
    - Fix mixed GEMM call between double and dcomplex types
    - Fix input basis case sensitivity when parsing element labels
    - Fix simultaneous file read failure from GPLHR_MPI tests
    - Fix include guards
    - Bump OpenBLAS -> v0.3.9
  <br>

  - 2020-04-25 0.3.2 (BETA)
    - Add ability to read user specified basis from input file or path to file
    - Add B3PW91 functional
    - Add CMake flag to enable linking to an external OpenMP
    - Make compatible with compilation on Mac with the default compilers
    - Fix sign error in DFT gradient evaluation
    - Add templates for issue and merge requests
  <br>

  - 2020-03-12 0.3.1 (BETA)
    - Add ability to restart interrupted RT jobs
    - Add ability to turn the SCF field off in RT jobs
    - Fix sign bug in X2C calculations
    - Fix CXXBLACS linking error for MPI builds
    - Add MPI build into continuous integration
    - Fix bug in reading basis sets using Fortran float notation
  <br>

  - 2019-05-06 0.3.0 (BETA)
    - Implementation of Explicit Magnus 2nd order step in RT module
    - Default parameter changes in GPLHR and SCF
    - CI through GitLab
    - Clang 9+ compatible
  <br>

  - 2018-11-28 0.2.1 (BETA)
    - Removed Boost depedency
      - Switched Boost::Test -> GTest for UT system
      - Created in house segregated storage engine to remove boost::simple_segregated_storage
      - Removed Boost::Filesystem for file search / copy
    - Fixed HDF5 link in the presense of static libraries
      - This allows CQ to be compiled statically
    - Fixed dependency tree to allow for parallel make in openblas builds    
    - Misc GCC 6 comptibility fixes
    - Parallel (SMP + MPI) GIAO Fock builds
    - Direct GIAO Fock builds in RT module
    - Bump Libint -> 2.5.0-beta
    - Bump CMake  -> 3.11
  <br>

  - 2018-07-13 0.2.0 (BETA)
    - Full integration of GIAO basis set into SCF and RT modules
    - Added RESPONSE module
      - Supports the PolarizationPropagator (TDDFT/TDHF) and ParticleParticlePropagator (pp-RPA/pp-TDA/hh-TDA)
      - RESIDUE -> eigen decomposition
      - (D)FDR -> (damped) frequency depenent response
      - GPLHR for partial diagonalization
        - Supports arbitrary energy domain for diagonalization
    - Full integration of MPI functionality throughout (CQ_ENABLE_MPI)
      - Using MXX for C++11 MPI bindings
      - Using CXXBLACS for C++ interface to BLACS type functionality
      - Integration of ScaLAPACK into RESPONSE module
    - Bump Libint -> 2.4.2
    - Bump Libxc  -> 4.0.4
    - Added support for coverage checks (CQ_ENABLE_COVERAGE)
    - Various logic checks / bug fixes
  <br>

  - 2017-09-01: 0.1.0 (BETA)
    - Complete overhaul of ChronusQ development stream (new repo)
    - Currently tested functionality:
      - Full support for Hartree-Fock and Kohn-Sham references
      - SCF general to Restricted (R), Unrestricted (U), and Generalized (G) references
      - Real-time propagation of R/U/G references
      - X2C Relativistic references (both HF and KS)
      - Full integration of OpenMP in performance critical code
      - Support for both INCORE (full ERI) and DIRECT integral contraction schemes
