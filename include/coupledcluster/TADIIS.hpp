/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2021 Li Research Group (University of Washington)
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
#pragma once

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/cqlinalg_config.hpp>

#include <util/files.hpp>
#include <memmanager.hpp>
#include <coupledcluster/EOMCCSDVector.hpp>

namespace ChronusQ {

    /**
     *   \brief A general DIIS class for extrapolations using Tiled Array
     *    based on a series of error metrics stored in memory.
     *
     */

template <typename T>
class DIISTA {

public:

    /**
       *  TADIIS Constructor. Constructs a TADIIS object
       *
       *  \param [in]  maxdiis     Maximum number of error vectors
       */
    DIISTA(int maxdiis);

    ~DIISTA();

    /**
       *  WriteVector. Set current solution
       *
       *  \param [in]  vector      array of tiled array object representing current solution
       */
    void WriteVector(const EOMCCSDVector<T> &vector);

    /**
       *  WriteErrorVector. Set current error vector
       *
       *  \param [in]  vector      array of tiled array object representing current error vector
       */
    void WriteErrorVector(const EOMCCSDVector<T> &vector);

    /**
       *  extrapolate. solve DiskDIIS equations and extrapolation solution
       *
       *  \param [in]  vector      Vector of tiled array solution for setting extrapolation
       *  \param [out] vector      Vector of extrapolated tiled array amplitudes
       */
    void Extrapolate(EOMCCSDVector<T> &vector);

    /// restart diis
    void restart();

protected:

    /// stores history of amplitude matrix
    std::vector<EOMCCSDVector<T>> hist;

    /// stores history of error vector
    std::vector<EOMCCSDVector<T>> histerr;

    /// determine diis expansion coefficients
    void DIISCoefficients(int nvec);

    /// Compute dot product of error vectors i and j
    T diis_dot(int i, int j);

    /// Compute norm of error vector i
    T diis_norm(int i);

    /// maximum number of diis vectors
    int maxdiis_;          

    /// dimension of each diis vector
    T * diisvec_;

    /// The error matrix
    T * errmtx;

    /// current number of diis vectors
    int diis_iter_;         

    /// diis vector to be replaced
    int replace_diis_iter_;
};

    template <typename T>
    DIISTA<T>::DIISTA(int maxdiis) {

        maxdiis_           = maxdiis;
        diis_iter_         = 0;
        replace_diis_iter_ = 1;
        diisvec_           = CQMemManager::get().malloc<T>(maxdiis_+1);
        errmtx             = CQMemManager::get().malloc<T>(maxdiis_*maxdiis_);
        memset((void*)errmtx,'\0',maxdiis_*maxdiis_*sizeof(T));

        hist.reserve(maxdiis_ + 1);
        histerr.reserve(maxdiis_ + 1);
    }

    template <typename T>
    DIISTA<T>::~DIISTA() {
        CQMemManager::get().free(diisvec_);
        CQMemManager::get().free(errmtx);
    }

// reset diis solver (without freeing memory).
    template <typename T>
    void DIISTA<T>::restart() {
        diis_iter_         = 0;
        replace_diis_iter_ = 1;
    }

// Store current solution vector to DIIS history
    template <typename T>
    void DIISTA<T>::WriteVector(const EOMCCSDVector<T> &vector){

        if ( diis_iter_ <= maxdiis_ ){
            hist.push_back(vector);
        }
        else{
            hist[replace_diis_iter_] = vector;
        }
    }

// Store current error vector to DIIS history
    template <typename T>
    void DIISTA<T>::WriteErrorVector(const EOMCCSDVector<T> &vector){

        if ( diis_iter_ <= maxdiis_ ){
            histerr.push_back(vector);
        }
        else{
            histerr[replace_diis_iter_] = vector;
        }
    }

// Perform DIIS extrapolation.
    template <typename T>
    void DIISTA<T>::Extrapolate(EOMCCSDVector<T> &vector){

        if ( diis_iter_ > 1 ) {

            // Compute coefficients for the extrapolation
            DIISCoefficients( std::min(diis_iter_, maxdiis_) );

            int max = diis_iter_;
            if (max > maxdiis_) max = maxdiis_;

            //zero out amplitudes in vector to write in extrapolated amplitudes
            vector.scale(0.0);

            for (int j = 1; j <= max; j++){

                // Accumulate extrapolated vector.
                vector.axpy(diisvec_[j-1], hist[j]);
            }
        }

        if (diis_iter_ <= maxdiis_){
            diis_iter_++;
        }
        /*else {
            // If we already have maxdiis_ vectors, choose the one with
            // the largest error as the one to replace.
            int jmax   = 0;
            T max = -1.0e99;
            for (int j = 0; j < maxdiis_; j++){
                T nrm = diis_norm(j+1);
                if ( nrm > max ) {
                    max  = nrm;
                    jmax = j+1;
                }
            }
            replace_diis_iter_ = jmax;
        }*/
        else if (replace_diis_iter_ < maxdiis_) {
            replace_diis_iter_++;
        }else {
            replace_diis_iter_ = 1;
        }
    }

// Evaluate extrapolation coefficients for DIIS.
    template <typename T>
    void DIISTA<T>::DIISCoefficients(int nvec){

        // Allocate memory for small matrices/vectors.
        int64_t * ipiv    = CQMemManager::get().malloc<int64_t>(nvec+1);
        T * A    = CQMemManager::get().malloc<T>((nvec+1)*(nvec+1));
        T * B    = CQMemManager::get().malloc<T>(nvec+1);
        memset((void*)A,'\0',(nvec+1)*(nvec+1)*sizeof(T));
        memset((void*)B,'\0',(nvec+1)*sizeof(T));
        B[nvec] = -1.0;

        // Reshape the error matrix, in case its dimension is less than maxdiis_.
        for (int i = 0; i < nvec; i++){
            for (int j = 0; j < nvec; j++){
                A[i*(nvec+1)+j] = errmtx[i*maxdiis_+j];
            }
        }

        if (nvec <= 3) {
            // At early iterations, just build the whole matrix.
            for (int i = 0; i < nvec; i++) {
                for (int j = i+1; j < nvec; j++){
                    T sum  = diis_dot(i+1, j+1);
                    A[i*(nvec+1)+j] = sum;
                    A[j*(nvec+1)+i] = sum;
                }
                T sum  = diis_dot(i+1, i+1);
                A[i*(nvec+1)+i] = sum;
            }
        }else {
            // At later iterations, don't build the whole matrix.
            // Just replace one row/column.

            // Which row/column will be replaced?
            int i = nvec < maxdiis_ ? nvec - 1 : replace_diis_iter_-1;

            for (int j = 0; j < nvec; j++){
                T sum  = diis_dot(i+1, j+1);
                A[i*(nvec+1)+j] = sum;
                A[j*(nvec+1)+i] = sum;
            }
        }

        int j = nvec;
        for (int i = 0; i < (nvec+1); i++){
            A[j*(nvec+1)+i] = -1.0;
            A[i*(nvec+1)+j] = -1.0;
        }
        A[(nvec+1)*(nvec+1)-1] = 0.0;

        // Write error matrix for next iteration.
        for (int i = 0; i < nvec; i++){
            for (int j = 0; j < nvec; j++){
                errmtx[i*maxdiis_+j] = A[i*(nvec+1)+j];
            }
        }

        // Solve the set of linear equations for the extrapolation coefficients
        int nrhs = 1;
        int lda  = nvec+1;
        int ldb  = nvec+1;
        lapack::gesv(nvec+1,nrhs,A,lda,ipiv,B,ldb);
        std::copy_n(B,nvec,diisvec_);

        CQMemManager::get().free(A);
        CQMemManager::get().free(B);
        CQMemManager::get().free(ipiv);
    }

    template <typename T>
    T DIISTA<T>::diis_dot(int i, int j){
        T res = histerr[i].dot(histerr[j]);
	TA::get_default_world().gop.fence();
	return res;
    }

    template <typename T>
    T DIISTA<T>::diis_norm(int i){
        return histerr[i].norm();
    }

} // end of namespace
