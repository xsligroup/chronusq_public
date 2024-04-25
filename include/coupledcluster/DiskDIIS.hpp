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
#pragma once

#include <cqlinalg/blas1.hpp>
#include <cqlinalg/solve.hpp>

#include <util/files.hpp>
#include <memmanager.hpp>

namespace ChronusQ {

  /**
   *   \brief A general DIIS class for extrapolations
   *    based on a series of error metrics stored on disk.
   *
   */

  template <typename T>
  class DiskDIIS{
  
    public:
 
      /**
       *  DiskDIIS Constructor. Constructs a DiskDIIS object
       *
       *  \param [in]  n           Dimension of extrapolation space. Two temporary arrays of length n will be created
       *  \param [in]  maxdiis     Maximum number of error vectors
       *  \param [in]  savFile     File manager
       */ 
      DiskDIIS(int n, int maxdiis, SafeFile savFile);
  
      ~DiskDIIS();
  
      /**
       *  getVector. Get current solution
       *
       *  \param [out] vector      Pointer to target array
       *  \param [in]  dim         Number of elements to be copied to vector
       *  \param [in]  off         Offset in full array (if multiple quantities extrapolated simultaneously)
       */ 
      void getVector(T * vector,int dim, int off);
  
      /**
       *  setVector. Set current solution
       *
       *  \param [in]  vector      Pointer to current solution 
       *  \param [in]  dim         Number of elements to be copied to DiskDIIS
       *  \param [in]  off         Offset in full array (if multiple quantities extrapolated simultaneously)
       */ 
      void setVector(T * vector,int dim, int off);
  
      /**
       *  setErrorVector. Set current error vector
       *
       *  \param [in]  vector      Pointer to current error vector
       *  \param [in]  dim         Number of elements to be copied to DiskDIIS
       *  \param [in]  off         Offset in full array (if multiple quantities extrapolated simultaneously)
       */ 
      void setErrorVector(T * vector,int dim, int off);
  
      /**
       *  extrapolate. solve DiskDIIS equations and extrapolation solution
       *
       *  \param [in]  vector      Pointer to current error vector
       *  \param [in]  dim         Number of elements to be copied to DiskDIIS
       *  \param [in]  off         Offset in full array (if multiple quantities extrapolated simultaneously)
       */ 
      void extrapolate();
  
    protected:
  
      /// file manager
      SafeFile savFile_;
  
      /// write current error vector to disk
      void write_error_vector(T * vector);
  
      /// write current solution to disk
      void write_vector(T * vector);
  
  
      /// DiskDIIS temporary storage 1
      T * tmp1_;
  
      /// DiskDIIS temporary storage 2
      T * tmp2_;
  
      /// determine DiskDIIS expansion coefficients
      void DIISCoefficients(int nvec);
  
      /// maximum number of vectors
      unsigned maxdiis_;
  
      /// array containing DiskDIIS expansion coefficients
      T * diisvec_;
  
      /// dimension of the solution/error vector
      unsigned dimdiis_;
  
      /// current number of DiskDIIS vectors
      unsigned diis_iter_;
  
      /// DiskDIIS vector to be replaced
      unsigned replace_diis_iter_;
  };
  
  template <typename T>
  DiskDIIS<T>::DiskDIIS(int n, int maxdiis, SafeFile savFile):
    savFile_(savFile)
  {
    dimdiis_           = n;
    maxdiis_           = maxdiis;
    diis_iter_         = 0;
    replace_diis_iter_ = 1;
    diisvec_           = CQMemManager::get().malloc<T>(maxdiis_+1);
    tmp1_              = CQMemManager::get().malloc<T>(dimdiis_);
    tmp2_              = CQMemManager::get().malloc<T>(dimdiis_);
  }
  
  template <typename T>
  DiskDIIS<T>::~DiskDIIS(){
    CQMemManager::get().free(diisvec_);
    CQMemManager::get().free(tmp1_);
    CQMemManager::get().free(tmp2_);
  }
  
  // Write current solution vector to disk
  template <typename T>
  void DiskDIIS<T>::write_vector(T * vector){
  
    // Name the entry in according to the current DiskDIIS iteration.
    // If we already have maxdiis_ vectors, then replace one.
  
    char * oldvector = CQMemManager::get().malloc<char>(1000);
  
    if ( diis_iter_ <= maxdiis_ ){
      sprintf(oldvector,"/DIIS/VECTOR%i",diis_iter_);
    }
    else{
      sprintf(oldvector,"/DIIS/VECTOR%i",replace_diis_iter_);
    }
  
    // Write the current solution vector.
    savFile_.safeWriteData(oldvector, vector, {dimdiis_});
  
    CQMemManager::get().free(oldvector);
  }
  
  template <typename T>
  void DiskDIIS<T>::getVector(T * vector,int dim, int off) {
    std::copy_n(tmp1_+off,dim,vector);
  }
  
  template <typename T>
  void DiskDIIS<T>::setVector(T * vector,int dim, int off) {
    std::copy_n(vector,dim,tmp1_+off);
  }
  
  template <typename T>
  void DiskDIIS<T>::setErrorVector(T * vector,int dim, int off) {
    std::copy_n(vector,dim,tmp2_+off);
  }
  
  // Write current error vector to disk 
  template <typename T>
  void DiskDIIS<T>::write_error_vector(T * vector){
  
    // Name the entry in according to the current DiskDIIS iteration.
    // If we already have maxdiis_ vectors, then replace one.
  
    char * evector = CQMemManager::get().malloc<char>(1000);
    if ( diis_iter_ <= maxdiis_ ){
      sprintf(evector,"/DIIS/ERROR%i",diis_iter_);
    }
    else{
      sprintf(evector,"/DIIS/ERROR%i",replace_diis_iter_);
    }
  
    if ( diis_iter_ == 0 ) {
      // On the first iteration, write an entry that will hold the error matrix
      T * temp = (T*)malloc(maxdiis_*maxdiis_*sizeof(T));
      memset((void*)temp,'\0',maxdiis_*maxdiis_*sizeof(T));
      savFile_.safeWriteData("/DIIS/ERROR_MATRIX", temp, {maxdiis_*maxdiis_});
      free(temp);
    }
  
    // Write current error vector.
    savFile_.safeWriteData(evector, vector, {dimdiis_});
  
    CQMemManager::get().free(evector);
  }
  
  // Perform DIIS extrapolation. solution in tmp1_
  template <typename T>
  void DiskDIIS<T>::extrapolate() {
    
    write_vector(tmp1_);
    write_error_vector(tmp2_);
   
    if ( diis_iter_ > 1 ) {
    
      // Compute coefficients for the extrapolation
      DIISCoefficients( diis_iter_ < maxdiis_ ? diis_iter_ : maxdiis_ );
   
      memset((void*)tmp1_,'\0',dimdiis_*sizeof(T));
    
      char * oldvector = CQMemManager::get().malloc<char>(1000);
    
      int max = diis_iter_;
      if (max > maxdiis_) max = maxdiis_;
    
      // Read each of the old vectors from disk.
      for (int j = 1; j <= max; j++){
    
        sprintf(oldvector,"/DIIS/VECTOR%i",j);
        savFile_.readData(oldvector,tmp2_);
   
        // Accumulate extrapolated vector.
        blas::axpy(dimdiis_,diisvec_[j-1],tmp2_,1,tmp1_,1);
      }
      CQMemManager::get().free(oldvector);
  
    }
    
    if (diis_iter_ <= maxdiis_){
      diis_iter_++;
    }
    //else {
    //    // If we already have maxdiis_ vectors, choose the one with
    //    // the largest error as the one to replace.
    //    int jmax   = 1;
    //    double max = -1.0e99;
    //    char * evector = CQMemManager::get().malloc<char>(1000);
    //    for (int j = 1; j <= maxdiis_; j++){
    //        sprintf(evector,"/DIIS/ERROR%i",j);
    //        savFile_.readData(evector,tmp2_);
    //        double nrm = sqrt(std::real(blas::dot(dimdiis_,tmp2_,1,tmp2_,1)));  // couldn't get TwoNorm template to work...
    //        if ( nrm > max ) {
    //            max  = nrm;
    //            jmax = j;
    //        }
    //    }
    //    replace_diis_iter_ = jmax;
    //    CQMemManager::get().free(evector);
    //}
    else if (replace_diis_iter_ < maxdiis_) {
      replace_diis_iter_++;
    }else {
      replace_diis_iter_ = 1;
    }
  
  }
  
  // Evaluate extrapolation coefficients for DiskDIIS.
  template <typename T>
  void DiskDIIS<T>::DIISCoefficients(int nvec){
  
    // Allocate memory for small matrices/vectors.
    int64_t * ipiv    = CQMemManager::get().malloc<int64_t>(nvec+1);
    
    T * temp = CQMemManager::get().malloc<T>(maxdiis_*maxdiis_);
    T * A    = CQMemManager::get().malloc<T>((nvec+1)*(nvec+1));
    T * B    = CQMemManager::get().malloc<T>(nvec+1);
    
    memset((void*)A,'\0',(nvec+1)*(nvec+1)*sizeof(T));
    memset((void*)B,'\0',(nvec+1)*sizeof(T));
    
    B[nvec] = -1.0;
    
    char * evector = CQMemManager::get().malloc<char>(1000);
    
    // Read in the previous error matrix, so we don't have 
    // to build the entire thing each iteration.
    
    savFile_.readData("/DIIS/ERROR_MATRIX",temp);
    
    // Reshape the error matrix, in case its dimension is less than maxdiis_.
    for (int i = 0; i < nvec; i++){
      for (int j = 0; j < nvec; j++){
        A[i*(nvec+1)+j] = temp[i*maxdiis_+j];
      }
    }
    
    if (nvec <= 3) {
      // At early iterations, just build the whole matrix.
      for (int i = 0; i < nvec; i++) {
        sprintf(evector,"/DIIS/ERROR%i",i+1);
        savFile_.readData(evector,tmp1_);
        for (int j = i+1; j < nvec; j++){
          sprintf(evector,"/DIIS/ERROR%i",j+1);
          savFile_.readData(evector,tmp2_);
          T sum  = blas::dot(dimdiis_,tmp1_,1,tmp2_,1);
          A[i*(nvec+1)+j] = sum;
          A[j*(nvec+1)+i] = sum;
        }
        T sum  = blas::dot(dimdiis_,tmp1_,1,tmp1_,1);
        A[i*(nvec+1)+i] = sum;
      }
    }else {
      // At later iterations, don't build the whole matrix.
      // Just replace one row/column.
    
      // Which row/column will be replaced?
      int i = nvec < maxdiis_ ? nvec - 1 : replace_diis_iter_ - 1;
    
      sprintf(evector,"/DIIS/ERROR%i",i+1);
      savFile_.readData(evector,tmp1_);
      for (int j = 0; j < nvec; j++){
        sprintf(evector,"/DIIS/ERROR%i",j+1);
        savFile_.readData(evector,tmp2_);
        T sum  = blas::dot(dimdiis_,tmp1_,1,tmp2_,1);
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
        temp[i*maxdiis_+j] = A[i*(nvec+1)+j];
      }
    }
    savFile_.safeWriteData("/DIIS/ERROR_MATRIX", temp, {maxdiis_*maxdiis_});
    
    CQMemManager::get().free(temp);
    CQMemManager::get().free(evector);
    
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

} // end of namespace


