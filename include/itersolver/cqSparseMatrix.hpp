#pragma once
#ifdef CQ_ENABLE_SPARSE

#include <vector>
#include <algorithm>
#ifdef CQ_HAS_PARALLEL_STL
  #include <execution>
#endif
#include <util/math.hpp>
#include <memmanager.hpp>
#include <cerr.hpp>
#include <boost/unordered/concurrent_flat_map.hpp>
#include <boost/unordered/unordered_flat_map.hpp>
#include <boost/sort/sort.hpp>

namespace ChronusQ {

   template <class T>
   class MyAlloc {
     public:
       // type definitions
       typedef T        value_type;
       typedef T*       pointer;
       typedef const T* const_pointer;
       typedef T&       reference;
       typedef const T& const_reference;
       typedef std::size_t    size_type;
       typedef std::ptrdiff_t difference_type;

       // rebind allocator to type U
       template <class U>
       struct rebind {
         typedef MyAlloc<U> other;
       };

       // return address of values
       pointer address (reference value) const {
         return &value;
       }
       const_pointer address (const_reference value) const {
         return &value;
       }

       /* constructors and destructor
        * - nothing to do because the allocator has no state
        */
       MyAlloc() throw() {}
       MyAlloc(const MyAlloc&) throw() {}
       template <class U>
       MyAlloc (const MyAlloc<U>&) throw() {}
       ~MyAlloc() throw() {}

       // return maximum number of elements that can be allocated
       size_type max_size () const throw() {
         return std::numeric_limits<std::size_t>::max() / sizeof(T);
       }

       // allocate but don't initialize num elements of type T
       pointer allocate (size_type num, const void* = 0) {
	 return CQMemManager::get().threadSafeMalloc<value_type>(num);
       }

       // initialize elements of allocated storage p with value value
       void construct (pointer p, const T& value) {
         // initialize memory with placement new
         new((void*)p)T(value);
       }

       // destroy elements of initialized storage p
       void destroy (pointer p) {
         // destroy objects by calling their destructor
         p->~T();
       }

       // deallocate storage p of deleted elements
       void deallocate (pointer p, size_type num) {
	 CQMemManager::get().threadSafeFree(p, num);
       }
   };

   // return that all specializations of this allocator are interchangeable
   template <class T1, class T2>
   bool operator== (const MyAlloc<T1>&, const MyAlloc<T2>&) throw() {
     return true;
   }
   template <class T1, class T2>
   bool operator!= (const MyAlloc<T1>&, const MyAlloc<T2>&) throw() {
     return false;
   }

//Forward declare LLSparseMatrix to use as friend class to HashSparseMatrix
template<typename MatsT>
class LLSparseMatrix;

/*
 * A class representing a sparse matrix using a thread-safe hash table. This allows for efficient random access
 */
template <typename MatsT>
class HashSparseMatrix {

  //private:
  public:
    std::vector<
      boost::unordered::concurrent_flat_map<size_t,
          MatsT,
          boost::hash<const size_t>,
          std::equal_to<const size_t>,
          MyAlloc<std::pair<const size_t, MatsT>>
      >
    > data;
    
  public:

    HashSparseMatrix(size_t nVec) {
      data.resize(nVec);
    }


    const MatsT get(size_t row, size_t col) const {
      MatsT ret = MatsT(0.);
      bool found = false;

      data[col].cvisit(row, [&](const auto& pair) {
        ret = pair.second;
        found = true;
    });

      return found? ret : MatsT(0.);
    }

    const auto& getCol(size_t col) {
      return data[col];
    }

    void update(size_t row, size_t col, MatsT value) {

      data[col].insert_or_visit(std::make_pair(row, value), [&value](auto& pair){pair.second += value;});
    }

    void divide(size_t row, size_t col, MatsT value) {

      data[col].insert_or_visit(std::make_pair(row, value), [&value](auto& pair){pair.second /= value;});
    }

    void add(HashSparseMatrix<MatsT>& other, double eps = 0) {
      for (auto j = 0ul; j < data.size(); ++j) {
        boost::unordered_flat_map serialOther = std::move(other.data[j]);
	for(auto it = serialOther.cbegin(); it != serialOther.cend(); it++) {
	  //if(std::abs(it->second) > eps)
            data[j].insert_or_visit(*it, [&it, &eps](auto& pair){pair.second += it->second; });
	    data[j].erase_if(it->first, [&eps](auto& pair){return std::abs(pair.second) < eps;});
	}
      }

    }


    //if(std::abs(pair.second) < eps){pair.second = MatsT(0.);}



    /*
    // davidson sparse preconditioner w/o Hamiltonian diagonal build
    void preCondition(HashSparseMatrix<MatsT>& other, 
    const LLSparseMatrix<MatsT>& myC,
    const CategoricalSpace& categoricalSpace, size_t iCat, size_t iCatOffset, const DeterminantFactory& detFactory, const IntegralsCollection& moints,
    const dcomplex * curEigenvalues,
    size_t nVec, double C_eps) {
   
	   const double small = 1e-12;
        const auto& hCore_tt = *(this->moints_.template getIntegral<DASOnePInts, MatsT>("hCore_tt"));
        const auto& antiSymmetricERI_ttuu = *(this->moints_.template getIntegral<OnePInts, MatsT>("antiSymmetricERI_ttuu"));


    // main loop
    const auto& nOrbs = this->detFactory_.nOrbitalsInEachSpace();
    const size_t nCorrE = this->detFactory_.nTotalCorrElectrons();

    //auto diagHLocalView = ketCategoricalSpace.createLocalCIVectorsView(diagH, 0ul, 1ul);

    const auto& cat = dynamic_cast<const FullDeterminantCategory&>(*ketCategoricalSpace.getCategory(iCat));
    //MatsT *catDiagH = diagHLocalView.getCategoryPointer(i);
    const auto nEs = braCategory.SpaceOccupations();
    const size_t nDets = braCategory.nDeterminants();
    const size_t nDetsPerThread = std::ceil(double(nDets) / GetNumThreads());

    //#pragma omp parallel default(shared)
    //{
      std::vector<size_t> occ(nCorrE);
      auto detCatGen = braCategory.template generator<uint64_t>();
      //size_t iBegin = nDetsPerThread * GetThreadID();
      //size_t iEnd   = std::min(nDets, iBegin + nDetsPerThread);
      //size_t iBegin = 0;
      //size_t iEnd   = nDets;

      for (auto j = 0ul; j < data.size(); ++j) {


      boost::unordered_flat_map serialOther = std::move(other.data[j]);
      for(auto it = serialOther.cbegin(); it != serialOther.cend(); it++) {
          //data[j].insert_or_visit(*it, [&it](auto& pair){pair.second += it->second;});

      detCatGen.visitDeterminants(it->first - iCatOffset, it->first - iCatOffset + 1,
          [&] (size_t addr, const auto& dets)  {
            determinantsToOccs(dets, nEs, nOrbs, occ);
            MatsT tmp = MatsT(0.);
            for(const auto& t : occ) {
              tmp += hCore_tt(t, 0);
              for(const auto& u : occ) {
                tmp += MatsT(0.5) * antiSymmetricERI_ttuu(t, u);
              }
            }

	    double P = std::real(tmp - std::real(curEigenvalues[j]);
	    // element needs to be dicarded; ASSUME HASH INTERMEDIATES ALREADY CONTAIN -S*EIG[j]
            if(std::abs(P) < small or std::abs( (it->second - std::real(curEigenvalues[j]) * myC.get(addr, j)) / P ) < eps) {
	      it->second = MatsT(0.)
	    }
	    else {
	      it->second =  (it->second - std::real(curEigenvalues[j]) * myC.get(addr, j)) / P;
	      other.data[j].insert_or_visit(*it, [&it](auto& pair){pair.second += it->second;});
	    }
          }
      );

      }  //it loop
      } // j loop
    //} // parallel region



  }

*/

    void reserve(size_t nnz, size_t j) {
        data[j].reserve(nnz);
    }

    void setMaxLoadFactor(double factor) {
      for(size_t i = 0; i < data.size(); i++)
        data[i].max_load_factor(factor);
    }
   
    void setZero() {
      #pragma omp parallel for
      for(size_t j = 0; j < data.size(); j++) {
        data[j].clear();
      }
    }

    friend class LLSparseMatrix<MatsT>;
};


template <typename MatsT>
class LLSparseMatrix {

  private:

    std::vector<std::vector<std::pair<size_t, MatsT>, MyAlloc<std::pair<size_t, MatsT>>  >> data;

  public:

    LLSparseMatrix(size_t nVec = 0) {
      data.resize(nVec);
    }
 
    void resize(size_t nVec) {
      data.resize(nVec);
    }

    size_t cols() const {
      return data.size();
    }

    auto& getCol(size_t col) {
      return data[col];
    }

    const auto& getCol(size_t col) const {
      return data[col];
    }

    size_t nonZeros() const {
      size_t nz = 0;
      #pragma omp parallel for reduction(+:nz)
      for(size_t j = 0; j < data.size(); j++) {
        nz += data[j].size();
      }
      return nz;
    }

    size_t nonZeros(size_t shift, size_t nVec = 1) const {
      size_t nz = 0;
      #pragma omp parallel for reduction(+:nz)
      for(size_t j = shift; j < shift + nVec; j++) {
        nz += data[j].size();
      }
      return nz;
    }

    void setZero() {
      #pragma omp parallel for
      for(size_t j = 0; j < data.size(); j++) {
	data[j].clear();
      }
    }

    void setZero(size_t shift, size_t nVec) {
      #pragma omp parallel for
      for(size_t j = shift; j < shift + nVec; j++) {
	data[j].clear();
      }
    }

    void prune(double eps) {
      #pragma omp parallel for
      for(auto& col : data) {
	std::erase_if(col, [&eps](const auto& element){return std::abs(element.second) < eps;});
	col.shrink_to_fit();
      }
    }
    void prune(double eps, size_t col) {
      std::erase_if(data[col], [&eps](const auto& element){return std::abs(element.second) < eps;});
      data[col].shrink_to_fit();
    }

    void copy(LLSparseMatrix& other, size_t shift, size_t nVec) const {
      #pragma omp parallel for
      for(size_t j = shift; j < shift + nVec; j++) {
        other.data[j] = data[j];
      }
    }

    void set_from_other(const LLSparseMatrix& other, size_t shiftA, size_t shiftB, size_t nVec) {
      #pragma omp parallel for
      for(size_t j = 0; j < nVec; j++) {
        data[j + shiftA] = other.data[j + shiftB];
      }
    }

    void swap(LLSparseMatrix& other, size_t shiftA, size_t shiftB, size_t nVec) {
      #pragma omp parallel for
      for(size_t j = 0; j < nVec; j++) {
        other.data[j + shiftA].swap(data[j + shiftB]);
      }
    }

    auto cbegin(size_t col) const {
      return data[col].cbegin();
    }
    auto cend(size_t col) const {
      return data[col].cend();
    }
    
    auto begin(size_t col) {
      return data[col].begin();
    }
    auto end(size_t col) {
      return data[col].end();
    }

    void sortedInsert(size_t row, size_t col, MatsT value) {
      if(data[col].size() > 0 and row < data[col].back().first) {
        CErr("unordered insert is not allowed");
      }
      if (value != MatsT(0.)) {
        data[col].push_back(std::make_pair(row, value));
      }
    }

    MatsT get(size_t row, size_t col) const {
      auto it = std::lower_bound(data[col].cbegin(), data[col].cend(), std::make_pair(row, MatsT(0.)), 
	  [](const auto& a, const auto& b){return a.first < b.first;}
      );

      if(it != data[col].cend() and (row == it->first))
	return it->second;

      return MatsT(0.); // Default value for non-existent elements
    }

    void conjugate(size_t shift, size_t nVec) {
      for(size_t j = 0; j < nVec; j++) {
	#pragma omp parallel for
        for (auto it = data[j + shift].begin(); it != data[j + shift].end(); ++it){
	  it->second = SmartConj(it->second);
	}
      }
    }

    void scale(MatsT scalar, size_t shift, size_t nVec) {
      if(scalar == MatsT(0.)) {
	setZero(shift, nVec);
      }
      else if(scalar != MatsT(1.)) {
        for(size_t j = 0; j < nVec; j++) {
	  #pragma omp parallel for
	  for (auto it = data[j + shift].begin(); it != data[j + shift].end(); ++it){
	    it->second *= scalar;
	  }
        }
      }
    }

    //add another matrix * beta to this
    void add(const LLSparseMatrix& B, MatsT beta, size_t shiftA, size_t shiftB, size_t nVec) {

      #pragma omp parallel for default(shared)
      for(size_t j = 0; j < nVec; j++) {
/*
	for(auto it = B.data[shiftB + j].cbegin(); it != B.data[shiftB + j].cend(); it++) {
	  auto itA = std::lower_bound(data[shiftA + j].begin(), data[shiftA + j].end(), *it, [](const auto& a, const auto& b){return a.first < b.first;});
	  
	  if(itA == data[shiftA + j].end() or itA->first != it->first) {
	    data[shiftA + j].insert(itA, std::make_pair(it->first, beta * it->second));
	  }
	  else {
	    itA->second += beta * it->second;
	  }
	}

	*/
	if(B.data[shiftB + j].size() == 0) { continue; }

	size_t endA = data[shiftA + j].size();	
	data[shiftA + j].resize(endA + B.data[shiftB + j].size());
        std::copy(CQ_EXEC_PAR B.data[shiftB + j].cbegin(), B.data[shiftB + j].cend(), data[shiftA + j].begin() + endA);
	
	#pragma omp parallel for default(shared)
	for(auto it =  data[shiftA + j].begin() + endA; it != data[shiftA + j].end(); it++) {
	  it->second *= beta;	
	}

	std::inplace_merge(CQ_EXEC_PAR data[shiftA + j].begin(), data[shiftA + j].begin() + endA, data[shiftA + j].end(),
	  [&] (const auto& i, const auto& j) {
            return i.first < j.first;
          }
        );

	#pragma omp parallel for default(shared)	
	for(size_t i = 0; i < data[shiftA + j].size() - 1; ++i) {
          if(data[shiftA + j][i].first == data[shiftA + j][i + 1].first) {
            data[shiftA + j][i].second += data[shiftA + j][i + 1].second;
            data[shiftA + j][i + 1].second = MatsT(0.);
          }
        }
	
	auto newEnd = std::remove_if(CQ_EXEC_PAR  data[shiftA + j].begin(),  data[shiftA + j].end(), [](const auto& x){ return x.second == MatsT(0.); });
        data[shiftA + j].erase(newEnd, data[shiftA + j].end());
	
      }
    }

    void kAdd(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
              MatsT alpha, MatsT const *B, int64_t ldb,
              LLSparseMatrix<MatsT>& C, size_t shiftC) const {

      if(alpha == MatsT(0.))
        return;

      for(size_t Ccol = 0; Ccol < n; Ccol++) {
        auto& res = C.data[shiftC + Ccol];
          for(size_t Acol = 0; Acol < k; Acol++) {
            switch(transB) {
              case blas::Op::NoTrans:
              {
                if(B[Acol + ldb * Ccol] != MatsT(0.)) {
                  res.resize(res.size() + data[shiftA + Acol].size());
                  std::transform(CQ_EXEC_PAR data[shiftA + Acol].cbegin(),  data[shiftA + Acol].cend(), res.end() - data[shiftA + Acol].size(),
                    [&](auto& x){return std::make_pair(x.first, alpha * B[Acol + ldb * Ccol] * (x.second));}
                  );
                }
                break;
              }
              case blas::Op::Trans:
              {
                if(B[Ccol + ldb * Acol] != MatsT(0.)) {
                  res.resize(res.size() + data[shiftA + Acol].size());
                  std::transform(CQ_EXEC_PAR data[shiftA + Acol].cbegin(),  data[shiftA + Acol].cend(), res.end() - data[shiftA + Acol].size(),
                    [&](auto& x){return std::make_pair(x.first, alpha * B[Ccol + ldb * Acol] * (x.second));}
                  );
                }
                break;
              }
              case blas::Op::ConjTrans:
              {
                if(B[Ccol + ldb * Acol] != MatsT(0.)) {
                  res.resize(res.size() + data[shiftA + Acol].size());
                  std::transform(CQ_EXEC_PAR data[shiftA + Acol].cbegin(),  data[shiftA + Acol].cend(), res.end() - data[shiftA + Acol].size(),
                    [&](auto& x){return std::make_pair(x.first, alpha * SmartConj(B[Ccol + ldb * Acol]) * (x.second));}
                  );
                }
              break;
            }
            default:
              CErr("multiply_matrix with this blas::Op is not implemented yet");
          }
        }

        boost::sort::block_indirect_sort(res.begin(), res.end(),
          [&] (const auto& i, const auto& j) {
          return i.first < j.first;
        }, GetNumThreads());

        //collapse duplicates
        size_t resIdx = 0;
        while(resIdx < res.size()) {
          size_t uniqueI = resIdx;
          resIdx++;

          while(resIdx < res.size() and res[resIdx].first == res[uniqueI].first) {
            res[uniqueI].second += res[resIdx].second;
            res[resIdx].second = MatsT(0.);
            resIdx++;
          }
        }

        auto newEnd = std::remove_if(CQ_EXEC_PAR res.begin(), res.end(), [](const auto& x){ return x.second == MatsT(0.); });
        res.erase(newEnd, res.end());
      }
    }


    MatsT dot(const LLSparseMatrix& B, size_t shiftA, size_t shiftB) const {
      
	auto itA = data[shiftA].begin();
        auto itB = B.data[shiftB].begin();
        MatsT dot = MatsT(0.);
 
	size_t iA = 0;
	size_t iB = 0;

	while (iA < data[shiftA].size() and iB < B.data[shiftB].size()) {
	  size_t deltaA = 1, deltaB = 1;

	  if(data[shiftA][iA].first < B.data[shiftB][iB].first) {
	    while(iA < data[shiftA].size() and data[shiftA][iA].first < B.data[shiftB][iB].first) {
	      iA += deltaA;
	      deltaA *= 2;
	    }

	    auto lowerBoundA = std::lower_bound(data[shiftA].begin() + (iA - deltaA / 2), data[shiftA].begin() + std::min(iA, data[shiftA].size()), B.data[shiftB][iB], 
	      [](const auto& i, const auto& j){return i.first < j.first; });
	    iA = std::distance(data[shiftA].begin(), lowerBoundA);
	  }
	  else if(B.data[shiftB][iB].first < data[shiftA][iA].first) {
	    while(iB < B.data[shiftB].size() and B.data[shiftB][iB].first < data[shiftA][iA].first) {
	      iB += deltaB;
	      deltaB *= 2;
	    }

	    auto lowerBoundB = std::lower_bound(B.data[shiftB].begin() + (iB - deltaB / 2), B.data[shiftB].begin() + std::min(iB, B.data[shiftB].size()), data[shiftA][iA], 
	      [](const auto& i, const auto& j){return i.first < j.first; });
	    iB = std::distance(B.data[shiftB].begin(), lowerBoundB);
	  }

	  if(iA < data[shiftA].size() and iB < B.data[shiftB].size() and data[shiftA][iA].first == B.data[shiftB][iB].first) {
	    dot += SmartConj(data[shiftA][iA].second) * B.data[shiftB][iB].second;
	    iA++;
	    iB++;
	  }
	}

	return dot;
    }

    double normSquared(size_t shift, size_t nVec) const {
      double v = 0.0;
      for(size_t j = 0; j < nVec; j++) {
        #pragma omp parallel for reduction(+:v)
        for(auto it = data[shift + j].begin(); it != data[shift + j].end(); it++)
          v += std::norm(it->second);
      }

      return v;
    }

    double maxNormElement(size_t shift, size_t nVec) const {
      double v = 0.0;
      for(size_t j = 0; j < nVec; j++) {
	#pragma omp parallel for reduction(max:v)
	for(auto it = data[shift + j].begin(); it != data[shift + j].end(); it++)
	  v = std::max(v, std::abs(it->second));
      }

      return v;
    }
    
    void setFromHashSparseMatrix(HashSparseMatrix<MatsT>& B, size_t shiftA, size_t shiftB, size_t nVec) {
      for(size_t col = 0; col < nVec; col++) {
        data[shiftA + col].resize(B.data[shiftB + col].size());

      boost::unordered_flat_map<size_t,
          MatsT,
          boost::hash<const size_t>,
          std::equal_to<const size_t>,
          MyAlloc<std::pair<const size_t, MatsT>> > serialBCol = std::move(B.data[shiftB + col]);



        std::copy(CQ_EXEC_PAR serialBCol.cbegin(), serialBCol.cend(), data[shiftA + col].begin());
	boost::sort::block_indirect_sort(data[shiftA + col].begin(), data[shiftA + col].end(), 
          [&] (const auto& i, const auto& j) {
            return i.first < j.first;
        }, GetNumThreads());
      }
    }
};

}//namespace ChronusQ
#endif
