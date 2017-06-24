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

#include <chronusq_sys.hpp>
#include <custom_storage.hpp>

#ifdef _OPENMP
#define MEM_IN_OMP_WARNING(func) \
    do { \
        if (omp_in_parallel()) { \
            std::cout << "Warning: " #func \
                      << " called from within parallel region. " \
                      << "This may lead to deadlock." << std::endl; \
        } \
    } while (false)
#else
#define MEM_IN_OMP_WARNING(func) do {} while (false)
#endif

//#define MEM_PRINT

namespace ChronusQ {

  enum class CQMemBackendType {
    PREALLOCATED,
    OS_DIRECT,
  };

  class CQMemManager {
  private:

    std::unique_ptr<CQMemBackend> mem_backend; ///< Memory backend
    size_t BlockSize_ = 2048; ///< Segregation block size
    size_t NAlloc_ = 0;       ///< Number of blocks currently allocated
    size_t NAllocHigh_ = 0;   ///< High-water mark of allocated blocks

    std::unordered_map<void*,std::pair<size_t,size_t>> AllocatedBlocks_;
      ///< Map from block pointer to the size of the block

    // Private default constructor
    CQMemManager() = default;

  public:

    // Disable copy and move construction and assignment
    CQMemManager(const CQMemManager &)            = delete;
    CQMemManager(CQMemManager &&)                 = delete;
    CQMemManager& operator=(const CQMemManager &) = delete;
    CQMemManager& operator=(CQMemManager &&)      = delete;

    static CQMemManager& get() {
      static CQMemManager instance;
      return instance;
    }

    /**
     *  \brief initializer.
     *
     *  Initialize a CQMemManager object with OS-direct memmanager or
     *  with preallocated chuck of memory. The latter has optionally defaulted values.
     *  If no paramemters are set, it will default to no memory (0 bytes)
     *  allocated and a 256 byte segregation length. If the requested
     *  memory and segregation legnth are not 0 , it will allocate and
     *  segregate the block
     *
     *  \param [in] type       Memory backend type
     *  \param [in] N          Total memory (in bytes) to be allocated
     *  \param [in] BlockSize  Segregation block size
     *
     */ 
    void initialize(CQMemBackendType type, size_t N = 0, size_t BlockSize = 2048) {
      MEM_IN_OMP_WARNING(CQMemManager::initialize);

      BlockSize_ = BlockSize;
      NAlloc_ = 0;
      NAllocHigh_ = 0;
      switch (type) {
        case CQMemBackendType::PREALLOCATED:
          mem_backend = std::make_unique<CustomMemManager>(N, BlockSize);
          break;
        case CQMemBackendType::OS_DIRECT:
          mem_backend = std::make_unique<OSDirectMemManager>();
          break;
        default:
          throw std::runtime_error("Unknown memory backend type");
      }
    };


     /**
      *  \brief Allocates a contiguous block of memory of a specified
      *  type.
      *
      *  Allocates the correct number of blocks from the segregated
      *  storage given a type and a number to allocate
      *
      *  \param [in] n  Number of items of type T to allocate
      *  \return        Pointer to contiguous memory block that contains
      *                 n items of type T
      */ 
     template <typename T>
     T* malloc(size_t n) {
       MEM_IN_OMP_WARNING(CQMemManager::malloc);
       // Determine the number of blocks to allocate
       size_t nBlocks = ( (n-1) * sizeof(T) ) / BlockSize_ + 1;
      
       #ifdef MEM_PRINT
         std::cerr << "Allocating " << n << " words of " << typeid(T).name()
                   << " data (" << nBlocks << " blocks): ";
       #endif

       // Get a pointer from boost::simple_segregated_storage
       void * ptr = mem_backend->malloc(nBlocks * BlockSize_);

       // Throw an error if boost returned a NULL pointer (many
       // possible causes)
       if(ptr == NULL) {
         std::bad_alloc excp;
         throw excp;
       }

       //// Update the number of allocated blocks
       NAlloc_ += nBlocks;
       NAllocHigh_ = std::max(NAlloc_,NAllocHigh_);
       
       #ifdef MEM_PRINT
         std::cerr << "  PTR = " << ptr << std::endl;
       #endif

       assert( AllocatedBlocks_.find(static_cast<void*>(ptr)) == AllocatedBlocks_.end() );

       // Keep a record of the block
       AllocatedBlocks_[ptr] = { n * sizeof(T), nBlocks }; 

       return static_cast<T*>(ptr); // Return the pointer
     }; // CQMemManager::malloc

    // malloc with default initialization
    template <typename T>
    T* calloc(size_t n) {
      T* ptr = malloc<T>(n);
      std::fill_n(ptr, n, T());
      return ptr;
    }; // CQMemManager::calloc

     
     /**
      *  Frees a contiguous memory block given a pointer previously
      *  returned by CQMemManager::malloc.
      *
      *  \warning Will throw an error if the pointer is not currently
      *  in the list of allocated pointers
      *
      *  \param [in] ptr Pointer to free 
      */ 
     template <typename T>
     void free( T* &ptr ) {
       MEM_IN_OMP_WARNING(CQMemManager::free);

       // Attempt to find the pointer in the list of 
       // allocated blocks
       auto it = AllocatedBlocks_.find(static_cast<void*>(ptr));

       // Kill the job if the pointer is not in the list of allocated
       // blocks
       assert( it != AllocatedBlocks_.end() );

       #ifdef MEM_PRINT
         std::cerr << "Freeing " << it->second.second 
                   << " blocks of data starting from " 
                   << static_cast<void*>(ptr) << std::endl;
       #endif

       NAlloc_ -= it->second.second; // deduct block size from allocated memory
  
       // deallocate the memory in an ordered fashion
       mem_backend->free(ptr);

       // Remove pointer from allocated list
       AllocatedBlocks_.erase(it);
                
       ptr = NULL; // NULL out the pointer
     }; // CQMemManager::free

     /**
      *  Parameter pack of CQMemManager::free. Allows for subsequent
      *  deallocation of an arbitrary number of memory blocks.
      *
      *  z.B.
      *
      *  free(X,Y,Z); 
      *
      *  is equivalant to
      *
      *  free(X);
      *  free(Y);
      *  free(Z);
      */ 
     template <typename T, typename... Targs>
     void free( T* &ptr, Targs... args) {
       // Free the first pointer then recurse down
       free(ptr); free(args...);
     }; // CQMemManager::free (parameter pack)


     /**
      *  Returns the size of an allocated memory block
      *
      *  \warning  Dies if pointer is not in the memory block
      *
      *  \param [in] ptr Pointer of interest
      *  \returns        Size of the allocated block in terms of type T
      */ 
     template <typename T>
     size_t getSize(T* ptr) const {
       // Attempt to find the pointer in the list of 
       // allocated blocks
       auto it = AllocatedBlocks_.find(static_cast<void*>(ptr));

       // Kill the job if the pointer is not in the list of allocated
       // blocks
       assert( it != AllocatedBlocks_.end() );

       return std::floor(it->second.first / sizeof(T));
     }; // CQMemManager::getSize


     /**
      *  Return the maximum number of objects allocatable
      *
      *  \prama [n_elem_each]   How many T type elements are in each object
      *  \prama [request_max]   Maximum number of objects to allocate
      *  \returns               Maximum number of objects that can be allocated
      *                         Return value <= request_max
      */
     template <typename T>
     size_t max_avail_allocatable(size_t n_elem_each, size_t request_max) {
       return mem_backend->max_avail_allocatable(sizeof(T) * n_elem_each, request_max);
     }




     /**
      *  Prints the CQMemManager allocation table to a specified output
      *  device.
      *
      *  \param [in] out Output device to print the table
      */ 
     void printAllocTable(std::ostream &out) const {
       out << "Allocation Table (unordered):\n\n";
       out << std::left;
       out << std::setw(15) << "Pointer" << std::setw(15) << "Size (Bytes)" 
           << std::endl;
       for( auto &block : AllocatedBlocks_ )
         out << std::setw(15) << static_cast<void*>(block.first) 
             << std::setw(15) <<  BlockSize_*block.second.second 
             << std::endl;

     }; // CQMemManager::printAllocTable


    /**
     *  Outputs a summary of the memory allocation in its current state
     *
     *  \param [in/out] out Output device.
     *  \param [in]     mem CQMemManager object of interest
     */ 
    friend inline std::ostream& operator<<(std::ostream &out , 
      const CQMemManager &mem) {

      out << "Memory Allocation Summary:" << std::endl << std::left;
      
      auto outputFunc = [&](std::string str, double size, size_t nBlocks) {
         out << std::setw(30) << str << std::setw(30) << std::fixed << size << "  B / " << std::endl; 
         out << std::setw(30) << ""  << std::setw(30) << std::fixed << size / 1e3 << " KB /" << std::endl; 
         out << std::setw(30) << ""  << std::setw(30) << std::fixed << size / 1e6 << " MB /" << std::endl; 
         out << std::setw(30) << ""  << std::setw(30) << std::fixed << size / 1e9 << " GB " << std::endl;
         out << std::setw(30) << ""  << "(" << nBlocks << " Blocks) ";
      };

      out << std::setw(30) << " - Block Size: ";
      out << std::setw(10) << mem.BlockSize_ << " B" << std::endl; 
      out << std::endl;

      if (const CustomMemManager *cmm = dynamic_cast<const CustomMemManager*>(mem.mem_backend.get())) {
        size_t N = cmm->getTotalAllocation();

        outputFunc(" - Total Memory Allocated:", N, N / mem.BlockSize_);
        out << std::endl;

        outputFunc(" - Free Memory:", N - mem.NAlloc_ * mem.BlockSize_,
                   N / mem.BlockSize_ - mem.NAlloc_);
        out << std::endl;
      }

      outputFunc(" - Reserved Memory:", mem.NAlloc_ * mem.BlockSize_, mem.NAlloc_);
      out << std::endl;

      if( mem.NAlloc_ ) {
        out << std::endl << std::endl;;
        
        mem.printAllocTable(out);
      }

      return out; // Return the ouput device

    }; // CQMemManager operator<<

    /**
     *  Prints the CQMemManager high-water mark to a specified output
     *  device.
     *
     *  \param [in] out Output device to print the table
     */
    void printHighWaterMark(std::ostream &out) const {
      out << std::endl << "MemManager high-water mark: "
          << std::fixed << std::setprecision(3);
      if (CustomMemManager *cmm = dynamic_cast<CustomMemManager*>(mem_backend.get()))
        out << cmm->alloc_span() / 1e9;
      else
        out << NAllocHigh_ * BlockSize_ / 1e9;
      out << " GB." << std::endl;

    }; // CQMemManager::printAllocTable

  }; // class CQMemManager

}; // namespace ChronusQ

