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

//#define MEM_PRINT
#define CHRONUSQ_CUSTOM_BACKEND

#ifdef CHRONUSQ_CUSTOM_BACKEND
#include <custom_storage.hpp>
#else
#include <boost/pool/simple_segregated_storage.hpp>
#endif


namespace ChronusQ {

#ifdef CHRONUSQ_CUSTOM_BACKEND
  typedef CustomMemManager mem_backend;
#else
  typedef boost::simple_segregated_storage<size_t> mem_backend;
#endif

  class CQMemManager : public mem_backend {

    size_t N_;            ///< Total bytes to be allocated
    size_t NAlloc_;       ///< Number of blocks currently allocated
    size_t BlockSize_;    ///< Segregation block size
    std::vector<char> V_; ///< Internal memort

    std::unordered_map<void*,std::pair<size_t,size_t>> AllocatedBlocks_;
      ///< Map from block pointer to the size of the block

    bool isAllocated_;

    /**
     *  \brief Ensures that the memory block (N_) is divisible by
     *  the segregation block size (BlockSize_)
     */ 
    inline void fixBlockNumber() {
      if(N_ % BlockSize_) {
        #ifdef MEM_PRINT
          std::cerr << "Memory Block not Divisible by BlockSize ("
                    << BlockSize_ << " Bytes). Increasing allocation by " 
                    << BlockSize_ - (N_ % BlockSize_) << " Bytes" << std::endl;
  
        #endif
        N_ += BlockSize_ - (N_ % BlockSize_);
      }
    };

    public:

    // Disable copy and move construction and assignment
    CQMemManager(const CQMemManager &)            = delete;
    CQMemManager(CQMemManager &&)                 = delete;
    CQMemManager& operator=(const CQMemManager &) = delete;
    CQMemManager& operator=(CQMemManager &&)      = delete;

    /**
     *  \brief Constructor.
     *
     *  Constructs a CQMemManager object with optionally defaulted values.
     *  If no paramemters are set, it will default to no memory (0 bytes)
     *  allocated and a 256 byte segregation length. If the requested
     *  memory and segregation legnth are not 0 , it will allocate and
     *  segregate the block
     *
     *  \param [in] N          Total memory (in bytes) to be allocated
     *  \param [in] BlockSize  Segregation block size
     *
     */ 
     CQMemManager(size_t N = 0, size_t BlockSize = 2048) :
       mem_backend(), N_(N), BlockSize_(BlockSize), isAllocated_(false),
       NAlloc_(0) {
       if( N_ and BlockSize_ ) allocMem();
     };

     /**
      *  Allocates the memory block
      */ 
     void allocMem() {
       assert(not isAllocated_);

       fixBlockNumber(); // fix buffer length
       V_.resize(N_);    // allocate the memory

       // segregate (uses functionality from boost::simple_segregated_storage
       this->add_ordered_block(&V_.front(),V_.size(),BlockSize_);

       isAllocated_ = true; // ensure no additional allocs can be made

       #ifdef MEM_PRINT
         std::cerr << "Creating Memory Partition of " << N_
                   << " bytes starting at " << (int*) &V_.front() 
                   << std::endl;
       #endif
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
       // Determine the number of blocks to allocate
       size_t nBlocks = ( (n-1) * sizeof(T) ) / BlockSize_ + 1;
      
       // Check to see if requested memory would overflow allocated
       // memory. Throw a bad alloc if so
       //if( (NAlloc_ + nBlocks) * BlockSize_ > N_ ) {
       //  std::bad_alloc excp;
       //  throw excp;
       //}

       //// Update the number of allocated blocks
       //NAlloc_ += nBlocks; 

       #ifdef MEM_PRINT
         std::cerr << "Allocating " << n << " words of " << typeid(T).name()
                   << " data (" << nBlocks << " blocks): ";
       #endif

       // Get a pointer from boost::simple_segregated_storage
       void * ptr = mem_backend::malloc_n(nBlocks,BlockSize_);

       // Throw an error if boost returned a NULL pointer (many
       // possible causes)
       if(ptr == NULL) {
         std::bad_alloc excp;
         throw excp;
       }

       #ifdef MEM_PRINT
         std::cerr << "  PTR = " << ptr << std::endl;
       #endif

       assert( AllocatedBlocks_.find(static_cast<void*>(ptr)) == AllocatedBlocks_.end() );

       // Keep a record of the block
       AllocatedBlocks_[ptr] = { n * sizeof(T), nBlocks }; 

       return static_cast<T*>(ptr); // Return the pointer
     }; // CQMemManager::malloc


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
       mem_backend::ordered_free_n(ptr,it->second.second,BlockSize_);

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
     size_t getSize(T* ptr) {
       // Attempt to find the pointer in the list of 
       // allocated blocks
       auto it = AllocatedBlocks_.find(static_cast<void*>(ptr));

       // Kill the job if the pointer is not in the list of allocated
       // blocks
       assert( it != AllocatedBlocks_.end() );

       return std::floor(it->second.first / sizeof(T));
     }; // CQMemManager::getSize


     /**
      *  Return the maximum allocatable size
      *
      */
     template <typename T>
     size_t max_avail_allocatable(size_t size = 1) {

       size_t mem_max = N_/(sizeof(T) * size);
       size_t mem_min = 0;
       size_t trial_mem = 0;
       T * trial_alloc;

       while(true) {
         try {
           trial_mem = (mem_max-mem_min)/2 + mem_min;
           trial_alloc = CQMemManager::template malloc<T>(trial_mem * size);
           free(trial_alloc);
           mem_min = trial_mem;
         } catch (std::bad_alloc& ba) {
           mem_max = trial_mem - 1;
         }
         if (mem_max - mem_min == 1 or mem_max == mem_min) break;
       }

       return mem_min;
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

      out << std::setw(30) << "Total Memory Allocated: ";
      out << std::setw(15) << mem.N_ << " B / "; 
      out << std::setw(8) << std::fixed << mem.N_ / 1e3 << " kB / "; 
      out << std::setw(8) << std::fixed << mem.N_ / 1e6 << " MB / "; 
      out << std::setw(8) << std::fixed << mem.N_ / 1e9 << " GB"; 

      out << std::endl;

      out << std::setw(30) << "Block Size: ";
      out << std::setw(15) << mem.BlockSize_ << " B"; 
      out << " (" << mem.N_ / mem.BlockSize_ << " Blocks)";

      out << std::endl;

      out << std::setw(30) << "Reserved Memory: ";
      out << std::setw(15) << mem.NAlloc_*mem.BlockSize_ << " B"; 
      out << " (" << mem.NAlloc_ << " Blocks)";

      out << std::endl;

      out << std::setw(30) << "Free Memory: ";
      out << std::setw(15) << mem.N_ - mem.NAlloc_*mem.BlockSize_ << " B"; 
      out << " (" << mem.N_ / mem.BlockSize_ - mem.NAlloc_ << " Blocks)";

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
          << std::fixed << std::setprecision(3)
          << alloc_span() / 1e9
          << " GB." << std::endl;

    }; // CQMemManager::printAllocTable

  }; // class CQMemManager

}; // namespace ChronusQ

