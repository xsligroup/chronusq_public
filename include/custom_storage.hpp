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

#include <new>
#include <utility>
#include <unordered_map>

#include <iostream>

namespace ChronusQ {

  class CQMemBackend {
  public:
    virtual void* malloc(size_t) = 0;
    virtual void free(void *) = 0;
    /**
      *  Return the maximum number of objects allocatable
      *
      *  \prama [obj_size]      Size of the object to allocate
      *  \prama [request_max]   Maximum number of objects to allocate
      *  \returns               Maximum number of objects that can be allocated
      *                         Return value <= request_max
      */
    virtual size_t max_avail_allocatable(size_t obj_size, size_t request_max) = 0;
    virtual ~CQMemBackend() = default;
  };

  class OSDirectMemManager : public CQMemBackend {

  public:

    OSDirectMemManager() = default;

    /**
     * \brief Allocates a contiguous block of memory for N values of the given
     * type.
     *
     * \param [in] N  Number of items of type T to allocate.
     * \return        Pointer to block large enough to hold N objects of type T
     */
    void* malloc(size_t N) override {
      void* trial_alloc = std::malloc(N);
      if ( !trial_alloc )
        throw std::bad_alloc();
      return trial_alloc;

    }; // malloc

    /**
     * \brief Frees a pointer previously allocated by OSMemManager
     *
     * \param [in] ptr  Pointer to block to free
     */
    void free(void * ptr) override {

      std::free(ptr);

    }; // free

    /**
      *  Return the maximum number of objects allocatable
      *
      *  \prama [obj_size]      Size of the object to allocate
      *  \prama [request_max]   Maximum number of objects to allocate
      *  \returns               Maximum number of objects that can be allocated
      *                         Return value <= request_max
      */
    size_t max_avail_allocatable(size_t obj_size, size_t request_max) override {

      // First try if request_max is allocatable
      void * trial_alloc;
      try {
        trial_alloc = malloc(request_max * obj_size);
        free(trial_alloc);
        return request_max;
      } catch (std::bad_alloc& ba) {}

      size_t mem_max = request_max;
      size_t mem_min = 0;
      size_t trial_mem = 0;

      while(true) {
        try {
          trial_mem = (mem_max-mem_min)/2 + mem_min;
          trial_alloc = malloc(trial_mem * obj_size);
          free(trial_alloc);
          mem_min = trial_mem;
        } catch (std::bad_alloc& ba) {
          mem_max = trial_mem - 1;
        }
        if (mem_max - mem_min == 1 or mem_max == mem_min) break;
      }

      return mem_min;
    }

  }; // OSMemManager

  class CustomMemManager : public CQMemBackend {

    std::vector<char> V_;     ///< Internal memory
    void * top;               ///< Pointer to the top of the free list
    void * origin = nullptr;  ///< Pointer to the historic origin of the free list
    void * bottom = nullptr;  ///< Pointer to the historic bottom of the allocated memory
    size_t align_size; 
      ///< Size to which added and allocated blocks are aligned
    std::unordered_map<void*,size_t> alloc_blocks;
      ///< Maps allocated pointers to their size

      
    /**
     *  \brief Ensures that the memory block (N) is divisible by
     *  the segregation block size (BlockSize)
     */
    inline void fixBlockNumber(size_t &N, size_t BlockSize) {
      if(N % BlockSize) {
#ifdef MEM_PRINT
        std::cerr << "Memory Block not Divisible by BlockSize ("
                    << BlockSize << " Bytes). Increasing allocation by " 
                    << BlockSize - (N % BlockSize) << " Bytes" << std::endl;

#endif
        N += BlockSize - (N % BlockSize);
      }
    };
    
    /**
     *  Allocates the memory block
     */
    void allocMem(size_t N, size_t BlockSize) {

      fixBlockNumber(N, BlockSize); // fix buffer length
      V_.resize(N);    // allocate the memory

      // segregate (uses functionality from boost::simple_segregated_storage
      add_ordered_block(&V_.front(),V_.size(),BlockSize);

#ifdef MEM_PRINT
      std::cerr << "Creating Memory Partition of " << N
                   << " bytes starting at " << (int*) &V_.front() 
                   << std::endl;
#endif
    };
    
    /**
     *  \brief Gets the pointer to the next element in the free list
     *
     *  \param [in]  ptr  Pointer to current element in free list
     *  \return           Pointer to next element in free list
     */
    static inline void * & get_next(void * ptr) {
      return *(static_cast<void **>(ptr));
    };
  
    /**
     *  \brief Gets the size of the current element of the free list
     *
     *  \param [in] ptr  Pointer to current element in free list
     *  \return          Size of current element of free list
     */
    static inline size_t & get_size(void * ptr) {
      return reinterpret_cast<size_t&>(*(static_cast<void **>(ptr) + 1));
    };
  
    /**
     *  \brief Throws an exception if the pointer is not aligned to align_size
     *
     *  \param [in] ptr  Pointer to check
     */
    inline void check_alignment(void * ptr) {
      if (reinterpret_cast<size_t>(ptr) % align_size != 0)
        throw "Bad alignment of pointer in MemManager";
    };
  

    /**
     *  \brief Gets the pointer in the free list before the passed pointer
     *
     *  Scans the free list, searching for a pointer larger than the passed
     *  pointer. When found, the previous pointer is returned. Returns nullptr
     *  if the passed pointer is smaller than top.
     *
     *  \param [in]  ptr  Pointer to compare
     *  \return           Pointer to element before the passed pointer 
     */
    void * get_prev(void * ptr) {
  
      // Handle edge case that all blocks are after ptr
      if ( top > ptr )
        return nullptr;
  
      void * ret = top;
  
      while ( get_next(ret) ) {
        if ( get_next(ret) > ptr ) 
          break;
        ret = get_next(ret);
      }
  
      return ret;
    };
  
  
    public:
  
    /**
     * \brief Constructor
     *
     * Constructs a CustomMemManager object and sets the align_size. If no
     * align_size is specified, it defaults to sizeof(size_t).
     *
     * \param [in] align_size  The alignment size.
     */
    CustomMemManager(size_t N = 0, size_t BlockSize = 2048,
                     size_t align_size_ = sizeof(size_t)) :
        top(nullptr), align_size(align_size_) {

      if( N and BlockSize ) allocMem(N, BlockSize);
    };
  
    /**
     * \brief Adds the block of memory to the free list in an ordered manner
     *
     * Searches the free list for the appropriate location to store this block
     * to maintain order. Merges blocks together if they are contiguous. Throws
     * if the block is smaller than the header or misaligned. Does not restore
     * order if free list is not ordered.
     *
     * \param [in] block       Pointer to block to add
     * \param [in] block_size  Size of block to add (in bytes)
     */
    void add_block(void * block, size_t block_size) {
      if ( block_size == 0 )
        return;
      if ( block_size < (sizeof(void*) + sizeof(size_t)) ) {
        throw "Blocks need to be large enough to hold the linked list header";
      };
      check_alignment(block);
    
      // Handle edge case with empty free list
      if ( !top ) {
        get_next(block) = top;
        get_size(block) = block_size;
        top = block;
        return;
      }
    
      void * prev = get_prev(block);
      void * next = prev ? get_next(prev) : top;
    
      // Append or link next to block
      if ( static_cast<char*>(block) + block_size == next ) {
        get_next(block) = get_next(next);
        get_size(block) = block_size + get_size(next);
      }
      else {
        get_next(block) = next;
        get_size(block) = block_size;
      }
    
      // Append or link block to prev
      if ( !prev ) {
        top = block;
      }
      else if ( static_cast<char*>(prev) + get_size(prev) == block ) {
        get_next(prev) = get_next(block);
        get_size(prev) += get_size(block);
      }
      else {
        get_next(prev) = block;
      }
    
    };
  
    /**
     * \brief Adds the block of memory to the top of the free list
     *
     * Pushes the block to the top of the free list. Throws if the block is
     * smaller than the header or misaligned. Can cause fragmentation of the
     * free list, and does not merge contiguous blocks.
     *
     * \param [in] block       Pointer to block to add
     * \param [in] block_size  Size of block to add (in bytes)
     */
    void add_block_fast(void * block, size_t block_size) {
      if ( block_size == 0 )
        return;
      if ( block_size < sizeof(void*) + sizeof(size_t) )
        throw "Blocks need to be large enough to hold the linked list header";
      check_alignment(block);
    
      get_next(block) = top;
      get_size(block) = block_size;
      top = block;
    };
  
  
    /**
     * \brief Allocates a contiguous block of memory for N values of the given
     * type.
     *
     * Searches the free list for a block with the size required to hold N
     * objects of type T. Throws bad_alloc if the free list doesn't contain a
     * block large enough. 
     *
     * \param [in] N  Number of items of type T to allocate.
     * \return        Pointer to block large enough to hold N objects of type T
     */
    void* malloc(size_t N) override {

      // Add padding and make sure block is large enough for the header
      // (for freeing later)
      size_t mod_size = N % align_size;
      size_t req_size = mod_size == 0 ?
                        N : N + align_size - mod_size;
      size_t min_size = sizeof(void*) + sizeof(size_t);
      req_size = req_size > min_size ? req_size : min_size;
    
      // Throw if no free list
      if ( !top )
        throw std::bad_alloc();
    
      // Search for block large enough
      void* ptr(top);
      void* prev(nullptr);
      size_t block_size(get_size(ptr));
      while ( ptr && block_size < req_size ) {
        prev = ptr;
        ptr = get_next(ptr);
        block_size = ptr ? get_size(ptr) : 0;
      };
    
      // Throw if no block large enough
      if (!ptr)
        throw std::bad_alloc();

      void * new_free = static_cast<void*>(static_cast<char*>(ptr) + req_size);
      bottom = std::max(bottom, new_free);
      if (req_size == block_size) {
        // Link prev to next if the found block is exactly large enough
        if ( prev ){
          get_next(prev) = get_next(ptr);
        }
        else if (ptr) {
          top = get_next(ptr);
        }
        else {
          top = nullptr;
        }
      }
      else {
        // Link prev to remaining memory in current block
        get_next(new_free) = get_next(ptr);
        get_size(new_free) = get_size(ptr) - req_size;
        if (top == ptr) {
          top = new_free;
        }
        else {
          get_next(prev) = new_free;
        }
      }
    
      // Record that this has been allocated
      alloc_blocks[ptr] = req_size;
    
      return ptr;
    }; 
  
  
    /**
     * \brief Frees a pointer previously allocated by CustomMemManager
     * (ordered)
     *
     * A wrapper for add_block that checks if ptr was previously allocated.
     * Maintains order on an ordered free list.
     *
     * \param [in] ptr  Pointer to block to free
     */
    void free(void * ptr) override {
      auto entry = alloc_blocks.find(ptr);
    
      if (entry == alloc_blocks.end())
        throw "Pointer was not allocated by MemManager";
    
      add_block(ptr, entry->second);
    
      alloc_blocks.erase(entry);
    };
  
    /**
     * \brief Frees a pointer previously allocated by CustomMemManager
     * (unordered)
     *
     * A wrapper for add_block_fast that checks if ptr was previously allocated
     * Can cause fragmentation of the free list.
     *
     * \param [in] ptr  Pointer to block to free
     */
    void free_fast(void * ptr) {
      auto entry = alloc_blocks.find(ptr);
    
      if (entry == alloc_blocks.end())
        throw "Pointer was not allocated by MemManager";
    
      add_block_fast(ptr, entry->second);
    
      alloc_blocks.erase(entry);
    };
  
    /**
     * \brief Defragments free list and restores order
     *
     * Sorts and remerges the free list so it is ordered and contains no
     * adjacent blocks
     */
    void defrag() {
      void * ptr(top);
      void * next;
      
      // Delete free list
      top = nullptr;
    
      // Traverse previous free list, adding blocks in an ordered manner
      while ( ptr ) {
        next = get_next(ptr);
        add_block(ptr, get_size(ptr));
        ptr = next;
      }
    };
  
  #ifdef MEMMANAGER_DEBUG
    void print_free() {
      void * ptr(top);
  
      std::cout << "                 MemManager Free List  \n";
      std::cout << "---------------------------------------------------------\n";
      while (ptr) {
        std::cout << ptr << " | " << get_size(ptr) << '\n';
        ptr = get_next(ptr);
      }
      std::cout << "---------------------------------------------------------\n";
    };
  #endif

    /**
     *  Mimic function for boost::simple_segregated_storage::add_ordered_block
     *
     * \param [in] block  Pointer to block to add
     * \param [in] nsz    Size of block to add (in bytes)
     * \param [in] dummy  Dummy to match call signature
     */
    void add_ordered_block(void * const block,
      const size_t nsz, const size_t dummy) {
      origin = block;
      add_block(block, nsz);
    };

    /**
     *  Return the span of the allocated memory
     */
    size_t alloc_span() const {
      return static_cast<char*>(bottom) - static_cast<char*>(origin);
    };


    size_t getTotalAllocation() const {
      return V_.size();
    };



    /**
      *  Return the maximum number of objects allocatable
      *
      *  \prama [obj_size]      Size of the object to allocate
      *  \prama [request_max]   Maximum number of objects to allocate
      *  \returns               Maximum number of objects that can be allocated
      *                         Return value <= request_max
      */
    size_t max_avail_allocatable(size_t obj_size, size_t request_max) override {

      size_t avail_max = 0;
      void * ptr(top);
      while (ptr) {
        avail_max = std::max(get_size(ptr), avail_max);
        ptr = get_next(ptr);
      }

      return std::min(avail_max/obj_size, request_max);
    }
  
  }; // MemManager

}; // namespace ChronusQ

