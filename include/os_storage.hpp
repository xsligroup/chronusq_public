/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

#include <iostream>

namespace ChronusQ {

  class OSMemManager {

    public:

    OSMemManager() = default;

    /**
     * \brief Allocates a contiguous block of memory for N values of the given
     * type.
     *
     * \param [in] N  Number of items of type T to allocate.
     * \return        Pointer to block large enough to hold N objects of type T
     */
    template <typename T>
    T* malloc(size_t N) {

      return static_cast<T*>(std::malloc(N * sizeof(T)));

    }; // malloc

    /**
     * \brief Frees a pointer previously allocated by OSMemManager
     *
     * \param [in] ptr  Pointer to block to free
     */
    void free(void * ptr) {

      std::free(ptr);

    }; // free

    //
    // Mimic functions for boost::simple_segregated_storage
    //

    /**
     *  Mimic function for boost::simple_segregated_storage::add_ordered_block
     *
     * \param [in] block  Pointer to block to add
     * \param [in] nsz    Size of block to add (in bytes)
     * \param [in] dummy  Dummy to match call signature
     */
    void add_ordered_block(void * const, const size_t, const size_t) {};

    /**
     *  Mimic function for boost::simple_segregated_storage::malloc_n
     *
     * \param [in] n           Number of blocks to allocate
     * \param [in] block_size  Size of block to allocate (in bytes)
     * \return                 Pointer to memory of size n*block_size
     */
    void * malloc_n(size_t n, size_t block_size) {
      return static_cast<void *>(malloc<char>(n * block_size));
    };

    /**
     *  Mimic function for boost::simple_segregated_storage::ordered_free_n
     *
     * \param [in] chunks  Pointer to block to free
     * \param [in] dummy1  Dummy to match call signature
     * \param [in] dummy2  Dummy to match call signature
     */
    void ordered_free_n(void * const chunks, const size_t, const size_t) {
      free(chunks);
    };

    /**
     *  Return the span of the allocated memory
     */
    size_t alloc_span() const {
      return 0;
    };


  }; // OSMemManager

}; // namespace ChronusQ
