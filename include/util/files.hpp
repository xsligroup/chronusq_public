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
#include <cxxapi/input.hpp>
#include <cerr.hpp>
#include <matrix/ndarray.hpp>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataType.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include <variant>

namespace ChronusQ {

// #define DEBUG_HDF5

namespace cqmatrix {
  template <typename MatsT>
  class PauliSpinorMatrices;
}

  class SafeFile {
  
    std::string fName_;
    bool        exists_;
  
    public:

      // Defaulted ctors
      SafeFile(const SafeFile &) = default;

      // String ctor
      SafeFile( const std::string &fName = "", 
        bool exists = false ) :
        fName_(fName), exists_(exists) { }
  

      // Member functions

      inline bool exists() const { return exists_; }
      inline std::string fName() const{ return fName_; }
      inline void setFile(const std::string &name) { fName_ = name; }

      inline void createFile() {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        exists_ = true;
      }; 

      inline void createGroup(const std::string &group) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        if (file.exist(group)){
          file.getGroup(group);
        } else { 
          file.createGroup(group); 
        }
      }

      bool exists(const std::string &dataSet) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        return file.exist(dataSet);
      };

      template <typename T>
      inline HighFive::DataSet createDataSet(const std::string &dataSet,
              const std::vector<size_t> &dims) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        if (file.exist(dataSet)){
          return file.getDataSet(dataSet);
        } else {
          // std::vector<size_t> max_dims;
          // std::vector<hsize_t> chunk_dims;
          // hsize_t DEFAULT_CHUNK = 64;
          // for (auto dim : dims){
          //     max_dims.push_back(std::numeric_limits<size_t>::max());
          //     chunk_dims.push_back(DEFAULT_CHUNK);
          // }
          // // todo 2d resizable datasets
          // auto dataspace = HighFive::DataSpace(dims, dims);
          // HighFive::DataSetCreateProps props;
          // props.add(HighFive::Chunking(chunk_dims));
          //TODO Should we do any chunking/compression?
          auto dataspace = HighFive::DataSpace(dims);
          HighFive::DataSet dataset =  file.createDataSet<T>(dataSet, dataspace);
          //HighFive::DataSet dataset =  file.createDataSet(dataSet, dataspace, HighFive::create_datatype<T>(), props);
          return dataset;
        }
      };

      /**
       *
       @brief Safely find or recreate a dataset with given name and dimensions.
       *        If the dataset exists and has the correct dimensions, it is returned.
       *        If it exists but has different dimensions, it is deleted and recreated.
       *        If it does not exist, it is created.
       *
       * @tparam T The data type of the dataset.
       * @param name The name of the dataset.
       * @param dims The desired dimensions of the dataset.
       * @return HighFive::DataSet The found or newly created dataset.
       */
      template <typename T>
      HighFive::DataSet safeFindOrRecreate(
          const std::string& name,
          const std::vector<size_t>& dims
      ) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);

        // Dataset does not exist → create
        if (!file.exist(name)) {
          return createDataSet<T>(name, dims);
        }

        HighFive::DataSet dset = file.getDataSet(name);

        if (dset.getDimensions() == dims) {
          return dset;
        }

        // Recreate path
        file.unlink(name);
        return createDataSet<T>(name, dims);
      }

      template <typename T>
      void readData(const std::string &dataSet, T* data) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        if (file.exist(dataSet)){
          auto datasetObj = file.getDataSet(dataSet);
          datasetObj.read<T>(data);
        } else {
          CErr("HDF5 CQ IO Issue. Attempting to read from nonexistent dataset.");
        }
      };

      template <typename T>
      void readData(const std::string &dataSet,
                    cqmatrix::PauliSpinorMatrices<T> &data) {
        readData(dataSet + "_SCALAR", data.S().pointer());
        if (data.hasZ())
          readData(dataSet + "_MZ", data.Z().pointer());
        if (data.hasXY()) {
          readData(dataSet + "_MY", data.Y().pointer());
          readData(dataSet + "_MX", data.X().pointer());
        }
      };

      std::variant< std::shared_ptr<cqmatrix::NDArray<double>>,
                    std::shared_ptr<cqmatrix::NDArray<dcomplex>>,
                    std::shared_ptr<cqmatrix::NDArray<int>> >
      readNDArray(const std::string &dataSet) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        H5T_class_t type;
        if (file.exist(dataSet)){
          auto datasetObj = file.getDataSet(dataSet);
          hid_t type_id = H5Dget_type(datasetObj.getId());
          type = H5Tget_class(type_id);
        } else {
          CErr("HDF5 CQ IO Issue. Attempting to read from nonexistent dataset.");
        }
        std::vector<size_t> dims = getDims(dataSet);
        if (type == H5T_FLOAT) {
          auto data = std::make_shared<cqmatrix::NDArray<double>>(dims);
          readData(dataSet, data->pointer());
          return data;
        } else if (type == H5T_COMPOUND) {
          auto data = std::make_shared<cqmatrix::NDArray<dcomplex>>(dims);
          readData(dataSet, data->pointer());
          return data;
        } else if (type == H5T_INTEGER) {
          auto data = std::make_shared<cqmatrix::NDArray<int>>(dims);
          readData(dataSet, data->pointer());
          return data;
        } else {
          ChronusQ::CErr("Unsupported data type for NDArray");
        }
      };

      template <typename T>
      void partialReadData(const std::string &dataSet, T* data,
                           const std::vector<size_t> &start, const std::vector<size_t> &dims,
                           const std::vector<size_t> &memStart = {},
                           const std::vector<size_t> &memDims = {} ) {
          HighFive::File file(fName_, HighFive::File::OpenOrCreate);
          if (file.exist(dataSet)){
            HighFive::DataSet dataset = file.getDataSet(dataSet);
            dataset.select(start, dims).read(data);
          } else { 
              CErr("HDF5 CQ IO Issue. Attempting to read from nonexistent dataset.");
          }
      };


      template <typename T>
      void partialWriteData(const std::string &dataSet, T* data,
          const std::vector<size_t> &start, const std::vector<size_t> &dims,
          const std::vector<size_t> &memStart = {},
          const std::vector<size_t> &memDims = {} ) {
          #ifdef DEBUG_HDF5
            std::cout << "In partialWriteData: " << std::endl;
            std::cout << "dataSet" << std::endl;
            std::cout << dataSet << std::endl;
            std::cout << "start" << std::endl;
            for (auto s : start){
               std::cout << s << ", ";
            }
            std::cout << std::endl;
            std::cout << "dims" << std::endl;
            for (auto s : dims){
               std::cout << s << ", ";
            }
            std::cout << std::endl;

            std::cout << "memStart" << std::endl;
            for (auto s : memStart){
               std::cout << s << ", ";
            }
            std::cout << std::endl;
            std::cout << "memDims" << std::endl;
            for (auto s : memDims){
               std::cout << s << ", ";
            }
            std::cout << std::endl;
          #endif
           // TODO last two parameters are not used but eliminate them later to provide compatability
          HighFive::File file(fName_, HighFive::File::OpenOrCreate);
          if (file.exist(dataSet)){
            HighFive::DataSet dataset = file.getDataSet(dataSet);
            //dataset.resize(dims);
            dataset.select(start, dims).write_raw(data);
          } else {
              CErr("HDF5 CQ IO Issue. Attempting to read from nonexistent dataset.");
          }
      };


      template <typename T>
      void safeWriteData(const std::string &dataSet, T* data,
        const std::vector<size_t> &dims) {
        #ifdef DEBUG_HDF5
          std::cout << "In safeWriteData: " << std::endl;
           std::cout << "dataSet" << std::endl;
            std::cout << dataSet << std::endl;
            std::cout << "dims" << std::endl;
            for (auto s : dims){
               std::cout << s << ", ";
            }
            std::cout << std::endl;
        #endif
          HighFive::DataSet dataset = safeFindOrRecreate<T>(dataSet, dims);
          dataset.write_raw(data);
      };

      template <typename T>
      void safeWriteData(const std::string &dataSet,
                         cqmatrix::PauliSpinorMatrices<T> &data) {
        size_t M = data.nRows();
        size_t N = data.nColumns();
        safeWriteData(dataSet + "_SCALAR", data.S().pointer(), {N,M}); // Note: N and M are swapped for column-major storage
        if (data.hasZ())
          safeWriteData(dataSet + "_MZ", data.Z().pointer(), {N,M});
        if (data.hasXY()) {
          safeWriteData(dataSet + "_MY", data.Y().pointer(), {N,M});
          safeWriteData(dataSet + "_MX", data.X().pointer(), {N,M});
        }
      };

      std::vector<size_t> getDims(const std::string &dataSet) {
        HighFive::File file(fName_, HighFive::File::OpenOrCreate);
        std::vector<size_t> dims;
        if (file.exist(dataSet)){
          auto dataset = file.getDataSet(dataSet);
          return dataset.getDimensions();
        } else {
          // Return empty vector if dataSet does not exist
        }
        return dims;
      }

  }; // class SafeFile

}; // namespace ChronusQ

