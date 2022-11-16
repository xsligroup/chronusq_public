FROM centos:7 AS compile_stage

# Container system setting
RUN yum -y update && \
    yum -y install centos-release-scl scl-utils epel-release wget make git patch && \
    yum-config-manager --enable epel 

# Install CQ dependency

## install gcc, g++ and gfortran
RUN yum -y install devtoolset-9-gcc devtoolset-9-gcc-c++ devtoolset-9-gcc-gfortran 
SHELL ["/bin/scl", "enable", "devtoolset-9"]
## install cmake, 
RUN wget -qO- "https://cmake.org/files/v3.18/cmake-3.18.2-Linux-x86_64.tar.gz" | tar --strip-components=1 -xz -C /usr/local
## install eigen3 thru epel
RUN yum -y install eigen3 

ARG PACKAGE_INSTALL_PATH="/usr/local"

### install hdf5
WORKDIR /opt/
RUN wget -q "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.12/hdf5-1.12.0/src/CMake-hdf5-1.12.0.tar.gz" && \
    tar -xzf CMake-hdf5-1.12.0.tar.gz && \
    rm CMake-hdf5-1.12.0.tar.gz && \
    cd CMake-hdf5-1.12.0 && \
    tar -xzf SZip.tar.gz && \ 
    tar -xzf ZLib.tar.gz && \
    cd /opt/CMake-hdf5-1.12.0/SZip && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH .. && cmake --build . --target install && \
    cd /opt/CMake-hdf5-1.12.0/ZLib && mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH .. && cmake --build . --target install && \
    cd /opt/CMake-hdf5-1.12.0 && ctest -S HDF5config.cmake,BUILD_GENERATOR=Unix,INSTALLDIR=$PACKAGE_INSTALL_PATH -C Release -V -O hdf5.log && \
    cd build && cmake --build . --target install && \
    cd /opt && rm -rf CMake-hdf5-1.12.0

# Compile libint2 seperately
WORKDIR /opt/
RUN git clone https://urania.chem.washington.edu/chronusq/libint-cq.git && \ 
    cd /opt/libint-cq/ && \
    git checkout 2.7.0-beta.6 && \
    mkdir build && cd build && cmake -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH -DCMAKE_UNITY_BUILD=ON .. && \ 
    cmake --build . --target install -j 5 && \
    cd /opt/ && rm -rf libint-cq 

# Compile openblas
WORKDIR /opt/
RUN git clone https://github.com/xianyi/OpenBLAS.git && \
    cd /opt/OpenBLAS/ && \
    git checkout v0.3.21 && \
    make -j 5 CFLAGS="-Wno-error=implicit-function-declaration ${CFLAGS}" && \
    make install $PACKAGE_INSTALL_PATH

# Compile chronusq
COPY . /opt/chronusq_public/
WORKDIR /opt/chronusq_public/
RUN mkdir build && cd build && \
    cmake -DOPENBLAS_DYNAMIC_ARCH=ON \
          -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH \
          -DCMAKE_CXX_FLAGS='-O3' \
          -DCMAKE_C_FLAGS='-O3' \
          -DCMAKE_Fortran_FLAGS='-O3' \
          .. && \
    cmake --build . --target install -j 5 && \
    cd /opt/chronusq_public/ && rm -rf build 

# Copy chornusq to a lighter container
FROM alpine:3.12
COPY --from=compile_stage /opt/chronusq_public/basis /opt/chronusq_public/basis 
COPY --from=compile_stage /lib64/ /lib64/
COPY --from=compile_stage /usr/local/lib/libhdf5* /lib64/
COPY --from=compile_stage /usr/local/lib/libopenblas* /lib64/
COPY --from=compile_stage /usr/local/bin/chronusq /usr/bin/
COPY --from=compile_stage /usr/local/lib64/* /lib64/


WORKDIR /home/chronusq/
ENTRYPOINT ["chronusq"]

