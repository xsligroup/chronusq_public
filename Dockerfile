FROM ubuntu:24.04 AS compile_stage

# Container system setting
RUN apt-get -y update && \
    apt-get -y install wget make git cmake

# Install CQ dependency

## install gcc, g++ and gfortran
RUN apt-get -y install gcc g++ gfortran 
## install eigen3 thru epel
RUN apt-get -y install libeigen3-dev
## install openblas
RUN apt-get -y install libopenblas-openmp-dev
## install libz
RUN apt-get -y install zlib1g-dev

ARG PACKAGE_INSTALL_PATH="/usr/local"

### install hdf5
WORKDIR /opt/
RUN wget -q "https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz" && \
    tar -xzf hdf5_1.14.6.tar.gz && \
    rm hdf5_1.14.6.tar.gz && \
    mkdir hdf5_1.14.6-build && \
    cd hdf5_1.14.6-build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLE=OFF -DHDF5_BUILD_CPP_LIB=ON -DHDF5_ENABLE_PARALLEL=OFF ../hdf5-hdf5_1.14.6 && \
    cmake --build . --target install -j 4 && \
    cd /opt && rm -rf hdf5-hdf5_1.14.6 hdf5_1.14.6-build

# Compile libint2 seperately
WORKDIR /opt/
RUN git clone https://github.com/xsligroup/libint-cq.git && \ 
    cd /opt/libint-cq/ && \
    git checkout 2.7.0-beta.6 && \
    mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH -DCMAKE_UNITY_BUILD=ON .. && \ 
    cmake --build . --target install -j 4 && \
    cd /opt/ && rm -rf libint-cq 

# Compile chronusq
RUN git clone https://github.com/xsligroup/chronusq_public.git
WORKDIR /opt/chronusq_public/
RUN mkdir build && cd build && \
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$PACKAGE_INSTALL_PATH -DCMAKE_PREFIX_PATH=$PACKAGE_INSTALL_PATH -DLibint2_ROOT=$PACKAGE_INSTALL_PATH .. && \
    cmake --build . --target install -j 4

# Copy chornusq to a lighter container
FROM ubuntu:24.04
COPY --from=compile_stage /opt/chronusq_public/basis /opt/chronusq_public/basis 
COPY --from=compile_stage /lib64/ /lib64/
COPY --from=compile_stage /lib/x86_64-linux-gnu/libgomp* /lib/x86_64-linux-gnu/
COPY --from=compile_stage /lib/x86_64-linux-gnu/libopenblas* /lib/x86_64-linux-gnu/
COPY --from=compile_stage /lib/x86_64-linux-gnu/libquadmath* /lib/x86_64-linux-gnu/
COPY --from=compile_stage /lib/x86_64-linux-gnu/libgfortran* /lib/x86_64-linux-gnu/
COPY --from=compile_stage /usr/local/lib/* /lib64/
COPY --from=compile_stage /usr/local/bin/chronusq /usr/bin/
ENV LD_LIBRARY_PATH=/lib64

WORKDIR /home/chronusq/
ENTRYPOINT ["chronusq"]

