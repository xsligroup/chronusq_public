<div align="center">
  <img src="cq_logo.png" height="220px"/>
</div>

Chronus Quantum 
===============

The Chronus Quantum (ChronusQ) Software Package [v. Beta] is a high-performance
computational chemistry software package with a strong emphasis on explicitly
time-dependent and post-SCF quantum mechanical methods.

* [Changelog](CHANGELOG.md)
* [Documentation](../../wikis/home)
* [Installation](#installation)



Installation
------------

### Docker Quickstart

ChronusQ is now available as a Docker image so you can quickly and easily test
it on your system! To do so, install
[Docker desktop](https://www.docker.com/get-started) and pull the latest
ChronusQ image:
```
docker pull uwligroup/chronusq
```

You can then run it on an input file by invoking `docker run`:
```
docker run -v ${PWD}:/home/chronusq uwligroup/chronusq <inputfile> 
```

For more details about using the Docker image, see the
[getting ChronusQ](../../wikis/getting-chronusq#docker-image) and
[running ChronusQ](../../wikis/running-chronusq#docker-image) wiki pages.

### Prerequisites for compilation

- C++14 compiler 
- C compiler (for LibXC)
- Fortran compiler (for LibXC)
- [CMake](http://cmake.org) build system (Version 3.11+).
- [HDF5](https://support.hdfgroup.org/HDF5/)
- [Eigen3](http://eigen.tuxfamily.org)

### Quick compilation

If you have all the prerequisites above, all you need to do is:

```
git clone https://urania.chem.washington.edu/chronusq/chronusq_public.git
mkdir chronusq_public/build && cd chronusq_public/build
cmake ..
cmake --build .
```
To install, you then run
```
cmake --build . --target install
```

For more details on installation requirements, and running ChronusQ, see the
[getting ChronusQ](../../wikis/getting-chronusq#compilation-from-source)
and [running ChronusQ](../../wikis/Running-ChronusQ#compiled-from-source)
wiki pages.


Citing ChronusQ
---------------
The following WIREs paper and software citation should be cited in publications using the ChronusQ package located in [CITE.txt](CITE.txt).


Found a bug or want a new feature?
----------------------------------
Please submit a bug report or feature request on the [issues](https://urania.chem.washington.edu/chronusq/chronusq_public/-/issues) page.


General Inquiries
-----------------
- Contact xsli at uw dot edu

