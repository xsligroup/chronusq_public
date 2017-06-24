#!/bin/bash

# A short script that generates the lcov html files
# It assumes that you have compiled with the CMake flag
# ENABLE_COVERAGE and run `make test`.

if ["${1}" == ""]; then
  build_dir="build"
else
  build_dir="${1}"
fi

echo "Looking in ${build_dir}"

lcov -c -d $build_dir -o "coverage.info"
mkdir -p coverage
genhtml -o coverage --demangle-cpp \
  --num-spaces 2 --legend \
  coverage.info
