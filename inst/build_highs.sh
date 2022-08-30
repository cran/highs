#!/bin/bash

if test -z "${MAKE}"; then MAKE=`which make` 2> /dev/null; fi
if test -z "${MAKE}"; then MAKE=`which /Applications/Xcode.app/Contents/Developer/usr/bin/make` 2> /dev/null; fi


if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake4` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake3` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake2` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which cmake` 2> /dev/null; fi
if test -z "${CMAKE_EXE}"; then CMAKE_EXE=`which /Applications/CMake.app/Contents/bin/cmake` 2> /dev/null; fi


if test -z "${CMAKE_EXE}"; then
    echo "Could not find 'cmake'!"
    exit 1
fi


R_MACHINE=`R -s -e 'cat(unname(Sys.info()["machine"]))'`

echo ""
echo "arch: $(arch)"
echo "R_MACHINE: ${R_MACHINE}"
echo ""


: ${R_CC=`R CMD config CC`}
if test "$(uname -s)" = "Darwin"; then
    if test "${R_MACHINE}" = "arm64"; then
        echo "Detected 'darwin' with 'M1'"
        USE_DARWIN_M1="TRUE"
        cp -f inst/patches/CMakeLists_darwin-M1-static.txt inst/HiGHS/CMakeLists.txt
    else
        echo "Detected 'darwin' without 'M1'"
        cp -f inst/patches/CMakeLists_darwin-static.txt inst/HiGHS/CMakeLists.txt
    fi
else
    echo "Detected non 'Darwin'"
    # Detect the compiler with which R was built.
    if [[ "${R_CC}" == 'gcc'* ]]; then
        echo "  use GCC"
        R_CC="gcc"
        R_CXX="g++"
        cp -f inst/patches/CMakeLists_gcc-static.txt inst/HiGHS/CMakeLists.txt
    else
        echo "  use Clang"
        R_CC="clang"
        R_CXX="clang++"
        cp -f inst/patches/CMakeLists_clang-static.txt inst/HiGHS/CMakeLists.txt
    fi
fi


R_HIGHS_PKG_HOME=`pwd`
HIGHS_SRC_DIR=${R_HIGHS_PKG_HOME}/inst/HiGHS
R_HIGHS_BUILD_DIR=${HIGHS_SRC_DIR}/build
R_HIGHS_LIB_DIR=${R_HIGHS_PKG_HOME}/src/highslib

mkdir -p ${R_HIGHS_BUILD_DIR}
mkdir -p ${R_HIGHS_LIB_DIR}
cd ${R_HIGHS_BUILD_DIR}

if test -z "${USE_DARWIN_M1}"; then
    echo "CMAKE USE DEFAULT CMAKE"
    ${CMAKE_EXE} -DCMAKE_INSTALL_PREFIX=${R_HIGHS_LIB_DIR} ..
else
    echo "CMAKE USE DARWIN M1"
    ${CMAKE_EXE} -DCMAKE_INSTALL_PREFIX=${R_HIGHS_LIB_DIR} -DCMAKE_OSX_ARCHITECTURES="x86_64;arm64" ..
fi

${MAKE} install
