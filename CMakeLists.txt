# SPDX-License-Identifier: Apache-2.0

cmake_minimum_required(VERSION 3.13)
project(SQIsign VERSION 2.0 LANGUAGES C ASM)

set(SQISIGN_SO_VERSION "0")
set(CMAKE_C_STANDARD 11)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CTest)

option(ENABLE_STRICT "Build with strict compile options." ON)
option(ENABLE_TESTS  "Enable compilation of tests." ON)
option(ENABLE_CT_TESTING  "Enable compilation for constant time testing." OFF)
option(ENABLE_SIGN "Build with sign functionality" ON)
option(SQISIGN_USE_MLLL_GRAM "Route quat_lideal_prime_norm_reduced_equivalent through the Gram-based MLLL (Compact paper Alg 4) instead of HNF." OFF)
set(MLLL_FP_KIND "2" CACHE STRING "MLLL μ/B float backend: 0=double, 1=mpf_t, 2=mpfr_t (paper-strict default), 3=dpe_t")
set(MLLL_FP_PREC "1024" CACHE STRING "MLLL μ/B float precision in bits (covers L5 paper bound ~1010b)")
set(GMP_LIBRARY "SYSTEM" CACHE STRING "Which version of GMP to use: SYSTEM, BUILD or MINI")
set(GF_RADIX "AUTO" CACHE STRING "Set the radix for the gf module (currently supported values: 32 or 64), or AUTO.")

if (NOT DEFINED SQISIGN_BUILD_TYPE)
  SET(SQISIGN_BUILD_TYPE "ref")
endif()

if (${SQISIGN_BUILD_TYPE} MATCHES "ref")
  add_compile_definitions(SQISIGN_BUILD_TYPE_REF SQISIGN_GF_IMPL_REF)
elseif (${SQISIGN_BUILD_TYPE} MATCHES "opt")
  add_compile_definitions(SQISIGN_BUILD_TYPE_OPT SQISIGN_GF_IMPL_REF)
elseif (${SQISIGN_BUILD_TYPE} MATCHES "broadwell")
  add_compile_definitions(SQISIGN_BUILD_TYPE_BROADWELL SQISIGN_GF_IMPL_BROADWELL)
elseif (${SQISIGN_BUILD_TYPE} MATCHES "arm64crypto")
  add_compile_definitions(SQISIGN_BUILD_TYPE_ARM64CRYPTO SQISIGN_GF_IMPL_REF)
endif()

SET(SVARIANT_S "lvl1;lvl3;lvl5")

include(.cmake/flags.cmake)
include(.cmake/sanitizers.cmake)
include(.cmake/target.cmake)
if(ENABLE_SIGN)
	include(.cmake/gmpconfig.cmake)
	add_compile_definitions(ENABLE_SIGN)
endif()

if(SQISIGN_USE_MLLL_GRAM)
	add_compile_definitions(SQISIGN_USE_MLLL_GRAM)
endif()

add_compile_definitions(MLLL_FP_KIND=${MLLL_FP_KIND} MLLL_FP_PREC=${MLLL_FP_PREC})
if(MLLL_FP_KIND EQUAL 2)
	find_library(MPFR_LIB NAMES mpfr REQUIRED)
	message(STATUS "MLLL: using mpfr_t backend, linking ${MPFR_LIB}")
	link_libraries(${MPFR_LIB})
endif()
message(STATUS "MLLL: MLLL_FP_KIND=${MLLL_FP_KIND} MLLL_FP_PREC=${MLLL_FP_PREC}")

set(BM_BINS "" CACHE INTERNAL "List of benchmark binaries")

set(SELECT_IMPL_TYPE ${PROJECT_SOURCE_DIR}/.cmake/impl_type.cmake)
set(SELECT_SQISIGN_VARIANT ${PROJECT_SOURCE_DIR}/.cmake/sqisign_variant.cmake)

set(INC_PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
add_subdirectory(apps)
add_subdirectory(test)

include(.cmake/bm.cmake)

#if(ENABLE_TESTS)
#    enable_testing()
#    add_subdirectory(test)
#endif()
