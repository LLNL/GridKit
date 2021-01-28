#
# Copyright (c) 2017, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
# Written by Slaven Peles <peles2@llnl.gov>.
# LLNL-CODE-718378.
# All rights reserved.
#
# This file is part of GridKit™. For details, see github.com/LLNL/GridKit
# Please also read the LICENSE file.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# - Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
# - Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
# - Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL
# SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
# OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISINGIN ANY
# WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.
#
# Lawrence Livermore National Laboratory is operated by Lawrence Livermore
# National Security, LLC, for the U.S. Department of Energy, National
# Nuclear Security Administration under Contract DE-AC52-07NA27344.
#
# This document was prepared as an account of work sponsored by an agency
# of the United States government. Neither the United States government nor
# Lawrence Livermore National Security, LLC, nor any of their employees
# makes any warranty, expressed or implied, or assumes any legal liability
# or responsibility for the accuracy, completeness, or usefulness of any
# information, apparatus, product, or process disclosed, or represents that
# its use would not infringe privately owned rights. Reference herein to
# any specific commercial product, process, or service by trade name,
# trademark, manufacturer, or otherwise does not necessarily constitute or
# imply its endorsement, recommendation, or favoring by the United States
# government or Lawrence Livermore National Security, LLC. The views and
# opinions of authors expressed herein do not necessarily state or reflect
# those of the United States government or Lawrence Livermore National
# Security, LLC, and shall not be used for advertising or product
# endorsement purposes.
#

cmake_minimum_required(VERSION 3.12)

project(gridkit)

set(component_models
  #generic_model
)

set(solver_libs
  # generic_solver
)

# Ipopt support is disabled by default
option(GRIDKIT_ENABLE_IPOPT "Enable Ipopt support" ON)

# SUNDIALS support is enabled by default
option(GRIDKIT_ENABLE_SUNDIALS "Enable SUNDIALS support" ON)

# Enable KLU
option(GRIDKIT_ENABLE_SUNDIALS_SPARSE "Enable SUNDIALS sparse linear solvers" OFF)


set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/config)

# TODO: Probably beter to set a debug interface target
set(CMAKE_CXX_FLAGS_DEBUG "-Wall -O0 -g")

set(CMAKE_CXX_STANDARD 11)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

option(BUILD_SHARED_LIBS:BOOL "Build shared libraries" ON)

if(GRIDKIT_ENABLE_IPOPT)
  include(FindIpopt)
  set(CMAKE_INSTALL_RPATH  ${IPOPT_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
endif()

include(FindSUNDIALS)
set(CMAKE_INSTALL_RPATH  ${SUNDIALS_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})

if(GRIDKIT_ENABLE_SUNDIALS_SPARSE)
  include(FindSuiteSparse)
  set(CMAKE_INSTALL_RPATH  ${SUITESPARSE_LIBRARY_DIR} ${CMAKE_INSTALL_RPATH})
endif()

include_directories(${SUITESPARSE_INCLUDE_DIR} ${SUNDIALS_INCLUDE_DIR} ${IPOPT_INCLUDE_DIR})

# Set up configuration header file
# configure_file(gridkit_config.h.in gridkit_config.h)

add_subdirectory(ComponentLib)
add_subdirectory(Solver)
add_subdirectory(Examples)

# install(FILES gridkit_config.hpp DESTINATION include)

# TESTING
enable_testing()
add_test(NAME AdjointSens COMMAND $<TARGET_FILE:adjoint> )
add_test(NAME Grid3Bus    COMMAND $<TARGET_FILE:grid3bus>)
if(GRIDKIT_ENABLE_IPOPT)
  add_test(NAME DynamicConOpt  COMMAND $<TARGET_FILE:dynconopt>   )
  add_test(NAME GenConstLoad   COMMAND $<TARGET_FILE:genconstload>)
  add_test(NAME GenInfiniteBus COMMAND $<TARGET_FILE:geninfbus>   )
  add_test(NAME ParameterEst   COMMAND $<TARGET_FILE:paramest> ${CMAKE_SOURCE_DIR}/Examples/ParameterEstimation/lookup_table.dat)
endif()
