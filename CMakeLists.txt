cmake_minimum_required(VERSION 3.14)

include(ExternalProject)
include(FetchContent)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(HOLYBUILD "Static compilation with holy-build-box" OFF)

set(CONDA_INCLUDE "")
if(DEFINED CONDAPREFIX)
  set(CONDA_INCLUDE "-I${CONDAPREFIX}/include")
  # include_target_directories("${CONDAPREFIX}/include")
  message(STATUS "Using conda: ${CONDA_INCLUDE}")
endif()

project(SVDSS VERSION 2.0.0)
add_executable(SVDSS assembler.cpp bam.cpp caller.cpp chromosomes.cpp clipper.cpp clusterer.cpp config.cpp ping_pong.cpp sfs.cpp smoother.cpp sv.cpp main.cpp)
set_target_properties(SVDSS PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
target_compile_options(SVDSS PRIVATE -Wall -Wextra)

##########################
### setup dependencies ###
##########################

# pthread, OpenMP and others libs
##################################
if (NOT HOLYBUILD)
  find_package(Threads REQUIRED)
  find_package(OpenMP REQUIRED)
else ()
  # ugly workaround
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  # libbz2 and liblzma are statically compiled for the portable binary
  set(BZIP_LIB_DIR ${CMAKE_SOURCE_DIR}/bzip2-1.0.6)
  find_library(LBZ2 libbz2.a HINTS ${BZIP_LIB_DIR} REQUIRED)
  set(LZMA_LIB_DIR ${CMAKE_SOURCE_DIR}/xz-5.2.2/src/liblzma/.libs/)
  find_library(LLZMA liblzma.a HINTS ${LZMA_LIB_DIR} REQUIRED)
endif ()

if(NOT DEFINED CONDAPREFIX)
  # LIBDEFLATE
  #############
  message(STATUS "libdeflate will be built from source")
  ExternalProject_Add(deflate
    GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git
    GIT_TAG 020133854ff73b8506fe59f92a9b5b622d360716 # v1.15
    BUILD_IN_SOURCE 1
    INSTALL_COMMAND "")
  ExternalProject_Get_Property(deflate SOURCE_DIR)
  SET(DEFLATE_SOURCE_DIR ${SOURCE_DIR})
  add_library(DEFLATE STATIC IMPORTED)
  set_target_properties(DEFLATE PROPERTIES IMPORTED_LOCATION ${DEFLATE_SOURCE_DIR}/libdeflate.a)
  add_dependencies(DEFLATE deflate)

  # HTSLIB
  #########
  SET(CONF_CMD autoheader && autoreconf -i && ./configure --disable-libcurl --disable-gcs --with-libdeflate "CC=${CMAKE_C_COMPILER}" "CFLAGS=-O3 -I${DEFLATE_SOURCE_DIR}" "LDFLAGS=-L${DEFLATE_SOURCE_DIR}")
  if (HOLYBUILD)
    # we have to include local libraries (libbz2 and liblzma)
    SET(CONF_CMD autoheader && autoreconf -i && ./configure --disable-libcurl --disable-gcs --with-libdeflate "CC=${CMAKE_C_COMPILER}" "CFLAGS=-O3 -I${DEFLATE_SOURCE_DIR} -I${BZIP_LIB_DIR}" "LDFLAGS=-L${DEFLATE_SOURCE_DIR} -L${BZIP_LIB_DIR} -L${LZMA_LIB_DIR}")
  endif ()
  message(STATUS "htslib will be built from source")
  ExternalProject_Add(htslib
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG 4e61c128238f3e7cbb3b1f4e9c0fdb4880aa9a10 # v1.16
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ${CONF_CMD}
    BUILD_COMMAND make
    "CC=${CMAKE_C_COMPILER}" libhts.a
    INSTALL_COMMAND "")
  ExternalProject_Get_Property(htslib SOURCE_DIR)
  SET(HTS_SOURCE_DIR ${SOURCE_DIR})
  SET(HTS_INCLUDE_DIR ${HTS_SOURCE_DIR})
  add_library(HTSLIB STATIC IMPORTED)
  set_target_properties(HTSLIB PROPERTIES IMPORTED_LOCATION ${HTS_SOURCE_DIR}/libhts.a)
  add_dependencies(HTSLIB htslib)
  add_dependencies(htslib DEFLATE)
else()
  message(STATUS "libdeflate from conda")
  message(STATUS "htslib from conda")
endif()

# abPOA
########
SET(PATCH_CMD "")
if(DEFINED CONDAPREFIX OR HOLYBUILD)
  # we have to remove the -native flag
  message(STATUS "abpoa will be build from source (no -native)")
  SET(PATCH_CMD sed -i "/native/d" CMakeLists.txt)
else ()
  message(STATUS "abpoa will be build from source (with -native)")
endif ()
ExternalProject_Add(abpoa
  GIT_REPOSITORY https://github.com/yangao07/abPOA.git
  GIT_TAG e6bb6fdfa40d573558e5e2b545bdf5769631eaaf # v1.5.3 with bug fixes
  PATCH_COMMAND ${PATCH_CMD}
  INSTALL_COMMAND ""
)
ExternalProject_Get_Property(abpoa SOURCE_DIR BINARY_DIR)
SET(ABPOA_SOURCE_DIR ${SOURCE_DIR})
SET(ABPOA_INCLUDE_DIR ${ABPOA_SOURCE_DIR}/include)
# NOTE: abpoa library will be set manually
# add_library(ABPOA STATIC IMPORTED)
# set_target_properties(ABPOA PROPERTIES IMPORTED_LOCATION ${ABPOA_SOURCE_DIR}/lib/libabpoa.a)

# ksw2
SET(KSW2_MAKE make)
if(DEFINED CONDAPREFIX)
  SET(KSW2_MAKE make "CC=${CMAKE_C_COMPILER}" "CFLAGS=-I${CONDAPREFIX}/include" "LDFLAGS=-L${CONDAPREFIX}/lib")
endif ()
message(STATUS "ksw2 will be built from source")
ExternalProject_Add(ksw2
    GIT_REPOSITORY https://github.com/lh3/ksw2.git
    GIT_TAG ""
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ${KSW2_MAKE}
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(ksw2 SOURCE_DIR)
SET(KSW2_INCLUDE_DIR ${SOURCE_DIR})

# rapidfuzz-cpp
################
message(STATUS "rapidfuzz-cpp will be built from source")
FetchContent_Declare(rapidfuzz
  GIT_REPOSITORY https://github.com/maxbachmann/rapidfuzz-cpp.git
  GIT_TAG 5412d5d877518e7754394bdbec76e45c8187c631) # v1.10.4
FetchContent_MakeAvailable(rapidfuzz)

# interval-tree
################
message(STATUS "interval-tree will be cloned")
ExternalProject_Add(intervaltree
    GIT_REPOSITORY https://github.com/5cript/interval-tree.git
    GIT_TAG 309b9c725191d4bb1d134f28a8a32ad2f68a8ffa
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(intervaltree SOURCE_DIR)
SET(INTERVALTREE_SOURCE_DIR ${SOURCE_DIR})
SET(INTERVALTREE_INCLUDE_DIR ${INTERVALTREE_SOURCE_DIR}/include/interval-tree)

# ropebwt3
##########
message(STATUS "ropebwt3 will be built from source")
ExternalProject_Add(ropebwt
    GIT_REPOSITORY https://github.com/lh3/ropebwt3.git
    GIT_TAG 0ea3919ed21f10d857d508b0ea728a27abbd8a35
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    PATCH_COMMAND "" # bash ${CMAKE_SOURCE_DIR}/patches/patch_ropebwt.sh <SOURCE_DIR>/mrope.c # I haven't been able to run sed here.. I could use my branch btw
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make
      "CC=${CMAKE_C_COMPILER}"
      "CFLAGS=-g -Wall -O2 -fopenmp ${CONDA_INCLUDE}"
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(ropebwt SOURCE_DIR)
SET(ROPEBWT_SOURCE_DIR ${SOURCE_DIR})
SET(ROPEBWT_INCLUDE_DIR ${ROPEBWT_SOURCE_DIR})
set(ropebwt-src ${ROPEBWT_SOURCE_DIR}/build.c ${ROPEBWT_SOURCE_DIR}/sais-ss.c ${ROPEBWT_SOURCE_DIR}/libsais16x64.c ${ROPEBWT_SOURCE_DIR}/fm-index.c ${ROPEBWT_SOURCE_DIR}/rld0.c ${ROPEBWT_SOURCE_DIR}/mrope.c ${ROPEBWT_SOURCE_DIR}/rope.c ${ROPEBWT_SOURCE_DIR}/io.c ${ROPEBWT_SOURCE_DIR}/rle.c ${ROPEBWT_SOURCE_DIR}/kthread.c ${ROPEBWT_SOURCE_DIR}/kalloc.c ${ROPEBWT_SOURCE_DIR}/misc.c ${ROPEBWT_SOURCE_DIR}/ssa.c ${ROPEBWT_SOURCE_DIR}/dawg.c ${ROPEBWT_SOURCE_DIR}/libsais16.c)
set_source_files_properties(${ropebwt-src} PROPERTIES GENERATED TRUE)
add_library(ROPEBWT OBJECT ${ropebwt-src})
add_dependencies(ROPEBWT ropebwt)

# spdlog
########
message(STATUS "spdlog will be built from source")
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG 7e635fca68d014934b4af8a1cf874f63989352b7) # v1.12.0
FetchContent_MakeAvailable(spdlog)

############################################
### final setup (includes and libraries) ###
############################################

add_dependencies(SVDSS abpoa rapidfuzz intervaltree ropebwt spdlog ksw2)
if(NOT DEFINED CONDAPREFIX)
  # we need to build them
  add_dependencies(SVDSS deflate htslib)
endif()

target_include_directories(SVDSS
  PRIVATE ${HTS_INCLUDE_DIR}
  PRIVATE ${ABPOA_INCLUDE_DIR}
  PRIVATE ${KSW2_INCLUDE_DIR}
  PRIVATE ${RAPIDFUZZ_INCLUDE_DIR}
  PRIVATE ${INTERVALTREE_INCLUDE_DIR}
  PRIVATE ${ROPEBWT_INCLUDE_DIR}
  PRIVATE ${SPDLOG_INCLUDE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

target_link_options(SVDSS PRIVATE -static-libgcc -static-libstdc++)

if(DEFINED CONDAPREFIX)
  # link to conda
  find_library(HTSLIB libhts.so REQUIRED)
  find_library(DEFLATE libdeflate.so REQUIRED)
  target_link_libraries(SVDSS
    PUBLIC hts
    PUBLIC deflate
  )
else()
  # link to local
  target_link_libraries(SVDSS
    PUBLIC HTSLIB
    PUBLIC DEFLATE
  )
endif()

target_link_libraries(SVDSS
  PUBLIC rapidfuzz::rapidfuzz
  PUBLIC ${BINARY_DIR}/lib/libabpoa.a
  PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extz2_sse.o
  PUBLIC ${KSW2_INCLUDE_DIR}/ksw2_extd2_sse.o
  PUBLIC ROPEBWT
  PUBLIC spdlog::spdlog
  PUBLIC z
)

if (NOT HOLYBUILD)
  # link to system
  target_link_libraries(SVDSS
    PUBLIC lzma
    PUBLIC bz2
    PUBLIC pthread
    PUBLIC OpenMP::OpenMP_CXX
  )
else ()
  # link to local. no openmp since we manually added -fopenmp
  target_link_libraries(SVDSS
    PUBLIC ${LLZMA}
    PUBLIC ${LBZ2}
  )
endif ()
