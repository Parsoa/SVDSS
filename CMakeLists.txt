set(CMAKE_VERBOSE_MAKEFILE on)

cmake_minimum_required(VERSION 3.14)

include(ExternalProject)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
#set(BUILD_SHARED_LIBS OFF)
#set(MALLOC_STATIC_BUILD_FLAG "--enable-static")
#link_libraries("-static")

set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG -Wall -Wno-comment -Wno-unused-function -Wno-narrowing -Wno-sign-compare -Wno-unused-variable ${CMAKE_CXX_FLAGS}")

project(TEST VERSION 1.0)

add_executable(SVDSS assembler.cpp bam.cpp bed_utils.cpp caller.cpp chromosomes.cpp clipper.cpp cluster.cpp config.cpp extender.cpp lprint.cpp ping_pong.cpp sfs.cpp smoother.cpp sv.cpp vcf.cpp main.cpp)
set_target_properties(SVDSS PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")

# pthread and OpenMP
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

# LIBDEFLATE
#############
#find_library(DEFLATE libdeflate.a)
#if (NOT DEFLATE)
message(STATUS "libdeflate will be built form source")
ExternalProject_Add(deflate
  GIT_REPOSITORY https://github.com/ebiggers/libdeflate.git
  GIT_TAG 047aa84e01b38d82f3612810e357bd40f14a3d39
  BUILD_IN_SOURCE 1
  UPDATE_COMMAND ""
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make PREFIX="."
  "CC=${CMAKE_C_COMPILER}"
  "CFLAGS=-fPIC -O3" libdeflate.a
  INSTALL_COMMAND "")
ExternalProject_Get_Property(deflate SOURCE_DIR)
SET(DEFLATE_SOURCE_DIR ${SOURCE_DIR})
add_library(DEFLATE STATIC IMPORTED)
set_target_properties(DEFLATE PROPERTIES IMPORTED_LOCATION ${DEFLATE_SOURCE_DIR}/libdeflate.a)
add_dependencies(DEFLATE deflate)
#else ()
#  message(STATUS "libdeflate found ${DEFLATE}")
#endif ()

# HTSLIB
#########
#find_library(HTSLIB libhts.a)
#if (NOT HTSLIB)
message(STATUS "htslib will be built form source")
ExternalProject_Add(htslib
  GIT_REPOSITORY https://github.com/samtools/htslib.git
  GIT_TAG 2cd99cf333938a9eb70393bc59379e1d6ea66b37
  UPDATE_COMMAND ""
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND autoheader && autoreconf && ./configure --disable-libcurl --disable-gcs --with-libdeflate
  "CC=${CMAKE_C_COMPILER}"
  "CFLAGS=-O3 -I${DEFLATE_SOURCE_DIR}"
  "LDFLAGS=-L${DEFLATE_SOURCE_DIR}"
  BUILD_COMMAND make libhts.a
  INSTALL_COMMAND "")
ExternalProject_Get_Property(htslib SOURCE_DIR)
SET(HTS_SOURCE_DIR ${SOURCE_DIR})
SET(HTS_INCLUDE_DIR ${HTS_SOURCE_DIR})
add_library(HTSLIB STATIC IMPORTED)
set_target_properties(HTSLIB PROPERTIES IMPORTED_LOCATION ${HTS_SOURCE_DIR}/libhts.a)
add_dependencies(HTSLIB htslib)
add_dependencies(htslib DEFLATE)
#else ()
#    message(STATUS "htslib found ${HTSLIB}")
#endif ()

# abPOA
########
#find_library(ABPOA libabpoa.a)
#if (NOT ABPOA)
message(STATUS "abPOA will be built form source")
ExternalProject_Add(abpoa
  GIT_REPOSITORY https://github.com/yangao07/abPOA.git
  GIT_TAG 42a09de0d2ef076c4b173ff4bd8a9812ded26707
  UPDATE_COMMAND ""
  BUILD_IN_SOURCE 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make
  INSTALL_COMMAND "")
ExternalProject_Get_Property(abpoa SOURCE_DIR)
SET(ABPOA_SOURCE_DIR ${SOURCE_DIR})
SET(ABPOA_INCLUDE_DIR ${ABPOA_SOURCE_DIR}/include)
add_library(ABPOA STATIC IMPORTED)
set_target_properties(ABPOA PROPERTIES IMPORTED_LOCATION ${ABPOA_SOURCE_DIR}/lib/libabpoa.a)
#else ()
#    message(STATUS "abPOA found ${ABPOA}")
#endif ()

# parasail
###########
#find_library(PARASAIL libparasail.a)
#if (NOT PARASAIL)
message(STATUS "parasail will be built form source")
ExternalProject_Add(parasail
  GIT_REPOSITORY https://github.com/jeffdaily/parasail.git
  GIT_TAG e236a6feec7bc84fb81923d783acb345e538a281
  UPDATE_COMMAND ""
  BUILD_IN_SOURCE 1
  CMAKE_ARGS -DBUILD_SHARED_LIBS=OFF
  INSTALL_COMMAND ""
  )
ExternalProject_Get_Property(parasail SOURCE_DIR)
SET(PARASAIL_SOURCE_DIR ${SOURCE_DIR})
SET(PARASAIL_INCLUDE_DIR ${PARASAIL_SOURCE_DIR})
add_library(PARASAIL STATIC IMPORTED)
set_target_properties(PARASAIL PROPERTIES IMPORTED_LOCATION ${PARASAIL_SOURCE_DIR}/libparasail.a)
#else ()
#    message(STATUS "parasail found ${PARASAIL}")
#endif ()

# rapidfuzz-cpp
################
ExternalProject_Add(rapidfuzz
    GIT_REPOSITORY https://github.com/maxbachmann/rapidfuzz-cpp.git
    GIT_TAG d1e82379395cafc6d439c1c1e2cbe7512eaf2518
    BUILD_IN_SOURCE 1
    UPDATE_COMMAND ""
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(rapidfuzz SOURCE_DIR)
SET(RAPIDFUZZ_SOURCE_DIR ${SOURCE_DIR})
SET(RAPIDFUZZ_INCLUDE_DIR ${RAPIDFUZZ_SOURCE_DIR})

# interval-tree
################
ExternalProject_Add(intervaltree
    GIT_REPOSITORY https://github.com/5cript/interval-tree.git
    GIT_TAG 4d7c66ce4fb542a29a7d74039f8e7c92fae3d6f8
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(intervaltree SOURCE_DIR)
SET(INTERVALTREE_SOURCE_DIR ${SOURCE_DIR})
SET(INTERVALTREE_INCLUDE_DIR ${INTERVALTREE_SOURCE_DIR})

# ropebwt2
##########
ExternalProject_Add(ropebwt
    GIT_REPOSITORY https://github.com/lh3/ropebwt2.git
    GIT_TAG bd8dbd3db2e9e3cff74acc2907c0742c9ebbf033
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE 1
    CONFIGURE_COMMAND ""
    BUILD_COMMAND "make"
    INSTALL_COMMAND ""
)
ExternalProject_Get_Property(ropebwt SOURCE_DIR)
SET(ROPEBWT_SOURCE_DIR ${SOURCE_DIR})
SET(ROPEBWT_INCLUDE_DIR ${ROPEBWT_SOURCE_DIR})
set(ropebwt ${ROPEBWT_SOURCE_DIR}/mrope.c ${ROPEBWT_SOURCE_DIR}/rope.c ${ROPEBWT_SOURCE_DIR}/rld0.c ${ROPEBWT_SOURCE_DIR}/rle.c)
set_source_files_properties(ropebwt PROPERTIES GENERATED TRUE)
add_library(ROPEBWT OBJECT ${ropebwt})

add_dependencies(SVDSS deflate htslib abpoa parasail rapidfuzz intervaltree ropebwt)

target_include_directories(SVDSS
  PRIVATE ${HTS_INCLUDE_DIR}
  PRIVATE ${ABPOA_INCLUDE_DIR}
  PRIVATE ${PARASAIL_INCLUDE_DIR}
  PRIVATE ${RAPIDFUZZ_INCLUDE_DIR}
  PRIVATE ${INTERVALTREE_INCLUDE_DIR}
  PRIVATE ${ROPEBWT_INCLUDE_DIR}
  PRIVATE ${PROJECT_SOURCE_DIR}
)

#find_library(LZ libz.a REQUIRED)
#if (NOT LZ)
#  message(FATAL_ERROR "libz required")
#else ()
#  message(STATUS "libz found ${LZ}")
#  #add_library(Z STATIC IMPORTED)
#  #set_target_properties(Z PROPERTIES IMPORTED_LOCATION ${LZ})
#endif ()

#find_library(LLZMA liblzma.a REQUIRED)
#if (NOT LLZMA)
#  message(FATAL_ERROR "liblzma required")
#else ()
#  message(STATUS "liblzma found ${LLZMA}")
#endif ()

#find_library(LBZ2 libbz2.a REQUIRED)
#if (NOT LBZ2)
#  message(FATAL_ERROR "libbz2 required")
#else ()
#  message(STATUS "libbz2 found ${LBZ2}")
#endif ()

#find_library(LPT libpthread.a REQUIRED)
#if (NOT LPT)
#  message(FATAL_ERROR "pthread required")
#else ()
#  message(STATUS "pthread found ${LPT}")
#endif ()

#find_library(LOMP libgomp.a)
#if (NOT LOMP)
#  message(FATAL_ERROR "OpenMP required")
#else ()
#  message(STATUS "OpenMP found ${LOMP}")
#endif ()

#target_link_options(SVDSS PRIVATE -static-libgcc -static-libstdc++)

target_link_libraries(SVDSS
  PUBLIC HTSLIB
  PUBLIC DEFLATE
  PUBLIC ABPOA
  PUBLIC PARASAIL
  PUBLIC ROPEBWT
  PUBLIC pthread
  PUBLIC z
  PUBLIC lzma
  PUBLIC bz2
  PUBLIC OpenMP::OpenMP_CXX
)
