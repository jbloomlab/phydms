# CMake script for bpp-seq unit tests
# Author: Julien Dutheil
# Created: 30/10/2010

MACRO(TEST_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
  #start:
  FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})

  SET(${lib_name}_NAMES ${lib_name} ${lib_name}.lib ${lib_name}.dll)
  FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})
  IF(${lib_name}_LIBRARY)
    MESSAGE("-- Library ${lib_name} found here:")
    MESSAGE("   includes: ${${lib_name}_INCLUDE_DIR}")
    MESSAGE("   dynamic libraries: ${${lib_name}_LIBRARY}")
    MESSAGE(WARNING "Library ${lib_name} is already installed in the system tree. Test will be built against it. This may lead to unexpected results. You may want to do 'make install' before 'make test', or remove the installed version.")
  ELSE()
    SET(${lib_name}_LIBRARY "-L../src -lbpp-seq")
    SET(${lib_name}_INCLUDE_DIR "../src/")
  ENDIF()
  INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
  SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(TEST_FIND_LIBRARY)

#Find the bpp-seq library library:
TEST_FIND_LIBRARY(LIBS bpp-seq Bpp/Seq/Alphabet/Alphabet.h)

ADD_EXECUTABLE(test_alphabets test_alphabets.cpp)
TARGET_LINK_LIBRARIES(test_alphabets ${LIBS})
ADD_TEST(test_alphabets "test_alphabets")

ADD_EXECUTABLE(test_sequences test_sequences.cpp)
TARGET_LINK_LIBRARIES(test_sequences ${LIBS})
ADD_TEST(test_sequences "test_sequences")

ADD_EXECUTABLE(test_io test_io.cpp)
TARGET_LINK_LIBRARIES(test_io ${LIBS})
ADD_TEST(test_io "test_io")

ADD_EXECUTABLE(test_containers test_containers.cpp)
TARGET_LINK_LIBRARIES(test_containers ${LIBS})
ADD_TEST(test_containers "test_containers")

ADD_EXECUTABLE(test_alignment_scores test_alignment_scores.cpp)
TARGET_LINK_LIBRARIES(test_alignment_scores ${LIBS})
ADD_TEST(test_alignment_scores "test_alignment_scores")

ADD_EXECUTABLE(test_walker test_walker.cpp)
TARGET_LINK_LIBRARIES(test_walker ${LIBS})
ADD_TEST(test_walker "test_walker")

ADD_EXECUTABLE(test_bowker test_bowker.cpp)
TARGET_LINK_LIBRARIES(test_bowker ${LIBS})
ADD_TEST(test_bowker "test_bowker")

IF(UNIX)
  SET_PROPERTY(TEST test_alphabets test_sequences test_io test_containers test_alignment_scores test_walker test_bowker PROPERTY ENVIRONMENT "LD_LIBRARY_PATH=$ENV{LD_LIBRARY_PATH}:../src")
ENDIF()

IF(APPLE)
  SET_PROPERTY(TEST test_alphabets test_sequences test_io test_containers test_alignment_scores test_walker test_bowker PROPERTY ENVIRONMENT "DYLD_LIBRARY_PATH=$ENV{DYLD_LIBRARY_PATH}:../src")
ENDIF()

IF(WIN32)
  SET(ENV{PATH} "$ENV{PATH};..\\src")
ENDIF()

