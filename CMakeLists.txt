# CMake script for Bio++ Core
# Author: Sylvain Gaillard and Julien Dutheil
# Created: 17/09/2010

# Global parameters
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(bpp-core CXX)
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_CXX_FLAGS "-Wall -Weffc++ -Wshadow -Wconversion")
IF(NOT NO_VIRTUAL_COV)
  SET(NO_VIRTUAL_COV FALSE CACHE BOOL
      "Disable covariant return type with virtual inheritance, for compilers that do not support it."
      FORCE)
ENDIF(NOT NO_VIRTUAL_COV)

IF(NO_VIRTUAL_COV)
  MESSAGE("-- Covariant return with virtual inheritance disabled.")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DNO_VIRTUAL_COV=1")
ENDIF(NO_VIRTUAL_COV)

# Libtool-like version number
# CURRENT:REVISION:AGE => file.so.(C-A).A.R
# current:  The most recent interface number that this library implements.
# revision: The implementation number of the current interface.
# age:      The difference between the newest and oldest interfaces that this
#           library implements.
# In other words, the library implements all the interface numbers in the
# range from number current - age to current.
SET(BPPCORE_VERSION_CURRENT "2")
SET(BPPCORE_VERSION_REVISION "3")
SET(BPPCORE_VERSION_AGE "0")

# Effective version number computation
MATH(EXPR BPPCORE_VERSION_MAJOR "${BPPCORE_VERSION_CURRENT} - ${BPPCORE_VERSION_AGE}")
SET(BPPCORE_VERSION_MINOR ${BPPCORE_VERSION_AGE})
SET(BPPCORE_VERSION_PATCH ${BPPCORE_VERSION_REVISION})
SET(BPPCORE_VERSION "${BPPCORE_VERSION_MAJOR}.${BPPCORE_VERSION_MINOR}.${BPPCORE_VERSION_PATCH}")

# Set the CMAKE_PREFIX_PATH for the find_library fonction when using non
# standard install location
IF(CMAKE_INSTALL_PREFIX)
  SET(CMAKE_PREFIX_PATH "${CMAKE_INSTALL_PREFIX}" ${CMAKE_PREFIX_PATH})
ENDIF(CMAKE_INSTALL_PREFIX)

# Subdirectories
ADD_SUBDIRECTORY(src)

# Doxygen
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
  ADD_CUSTOM_TARGET (apidoc cp Doxyfile ${CMAKE_BINARY_DIR}/Doxyfile-build
    COMMAND echo "OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}" >> ${CMAKE_BINARY_DIR}/Doxyfile-build
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile-build
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
  ADD_CUSTOM_TARGET (apidoc-stable cp Doxyfile ${CMAKE_BINARY_DIR}/Doxyfile-stable
    COMMAND echo "OUTPUT_DIRECTORY=${CMAKE_BINARY_DIR}" >> ${CMAKE_BINARY_DIR}/Doxyfile-stable
    COMMAND echo "HTML_HEADER=header.html" >> ${CMAKE_BINARY_DIR}/Doxyfile-stable
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile-stable
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
ENDIF (DOXYGEN_FOUND)

# Packager
SET(CPACK_PACKAGE_NAME "libbpp-core")
SET(CPACK_PACKAGE_VENDOR "Bio++ Development Team")
SET(CPACK_PACKAGE_VERSION "2.2.0")
SET(CPACK_PACKAGE_VERSION_MAJOR "2")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "The Bio++ Core library")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/COPYING.txt")
SET(CPACK_RESOURCE_FILE_AUTHORS "${CMAKE_SOURCE_DIR}/AUTHORS.txt")
SET(CPACK_RESOURCE_FILE_INSTALL "${CMAKE_SOURCE_DIR}/INSTALL.txt")
SET(CPACK_SOURCE_GENERATOR "TGZ")
SET(CPACK_SOURCE_IGNORE_FILES
 "CMakeFiles"
 "Makefile"
 "_CPack_Packages"
 "CMakeCache.txt"
 ".*\\\\.cmake"
 ".*\\\\.git"
 ".*\\\\.gz"
 ".*\\\\.deb"
 ".*\\\\.rpm"
 ".*\\\\.dmg"
 ".*\\\\..*\\\\.swp"
 "src/\\\\..*"
 "src/libbpp*"
 "debian/tmp"
 "debian/libbpp.*/"
 "debian/libbpp.*\\\\.so.*"
 "debian/libbpp.*\\\\.a"
 "debian/libbpp.*\\\\.substvars"
 "debian/libbpp.*\\\\.debhelper"
 "debian/debhelper\\\\.log"
 "html"
 "Core.tag"
 "Testing"
 "build-stamp"
 "install_manifest.txt"
 "DartConfiguration.tcl"
 ${CPACK_SOURCE_IGNORE_FILES}
)
IF (APPLE)
  SET(CPACK_BUNDLE_NAME "libbpp-core")
  SET(CPACK_BUNDLE_PLIST "bpp-core-Info.plist")
  SET(CPACK_BUNDLE_ICON "BppIcon.icns")
  SET(CPACK_GENERATOR "Bundle")
ENDIF()

SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
SET(CPACK_DEBSOURCE_PACKAGE_FILE_NAME "lib${CMAKE_PROJECT_NAME}_${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}.orig")
INCLUDE(CPack)

#This adds the 'dist' target
ADD_CUSTOM_TARGET(dist COMMAND ${CMAKE_MAKE_PROGRAM} package_source)
# 'clean' is not (yet) a first class target. However, we need to clean the directories before building the sources:
IF("${CMAKE_GENERATOR}" MATCHES "Make")
  ADD_CUSTOM_TARGET(make_clean
  COMMAND ${CMAKE_MAKE_PROGRAM} clean
  WORKING_DIRECTORY ${CMAKE_CURRENT_DIR}
  )
  ADD_DEPENDENCIES(dist make_clean)
ENDIF()

IF (UNIX)
#This creates deb packages:
ADD_CUSTOM_TARGET(origdist COMMAND cp ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz ../${CPACK_DEBSOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(origdist dist)
ADD_CUSTOM_TARGET(deb dpkg-buildpackage -uc -us -i${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(deb origdist)

#This creates rpm packages:
ADD_CUSTOM_TARGET(rpm rpmbuild -ta ${CPACK_SOURCE_PACKAGE_FILE_NAME}.tar.gz)
ADD_DEPENDENCIES(rpm dist)

ENDIF()

SET(CTEST_UPDATE_TYPE git)
SET(UPDATE_COMMAND "git")
SET(UPDATE_OPTIONS "")

ENABLE_TESTING()
INCLUDE(CTest)
IF (BUILD_TESTING)
  ADD_SUBDIRECTORY(test)
ENDIF(BUILD_TESTING)
