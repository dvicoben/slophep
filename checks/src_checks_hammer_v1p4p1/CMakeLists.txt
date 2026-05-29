# set minimum required cmake version
cmake_minimum_required (VERSION 3.5 FATAL_ERROR)

project (UserProject LANGUAGES CXX)

set(Hammer_CompileOptions ;-D_FILE_OFFSET_BITS=64;-Wno-unknown-pragmas)
set(Hammer_PIC TRUE)
set(Hammer_Include $ENV{CONDA_PREFIX}/include)
set(Hammer_Libs )
SET(CMAKE_CXX_STANDARD 20)
set(Hammer_LIBTYPE "STATIC")

# don't use configuration-specific subdirectories in Xcode
if (CMAKE_GENERATOR STREQUAL "Xcode")
   set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PROJECT_SOURCE_DIR}/lib)
   set (CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PROJECT_SOURCE_DIR}/lib)
endif (CMAKE_GENERATOR STREQUAL "Xcode")

foreach(libname HammerCore Amplitudes FormFactors HammerTools HammerBase Rates)
    message(STATUS "Adding Dependency: ${libname}")
    if(Hammer_LIBTYPE STREQUAL "STATIC")
        add_library(${libname} STATIC IMPORTED)
        set_target_properties(${libname} PROPERTIES
            IMPORTED_LOCATION "$ENV{CONDA_PREFIX}/lib/Hammer/${CMAKE_STATIC_LIBRARY_PREFIX}${libname}${CMAKE_STATIC_LIBRARY_SUFFIX}"
            INTERFACE_INCLUDE_DIRECTORIES "$ENV{CONDA_PREFIX}/include"
        )
    else()
        add_library(${libname} SHARED IMPORTED)
        set_target_properties(${libname} PROPERTIES
            IMPORTED_LOCATION "$ENV{CONDA_PREFIX}/lib/Hammer/${CMAKE_SHARED_LIBRARY_PREFIX}${libname}${CMAKE_SHARED_LIBRARY_SUFFIX}"
            INTERFACE_INCLUDE_DIRECTORIES "$ENV{CONDA_PREFIX}/include"
        )
    endif()
endforeach()

add_library(ThirdParty::YamlCpp INTERFACE IMPORTED)
set_property(TARGET ThirdParty::YamlCpp PROPERTY INTERFACE_LINK_LIBRARIES $ENV{CONDA_PREFIX}/lib/libyaml-cpp.so)

SET(BOOST_INCLUDEDIR $ENV{CONDA_PREFIX}/include)
SET(BOOST_LIBRARYDIR $ENV{CONDA_PREFIX}/lib)
find_package(Boost 1.50.0 QUIET COMPONENTS thread)
if(NOT Boost_FOUND)
	message (FATAL_ERROR "Boost not found. You should install it before installing Hammer.")
else()
	message (STATUS "Found Boost version ${Boost_VERSION}.")
	if(NOT TARGET Boost::headers)
        add_library(Boost::headers INTERFACE IMPORTED)
		set_property(TARGET Boost::headers PROPERTY INTERFACE_INCLUDE_DIRECTORIES ${Boost_INCLUDE_DIR})
	endif()
	if(NOT TARGET Boost::thread)
		add_library(Boost::thread INTERFACE IMPORTED)
		set_property(TARGET Boost::thread PROPERTY INTERFACE_LINK_LIBRARIES ${Boost_LIBRARIES})
	endif()
endif()

if(TRUE)
    add_library(ThirdParty::ROOT INTERFACE IMPORTED)
    set_property(TARGET ThirdParty::ROOT PROPERTY INTERFACE_INCLUDE_DIRECTORIES $ENV{CONDA_PREFIX}/include)
    set_property(TARGET ThirdParty::ROOT PROPERTY INTERFACE_LINK_LIBRARIES -L$ENV{CONDA_PREFIX}/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -Wl,-rpath,$ENV{CONDA_ENV}/lib -pthread -lm -ldl -rdynamic)
endif()
if(TRUE)
    add_library(ThirdParty::HepMC INTERFACE IMPORTED)
    set_property(TARGET ThirdParty::HepMC PROPERTY INTERFACE_INCLUDE_DIRECTORIES $ENV{CONDA_PREFIX}/include)
    set_property(TARGET ThirdParty::HepMC PROPERTY INTERFACE_LINK_LIBRARIES $ENV{CONDA_PREFIX}/lib/libHepMC.so)
endif()

function(new_example name)
    add_executable(${name}.exe
                   ${name}.cc)

    target_compile_options(${name}.exe PRIVATE ${Hammer_CompileOptions})

    set_target_properties(${name}.exe PROPERTIES
                                  POSITION_INDEPENDENT_CODE ${Hammer_PIC})

    set(LIBS HammerCore Amplitudes FormFactors HammerTools HammerBase Rates ThirdParty::YamlCpp Boost::thread)
    if(TRUE)
        set(LIBS ${LIBS} ThirdParty::ROOT)
    endif()
    if(TRUE)
        set(LIBS ${LIBS} ThirdParty::HepMC)
    endif()
    if(OFF)
        add_sanitizers(${name}.exe)
    endif()
    target_link_libraries(${name}.exe PRIVATE ${LIBS})
endfunction()


file(GLOB CustomExamples "${PROJECT_SOURCE_DIR}/*.cc")
foreach(examplesourcefile ${CustomExamples})
  STRING(REGEX REPLACE "(.*)/(.*).cc" "\\2" examplename "${examplesourcefile}")
  MESSAGE(STATUS "Process Analysis ${examplename}")
  new_example(${examplename})
endforeach()
