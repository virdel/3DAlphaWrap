cmake_minimum_required(VERSION 3.15)
project(alpha_wrap CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

add_executable(alpha_wrap src/alpha_wrap.cpp)

if(MSVC)
  target_compile_options(alpha_wrap PRIVATE /W3 /O2 /EHsc)
else()
  target_compile_options(alpha_wrap PRIVATE -O2 -Wall -Wextra)
endif()
