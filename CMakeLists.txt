cmake_minimum_required(VERSION 3.5)

project(momoko LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Boost REQUIRED)

add_executable(momoko
  src/main.cpp
  src/ideal_lattice.hpp
  src/ideal_lattice.cpp
  src/ideal_lattice_element.hpp
  src/ideal_lattice_element.cpp
  src/tools.hpp
  src/tools.cpp
  src/gaussian_dist_sampler.hpp
  src/gaussian_dist_sampler.cpp
  src/bernoulli_sampler.hpp
  src/bernoulli_sampler.cpp
  src/CDT_sampler.hpp
  src/CDT_sampler.cpp
  src/descrete_gaussian_calculator.hpp
  src/descrete_gaussian_calculator.cpp
  src/bit_matrix.hpp
  src/bit_matrix.cpp
  src/knuth_yao_sampler.hpp
  src/knuth_yao_sampler.cpp
  src/pksystem.hpp src/pksystem.cpp
  src/pke.h src/pke.cpp
  src/background_sampler.h src/background_sampler.cpp
  src/background_sampler_lockfree.h src/background_sampler_lockfree.cpp
  src/signature.h src/signature.cpp
  )

target_include_directories(momoko PUBLIC ${Boost_INCLUDE_DIRS})
target_link_directories(momoko PUBLIC ${Boost_LIBRARIES})
