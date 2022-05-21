cmake_minimum_required(VERSION 3.5)

project(momoko LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (EMSCRIPTEN)
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
  add_compile_options(-pthread)
  add_link_options(--emrun -pthread -sPTHREAD_POOL_SIZE=5 -sPROXY_TO_PTHREAD -sALLOW_MEMORY_GROWTH)
endif()

option(USE_BOOST_QUEUE "Use boost queue instead of custom queue." ON)
add_library(momoko
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
  src/bit_matrix.hpp
  src/bit_matrix.cpp
  src/knuth_yao_sampler.hpp
  src/knuth_yao_sampler.cpp
  src/pksystem.hpp src/pksystem.cpp
  src/pke.h src/pke.cpp
  src/pke_blwe.h src/pke_blwe.cpp
  src/background_sampler.h src/background_sampler.cpp
  src/background_sampler_lockfree.h src/background_sampler_lockfree.cpp
  src/signature.h src/signature.cpp
  src/lockfree_cache_queue.h src/lockfree_cache_queue.cpp
  )
add_executable(momoko_test src/main.cpp src/benchmark_counter.h src/benchmark_counter.cpp)
target_link_libraries(momoko_test PUBLIC momoko)
if(USE_BOOST_QUEUE)
find_package(Boost REQUIRED)
target_include_directories(momoko PUBLIC ${Boost_INCLUDE_DIRS})
target_link_directories(momoko PUBLIC ${Boost_LIBRARIES})
target_compile_definitions(momoko PUBLIC BOOST_QUEUE_ENABLED)
endif()
