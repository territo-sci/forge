include_directories(${CMAKE_SOURCE_DIR}/inc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(Boost_INCLUDE_DIR ${ForgeApp_SOURCE_DIR}/../common_external/boost_1_53_0)
set(Boost_LIBRARY_DIR ${ForgeApp_SOURCE_DIR}/../common_external/boost_1_53_0/stage/lib)
find_package(Boost COMPONENTS system filesystem REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_executable(ForgeApp
               ForgeApp.cpp
               Forge.cpp)

target_link_libraries(ForgeApp ${Boost_LIBRARIES})
