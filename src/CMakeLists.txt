include_directories(${CMAKE_SOURCE_DIR}/inc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

add_executable(ForgeApp
               ForgeApp.cpp
               Forge.cpp
               BricksHeader.cpp)

