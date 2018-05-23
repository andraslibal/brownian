
cmake_minimum_required(VERSION 3.2)
project(brownian C)

set(CMAKE_C_COMPILER "icc")
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "-xAVX -O3 -DNDEBUG -fp-model fast=2")
#set(CMAKE_C_FLAGS "-mavx -O3 -DNDEBUG -ffast-math")


add_executable(brownian main.c particle.h data.h pinning.h vertex.h)

target_link_libraries(brownian m)
