
cmake_minimum_required(VERSION 3.2)
project(brownian C)

#set(CMAKE_C_COMPILER "icc")
set(CMAKE_C_STANDARD 99)
#set(CMAKE_C_FLAGS "-xAVX -O3")
#set(CMAKE_C_FLAGS_RELEASE  "-DNDEBUG -fp-model fast=2")
set(CMAKE_C_FLAGS "-mavx -O3")
set(CMAKE_C_FLAGS_RELEASE  "-DNDEBUG -fp-model fast=2")


add_executable(brownian602 main.c particle.h data.h pinning.h vertex.h)

target_link_libraries(brownian602 m)
