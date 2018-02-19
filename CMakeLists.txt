cmake_minimum_required(VERSION 3.9)
project(brownian C)

set(CMAKE_C_STANDARD 99)

add_executable(brownian main.c particle.h)

target_link_libraries(brownian m)