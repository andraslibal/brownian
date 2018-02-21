cmake_minimum_required(VERSION 3.5)
project(brownian C)

set(CMAKE_C_STANDARD 99)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")

add_executable(brownian main.c thpool.c)

target_link_libraries(brownian m)