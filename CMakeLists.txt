cmake_minimum_required(VERSION 3.11)
project(example-project)

set(SRCDIR src)
set(INCDIR include)
set(OBJDIR obj)
set(BIN main)

file(GLOB_RECURSE SRCS "${SRCDIR}/*.cpp")
file(GLOB SRCS "main.cpp")
include_directories(${INCDIR})
add_executable(${BIN} ${SRCS})
