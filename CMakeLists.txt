cmake_minimum_required(VERSION 3.23)
project(HMM)

set(CMAKE_CXX_STANDARD 11)

add_executable(HMM main.cpp Viterbi.cpp Viterbi.h)
