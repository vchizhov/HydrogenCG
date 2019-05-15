cmake_minimum_required(VERSION 3.0)

project(HydrogenCG_Part0)


file(GLOB sources *.cpp *.hpp)
add_executable (main ${sources})
target_compile_features(main PUBLIC cxx_std_17)