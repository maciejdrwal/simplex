cmake_minimum_required(VERSION 3.14)

project(Simplex VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(SOURCE_PATH "${PROJECT_SOURCE_DIR}/src")
set(INCLUDE_PATH "${PROJECT_SOURCE_DIR}/include")
file(GLOB_RECURSE HPPS "${SOURCE_PATH}/*.h")
file(GLOB_RECURSE CPPS "${SOURCE_PATH}/*.cpp")
list(REMOVE_ITEM CPPS "${SOURCE_PATH}/main.cpp")
						  
add_library(${PROJECT_NAME}_objects OBJECT ${CPPS} ${HPPS})
add_executable(${PROJECT_NAME} $<TARGET_OBJECTS:${PROJECT_NAME}_objects> "${SOURCE_PATH}/main.cpp")

target_include_directories(${PROJECT_NAME}_objects PUBLIC "${INCLUDE_PATH}")
target_include_directories(${PROJECT_NAME} PUBLIC "${INCLUDE_PATH}")
target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_17)

# Build UnitTests

include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

set(TEST_SOURCE_PATH ${PROJECT_SOURCE_DIR}/tests)
add_executable(UnitTests ${TEST_SOURCE_PATH}/UnitTests.cpp $<TARGET_OBJECTS:${PROJECT_NAME}_objects>)

target_include_directories(UnitTests PUBLIC "${INCLUDE_PATH}")
target_compile_features(UnitTests PRIVATE cxx_std_17)


target_link_libraries(
  UnitTests
  ${PROJECT_NAME}_objects
  GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(UnitTests)
