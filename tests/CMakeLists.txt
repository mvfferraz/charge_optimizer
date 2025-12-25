# Simple test framework
add_executable(test_basic test_basic.cpp)
target_link_libraries(test_basic PRIVATE Eigen3::Eigen)

enable_testing()
add_test(NAME BasicTest COMMAND test_basic)
