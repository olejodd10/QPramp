add_library(linalg OBJECT vector.c matrix.c)
add_library(qpramp STATIC iterable_set.c qp_ramp.c)

add_executable(qp_ramp_example main.c csv.c timing.c lti.c)

target_link_libraries(linalg PRIVATE compiler_flags)
target_link_libraries(qpramp PRIVATE linalg compiler_flags)
target_link_libraries(qp_ramp_example PRIVATE qpramp linalg compiler_flags)

target_include_directories(linalg PUBLIC ${QPramp_SOURCE_DIR}/include)
target_include_directories(qpramp PUBLIC ${QPramp_SOURCE_DIR}/include)
