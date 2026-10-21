foreach(example geometry_demo bracket_assembly solver_diagnostics)
  add_executable(${example} ${example}.cpp)
  target_link_libraries(${example} PRIVATE geoconstraint geoconstraint_warnings)
  # Each example exits non-zero if its own self-checks fail, so running it is a test.
  add_test(NAME example_${example} COMMAND ${example})
endforeach()
