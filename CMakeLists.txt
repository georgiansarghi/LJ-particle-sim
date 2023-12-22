CMAKE_MINIMUM_REQUIRED (VERSION 3.10.2)
PROJECT(Practical-Lab CXX)

IF(${CMAKE_CXX_COMPILER_ID} MATCHES "(GNU|Clang)")
  ADD_COMPILE_OPTIONS(-Wall -Wextra -pedantic)
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES "(GNU|Clang)")

FIND_PACKAGE(MPI REQUIRED COMPONENTS CXX MPICXX)

ADD_EXECUTABLE(main
  world.cpp
  velocityverlet.cpp
  ljpotential.cpp
  observer.cpp
  cell.cpp
  timediscretization.cpp
  main.cpp
  thermostat.cpp
  subdomain.cpp
  )
SET_PROPERTY(TARGET main PROPERTY CXX_STANDARD 98)
TARGET_LINK_LIBRARIES(main PUBLIC MPI::MPI_CXX)
