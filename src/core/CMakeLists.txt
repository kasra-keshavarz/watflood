project(area Fortran)

set(CORE_SOURCES
   area17.f
   Areacg.f
   area_debug.f90
   area_watflood.f
)
message(STATUS "Sources for CORE: ${CORE_SOURCES}")

add_library(${PROJECT_NAME} STATIC ${CORE_SOURCES})

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION include)

