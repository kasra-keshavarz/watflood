project(iso Fortran)

set(ISO_SOURCES
   craig_gordon_2H.f
   craig_gordon.f
   ISO2Hlake.f
   ISO2Hland.f
   ISO2Hriver.f
   ISO2Hsnow.f
   ISO2Hwetland.f
   ISOlake.f
   ISOland.f
   ISOriver.f
   ISOsnow.f
   isotopes.f
   ISOwetland.f
)

message(STATUS "Sources for CORE: ${ISO_SOURCES}")

add_library(${PROJECT_NAME} STATIC ${ISO_SOURCES})

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION include)

