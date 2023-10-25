project(utils Fortran)

set(UTILS_SOURCES
   arrange.f
   bsn.f
   distr.f
   fetch.f
   grade.f
   inpgrd.f
   make_evt.f
   moist.f
   non_ca.f
   ragmet.f
   rdtmp.f
   read_psm.f
   read_rag.f
   read_tmp_ef.f
   snw.f
   tmp.f
   weight.f
   write_event.f
)
message(STATUS "Sources for UTILS: ${UTILS_SOURCES}")

add_library(${PROJECT_NAME} STATIC ${UTILS_SOURCES})

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION include)

