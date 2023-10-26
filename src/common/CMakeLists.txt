project(common Fortran)

set(COMMON_SOURCES
   date_time.f
   EF_Module.f
   EF_ParseUtilities.f
   filecode.f
   flowinit.f
   gregorian_calendar.f90
   julian_day.f
   par_assign.f
   par_init.f
   rdpar.f
   rdsdc.f
   read_evt.f
   read_flow_ef.F90
   read_ice_factor.f
   read_nca.f90
   read_par.f
   read_par_parser.f90
   read_pt2.f
   read_r2c.f
   read_r2s_beta.f
   read_r2s.f
   read_resv_ef.f
   read_shed_ef.f
   read_shed_hype.f90
   read_static_r2c.f
   read_tb0.f90
   read_ts5.f
   write_diversion_tb0.f
   write_flowinit.f
   write_lzsinit.f
   write_par_10.f
   write_par.f
   write_r2c.f
   write_r2s.f
   write_tb0.f
)

message(STATUS "Sources for COMMON: ${COMMON_SOURCES}")

add_library(${PROJECT_NAME} STATIC ${COMMON_SOURCES})

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION include)

