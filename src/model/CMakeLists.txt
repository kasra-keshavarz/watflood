project(model Fortran)

set(PROJECT_EXECUTABLE
   CHARM
   )

set(MODEL_SOURCES
   aet.f90
   baseflow.f
   dds_code.f
   disaggregate.f90
   errormsg.f
   etharg_beta.f
   etharg.f
   etin.f
   etpriest.f
   flowreset.f90
   header.f
   ice_factor.f
   intcept.f
   lake_evap.f
   lake_ice.f
   lake_lvl_update.for
   lst.f
   melt.f
   opt.f
   options.f
   options.f90
   outevt.f
   precip_adjust.f
   process_hum.f
   process_rain.f
   process_temp.f
   rdshed.for
   rdsnow.f
   rd_wfo_spec.f
   read_divert.F90
   read_drn.f
   read_dsn.f
   read_flowinit.f
   read_gsm_ef.f
   read_gsn.f
   read_hum.f
   read_level.f90
   read_lvl_old.f
   read_outfiles.f
   read_rain_ef.f
   read_resume.f90
   read_resvin_ef.f
   read_soilinit.f
   read_swe.f
   read_sweinit.f
   read_temp_ef.f
   Rel_hum.f
   rerout.f
   route.f
   rt_channel.f90
   rt_pond.f90
   rt_wetland.f90
   rules_sl.f
   rules_tl.f90
   runof6.f
   sed.f
   sedfl.f
   sedroute.f
   sensitivity.f
   soilinit.f
   stats.f
   sub.f90
   synflw.f
   temper.f
   timer.f
   tracer0.f
   tracer1.f
   tracer2.f
   tracer3.f
   tracer4.f
   tracer5.f
   tracer.f
   tracerGW.f
   tracewet.f
   watbal.f
   wfocode.f
   withdraw.f
   wqnut.f
   wqread.f
   wqroute.f
   wqsed.f
   write_modflow.f
   write_resume.f
   write_river1d_tb0.f
   write_sensitivity.f
   write_soilinit.f
   write_wfo.f
   CHARM.f90
)

message(STATUS "Sources for MODEL: ${MODEL_SOURCES}")

add_library(${PROJECT_NAME} STATIC ${MODEL_SOURCES})
target_link_libraries(
   ${PROJECT_NAME}
   common
   iso
   )

add_executable(
   ${PROJECT_EXECUTABLE}
   CHARM.f90
   )
target_link_libraries(
   ${PROJECT_EXECUTABLE}
   core
   common
   utils
   iso
   model
   )

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION include)
install(TARGETS ${PROJECT_EXECUTABLE}
	RUNTIME DESTINATION bin
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION include)
