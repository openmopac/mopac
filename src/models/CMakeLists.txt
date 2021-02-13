target_sources(mopac-all PRIVATE
# modules
        ${CMAKE_CURRENT_SOURCE_DIR}/journal_references_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_AM1_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_AM1_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_INDO_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM3_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM3_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM6_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM6_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM7_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM7_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM7_TS_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_PM8_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_RM1_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_RM1_Sparkles_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_mndo_C.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/parameters_for_mndod_C.F90
# source
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_am1.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_mndo.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_mndod.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/alpb_and_xfac_pm3.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/calpar.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/refer.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/sp_two_electron.F90
        ${CMAKE_CURRENT_SOURCE_DIR}/switch.F90
)
