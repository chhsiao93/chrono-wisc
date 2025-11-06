# ------------------------------------------------------------------------------
# List of demos
# ------------------------------------------------------------------------------

set(DEMOS ${DEMOS}
   demo_FSI_DamBreak
   demo_FSI_Compressibility
   demo_FSI_Poiseuille_flow
   demo_FSI_SingleWheelTest
   demo_FSI_BCE
   demo_FSI_Kernels
   demo_FSI_ObjectDrop
   demo_FSI_CylindricalTank
   demo_FSI_WaveTank
   demo_FSI_EmbeddedObject
   demo_FSI_BaffleFlow
   demo_FSI_CouetteFlow
   demo_FSI_Flap
   demo_FSI_AngleRepose
   demo_FSI_Cratering
   benchmark_Acc_FSI_Cratering
   benchmark_Acc_FSI_ConePenetration
   benchmark_Acc_FSI_ConePenetrometer_GRC1
   benchmark_Acc_FSI_NormalBevameter_GRC1
   benchmark_Acc_FSI_Viper_SingleWheel
   benchmark_Acc_FSI_Rassor_SingleDrum
   GNS_benchmark_Acc_FSI_Cratering
   GNS_demo_FSI_AngleRepose
)

if(CH_ENABLE_MODULE_PARDISO_MKL)
   set(DEMOS ${DEMOS}
      demo_FSI_FEAdirections
      demo_FSI_Flexible_Plate
      demo_FSI_Flexible_Cable
      demo_FSI_Flexible_Toroidal_Tire
   )
endif()

#--------------------------------------------------------------

list(APPEND LIBS Chrono_core)
list(APPEND LIBS "Chrono_fsisph")

if(CH_ENABLE_MODULE_VSG)
    list(APPEND LIBS "Chrono_fsisph_vsg")
endif()

if(CH_ENABLE_MODULE_PARDISO_MKL)
    list(APPEND LIBS "Chrono_pardisomkl")
endif()

#--------------------------------------------------------------

message(STATUS "Add demo programs for FSI module")

# Find yaml-cpp for GNS_benchmark_Acc_FSI_Cratering
find_package(yaml-cpp QUIET)
if(NOT yaml-cpp_FOUND)
    find_library(YAML_CPP_LIBRARY NAMES yaml-cpp)
    if(YAML_CPP_LIBRARY)
        set(yaml-cpp_FOUND TRUE)
        message(STATUS "Found yaml-cpp library: ${YAML_CPP_LIBRARY}")
    endif()
endif()

foreach(PROGRAM ${DEMOS})
    add_executable(${PROGRAM} ${PROGRAM}.cpp)
    source_group("" FILES ${PROGRAM}.cpp)

    if(APPLE)
      set_target_properties(${PROGRAM} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${MACOSX_BUNDLE_INFO_PLIST})
    endif()

    if(MSVC)
      set_target_properties(${PROGRAM} PROPERTIES MSVC_RUNTIME_LIBRARY ${CH_MSVC_RUNTIME_LIBRARY})
    endif()

    set_property(TARGET ${PROGRAM} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "$<TARGET_FILE_DIR:${PROGRAM}>")
    target_link_libraries(${PROGRAM} ${LIBS})

    # Add yaml-cpp for GNS_benchmark_Acc_FSI_Cratering
    if(${PROGRAM} STREQUAL "GNS_benchmark_Acc_FSI_Cratering" AND yaml-cpp_FOUND)
        if(TARGET yaml-cpp)
            target_link_libraries(${PROGRAM} yaml-cpp)
        else()
            target_link_libraries(${PROGRAM} ${YAML_CPP_LIBRARY})
        endif()
        message(STATUS "Linking yaml-cpp to ${PROGRAM}")
    endif()

    install(TARGETS ${PROGRAM} DESTINATION ${CH_INSTALL_DEMO})
endforeach()
