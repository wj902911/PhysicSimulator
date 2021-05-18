#
# Core of libpsim
# 
aux_source_directory(core/ PSIM_SRC)
lsdir(subs ${CMAKE_CURRENT_LIST_DIR}/core)
foreach(sub ${subs})
	message(STATUS "adding ${sub}")
	aux_source_directory(core/${sub} PSIM_SRC)
endforeach()

#
# Visualizer
# 
aux_source_directory(vis/ PSIM_VIS_SRC)
lsdir(subs ${CMAKE_CURRENT_LIST_DIR}/vis)
foreach(sub ${subs})
	aux_source_directory(vis/${sub} PSIM_VIS_SRC)
endforeach()

#
# Create the library
#
add_library(psim STATIC ${PSIM_SRC} ${PSIM_VIS_SRC} "vis/fluid/FluidHook.h" "vis/fluid/FluidHook.cpp" "vis/fluid/FluidVisualizer.h" "vis/fluid/FluidVisualizer.cpp")
target_link_libraries(psim PRIVATE imgui glad glfw)
target_include_directories(psim PRIVATE "${PROJECT_SOURCE_DIR}/libigl/external/eigen"
                                       "${PROJECT_SOURCE_DIR}/libigl/include")
target_include_directories(psim BEFORE PRIVATE ${CMAKE_CURRENT_LIST_DIR})
#target_link_libraries(psim PUBLIC psim_std_dep)

#target_link_libraries(psim PUBLIC psim_gui_dep)

set_relative_rpath(psim)
