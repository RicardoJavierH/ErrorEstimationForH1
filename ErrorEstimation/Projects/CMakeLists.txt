# Create Tools library
add_subdirectory(Tools)

# Created executable targets
add_subdirectory(ErrorNaca)
add_subdirectory(ErrorH1HDiv)
add_subdirectory(ErrorEstimationH1)
add_subdirectory(ErrorEstimationH1HybridH1)
add_subdirectory(ErrorEstimationH1New)
add_subdirectory(ErrorEstimationHDiv)
add_subdirectory(HybridH1)
add_subdirectory(ErrorEstimationMHM)

# This project is marked as optional since it requires OpenCV external library
option(BUILD_MHM_MAZE_PROJECT "Whether to build MHM Maze project" OFF)
if(BUILD_MHM_MAZE_PROJECT)
    add_subdirectory(MHM_Maze)
endif()

# This project is marked as optional since it requires libInterpolate external library
option(BUILD_UNISIM_PROJECT "Whether to build estimation on UNISIM mesh" OFF)
if(BUILD_UNISIM_PROJECT)
    add_subdirectory(UNISIM)
endif()

# This project is marked as optional since it requires libInterpolate external library
option(BUILD_SPE10_PROJECT "Whether to build estimation on SPE10 mesh" OFF)
if(BUILD_SPE10_PROJECT)
    add_subdirectory(SPE10)
endif()

