set(GLM_VERSION "0.9.9")
set(GLM_INCLUDE_DIRS "C:/Users/A1234/Documents/INF584/MyRenderer_basic_raytracer/External/glm")

if (NOT CMAKE_VERSION VERSION_LESS "3.0")
    include("${CMAKE_CURRENT_LIST_DIR}/glmTargets.cmake")
endif()
