cmake_minimum_required(VERSION 3.21)
project(NBodySimulation LANGUAGES CXX VERSION 1.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---------------------------------------------------------
# Выходные директории
# ---------------------------------------------------------
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# ---------------------------------------------------------
# RPATH (Linux и MinGW)
# ---------------------------------------------------------
# включаем использование относительных rpath
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)

# Для Linux — lib лежит на уровень выше bin
if(UNIX AND NOT WIN32)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib")
    set(CMAKE_BUILD_RPATH   "$ORIGIN/../lib")
endif()

# Для MSYS2 / MinGW64 — rpath не используется, DLL ищутся автоматически через PATH
if(MINGW)
    # Ничего не делаем — MSYS2 сам подхватывает DLL из /mingw64/bin
endif()

# ---------------------------------------------------------
# Dependencies
# ---------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(GLEW   REQUIRED)
find_package(glfw3  REQUIRED)
find_package(glm    REQUIRED)
find_package(HDF5   REQUIRED COMPONENTS C CXX HL)

# Linux system libs
if(UNIX AND NOT WIN32)
    set(SYSTEM_LIBS X11 Xrandr Xi Xxf86vm Xcursor pthread)
endif()

# ---------------------------------------------------------
# Core library
# ---------------------------------------------------------
set(CORE_SRC
    src/object.cpp
    src/renderer.cpp
    src/control.cpp
    src/data.cpp
    src/bodysystem.cpp
    src/barnes_hut.cpp
)

add_library(Core ${CORE_SRC})

target_include_directories(Core PUBLIC
    include
    ${HDF5_INCLUDE_DIRS}
    ${GLEW_INCLUDE_DIRS}
    ${glm_INCLUDE_DIRS}
)

target_link_libraries(Core PUBLIC
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    OpenGL::GL
    GLEW::GLEW
    glfw
    ${SYSTEM_LIBS}
)

# ---------------------------------------------------------
# Executables
# ---------------------------------------------------------
add_executable(Simulate     src/main.cpp)
add_executable(DataCreator  src/datacreator.cpp)

foreach(tgt Simulate DataCreator)
    target_link_libraries(${tgt} PRIVATE Core)
endforeach()

# ---------------------------------------------------------
# Install rules
# ---------------------------------------------------------
install(TARGETS Simulate DataCreator Core
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

install(DIRECTORY include/ DESTINATION include)

install(DIRECTORY data/ DESTINATION bin/data)

