cmake_minimum_required(VERSION 3.21)
project(NBodySimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

# ---------------------------------------------------------
# 1. Поиск зависимостей
# ---------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(GLEW   REQUIRED)
find_package(glfw3  REQUIRED)
find_package(glm    REQUIRED)
find_package(HDF5   REQUIRED COMPONENTS C CXX HL)

if(WIN32)
    add_definitions(-DGLEW_STATIC)
endif()

# ---------------------------------------------------------
# 2. Определение CUDA
# ---------------------------------------------------------
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    add_compile_definitions(USE_CUDA)
    file(GLOB CUDA_SOURCES "src/*.cu")
else()
    set(CUDA_SOURCES "")
endif()

# ---------------------------------------------------------
# 3. Сборка ядра (Core Library)
# ---------------------------------------------------------
set(CORE_SRC
    src/object.cpp
    src/renderer.cpp
    src/control.cpp
    src/bodysystem.cpp
    src/data.cpp
    src/physics.cpp
    src/generators.cpp
    src/barnes_hut.cpp
    src/brut_force.cpp
    ${CUDA_SOURCES}
) 

add_library(Core STATIC ${CORE_SRC})

target_include_directories(Core PUBLIC
    include
    src
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
    glm::glm
) 

if(UNIX AND NOT WIN32)
    target_link_libraries(Core PUBLIC X11 Xrandr Xi Xxf86vm Xcursor pthread dl)
endif()

if(CMAKE_CUDA_COMPILER)
    target_link_libraries(Core PUBLIC CUDA::cudart)
    set_property(TARGET Core PROPERTY CUDA_ARCHITECTURES native)
endif()

# ---------------------------------------------------------
# 4. Исполняемые файлы
# ---------------------------------------------------------
add_executable(Simulate src/main.cpp)
target_link_libraries(Simulate PRIVATE Core)

add_executable(Replay src/replay.cpp)
target_link_libraries(Replay PRIVATE Core)

add_executable(DataCreator src/datacreator.cpp)
target_link_libraries(DataCreator PRIVATE Core)

# ---------------------------------------------------------
# 5. Установка (INSTALL)
# ---------------------------------------------------------

# Устанавливаем EXE в bin/
install(TARGETS Simulate Replay DataCreator
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Устанавливаем папку data/ в папку с исполняемыми файлами
install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_BINDIR})

# ---------------------------------------------------------
# 6. Windows DLL Deployment
# ---------------------------------------------------------
if(WIN32)
    # Копирование DLL при установке (cmake --install)
    install(FILES $<TARGET_RUNTIME_DLLS:Simulate> DESTINATION ${CMAKE_INSTALL_BINDIR})
    
    # Копирование папки data и DLL сразу после сборки (для запуска из папки build)
    add_custom_command(TARGET Simulate POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/data"
            "$<TARGET_FILE_DIR:Simulate>/data"
        COMMAND ${CMAKE_COMMAND} -E copy -t "$<TARGET_FILE_DIR:Simulate>"
            $<TARGET_RUNTIME_DLLS:Simulate>
        COMMENT "Copying data and DLLs to executable directory..."
        COMMAND_EXPAND_LISTS
    )
endif()