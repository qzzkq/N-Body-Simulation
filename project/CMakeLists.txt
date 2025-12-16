cmake_minimum_required(VERSION 3.21)
project(NBodySimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Для стандартных путей установки (bin, lib, etc.)
include(GNUInstallDirs)

# ---------------------------------------------------------
# 1. Поиск зависимостей
# ---------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(GLEW   REQUIRED)
find_package(glfw3  REQUIRED)
find_package(glm    REQUIRED)
find_package(HDF5   REQUIRED COMPONENTS C CXX HL)

# ---------------------------------------------------------
# 2. Определение CUDA (Автоматически)
# ---------------------------------------------------------
include(CheckLanguage)
check_language(CUDA)

if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "CUDA found! GPU acceleration enabled.")
    
    add_compile_definitions(USE_CUDA)
    set(CUDA_SOURCES src/barnes_hut_cuda.cu)
else()
    message(STATUS "CUDA NOT found. Building CPU-only version.")
    set(CUDA_SOURCES "")
endif()

# ---------------------------------------------------------
# 3. Настройка системных библиотек (Linux)
# ---------------------------------------------------------
if(UNIX AND NOT WIN32)
    set(SYSTEM_LIBS X11 Xrandr Xi Xxf86vm Xcursor pthread dl)
endif()

# ---------------------------------------------------------
# 4. Сборка ядра (Core Library)
# ---------------------------------------------------------
set(CORE_SRC
    src/object.cpp
    src/renderer.cpp
    src/control.cpp
    src/data.cpp
    src/bodysystem.cpp
    src/barnes_hut.cpp
    ${CUDA_SOURCES}) 

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
    ${SYSTEM_LIBS}
) 

if(CMAKE_CUDA_COMPILER)
    target_link_libraries(Core PUBLIC CUDA::cudart)
    set_property(TARGET Core PROPERTY CUDA_ARCHITECTURES native)
endif()

# ---------------------------------------------------------
# 5. Исполняемые файлы
# ---------------------------------------------------------

# Симуляция
add_executable(Simulate src/main.cpp)
target_link_libraries(Simulate PRIVATE Core)

# Проигрыватель
add_executable(Replay src/replay.cpp)
target_link_libraries(Replay PRIVATE Core)

# Генератор данных
add_executable(DataCreator src/datacreator.cpp)
target_link_libraries(DataCreator PRIVATE Core)

# ---------------------------------------------------------
# 6. Установка (INSTALL)
# ---------------------------------------------------------

# 1. Установка бинарных файлов и библиотеки
install(TARGETS Simulate Replay DataCreator Core
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR} 
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

# 2. Установка папки data
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/data")
    install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_BINDIR}) 
else()
    message(WARNING "Folder 'data' not found in source directory! It won't be installed.")
endif() 