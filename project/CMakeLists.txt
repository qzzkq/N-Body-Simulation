cmake_minimum_required(VERSION 3.21)
project(NBodySimulation LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(DEFINED CMAKE_TOOLCHAIN_FILE AND CMAKE_TOOLCHAIN_FILE MATCHES "vcpkg")
    set(X_VCPKG_APPLOCAL_DEPS_INSTALL ON CACHE BOOL "Automatically copy dependencies on install")
endif()

include(GNUInstallDirs)

# По умолчанию ВЫКЛЮЧЕНО, чтобы vcpkg сам копировал DLL без ошибок доступа
option(USE_MANUAL_DLL_COPY "Включить ручное копирование DLL после сборки" OFF)

# ---------------------------------------------------------
# 1. Поиск зависимостей
# ---------------------------------------------------------
find_package(OpenGL REQUIRED)
find_package(GLEW   REQUIRED)
find_package(glfw3  REQUIRED)
find_package(glm    REQUIRED)
find_package(HDF5   REQUIRED COMPONENTS C CXX HL)
find_package(OpenMP REQUIRED)

if(WIN32)
    add_compile_definitions(NOMINMAX)
endif()

# ---------------------------------------------------------
# 2. Определение CUDA
# ---------------------------------------------------------
if(WIN32 AND NOT MSVC)
    message(WARNING "Сборка с CUDA на Windows требует MSVC (Visual Studio). MinGW/MSYS2 не поддерживается. Поддержка CUDA отключена.")
    set(CUDA_SOURCES "")
else()
    include(CheckLanguage)
    check_language(CUDA)

    if(CMAKE_CUDA_COMPILER)
        enable_language(CUDA)
        find_package(CUDAToolkit REQUIRED)
        add_compile_definitions(USE_CUDA)
        file(GLOB CUDA_SOURCES CONFIGURE_DEPENDS "src/*.cu" "algo/*.cu")
    else()
        message(STATUS "CUDA компилятор не найден. Сборка продолжится без поддержки GPU.")
        set(CUDA_SOURCES "")
    endif()
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
    algo/barnes_hut.cpp
    algo/brut_force.cpp
    algo/whfast.cpp
    ${CUDA_SOURCES}
) 

add_library(Core STATIC ${CORE_SRC})

target_include_directories(Core 
    PUBLIC 
        include 
        algo 
        ${HDF5_INCLUDE_DIRS}
        ${GLEW_INCLUDE_DIRS}
        ${glm_INCLUDE_DIRS}
    PRIVATE 
        src  
)

target_link_libraries(Core PUBLIC
    ${HDF5_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
    OpenGL::GL
    GLEW::GLEW
    glfw
    glm::glm
    OpenMP::OpenMP_CXX
)

if(UNIX AND NOT WIN32)
    target_link_libraries(Core PUBLIC X11 Xrandr Xi Xxf86vm Xcursor pthread dl)
endif()

if(CMAKE_CUDA_COMPILER)
    target_link_libraries(Core PUBLIC CUDA::cudart)
    set_property(TARGET Core PROPERTY CUDA_ARCHITECTURES "80;86;89")
    # Fast-math для CUDA: rsqrtf/sqrtf/etc. — hardware-fast, без strict IEEE.
    # Безопасно для нашего force-kernel (точность всё равно FP32 внутри FP64 акк.).
    target_compile_options(Core PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>
    )
endif()

if(MSVC)
    target_compile_options(Core PRIVATE 
        $<$<COMPILE_LANGUAGE:CXX>:/Zc:preprocessor>
        $<$<COMPILE_LANGUAGE:CXX>:-openmp:llvm>
        $<$<COMPILE_LANGUAGE:CXX>:-openmp:experimental>
    )
    
    if(CMAKE_CUDA_COMPILER)
        target_compile_options(Core PRIVATE 
            $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
            $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
        )
        target_compile_definitions(Core PRIVATE 
            $<$<COMPILE_LANGUAGE:CUDA>:CCCL_IGNORE_MSVC_TRADITIONAL_PREPROCESSOR_WARNING>
        )
    endif()
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=generic")
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

install(TARGETS Simulate Replay DataCreator
    RUNTIME DESTINATION build_cmake
)

install(DIRECTORY data DESTINATION build_cmake)

if(WIN32)
    set(VCPKG_BIN_DIR "C:/vcpkg/installed/x64-windows/bin")

    install(FILES 
        "${VCPKG_BIN_DIR}/hdf5.dll"
        "${VCPKG_BIN_DIR}/hdf5_cpp.dll"
        "${VCPKG_BIN_DIR}/glew32.dll"
        "${VCPKG_BIN_DIR}/glfw3.dll"
        "${VCPKG_BIN_DIR}/zlib1.dll"
        "${VCPKG_BIN_DIR}/szip.dll"
        DESTINATION build_cmake
        OPTIONAL
    )
endif()

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/gui")
    install(DIRECTORY gui DESTINATION .)
endif()

# ---------------------------------------------------------
# 6. Windows DLL Deployment
# ---------------------------------------------------------
if(WIN32)
    add_custom_command(TARGET Simulate POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_directory
            "${CMAKE_CURRENT_SOURCE_DIR}/data"
            "$<TARGET_FILE_DIR:Simulate>/data"
        COMMENT "Copying data folder to Simulate directory..."
    )

    if(USE_MANUAL_DLL_COPY)
        install(FILES $<TARGET_RUNTIME_DLLS:Simulate> DESTINATION ${CMAKE_INSTALL_BINDIR})
        
        foreach(TARGET_NAME Simulate Replay DataCreator)
            add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy -t "$<TARGET_FILE_DIR:${TARGET_NAME}>"
                    $<TARGET_RUNTIME_DLLS:${TARGET_NAME}>
                COMMENT "Copying DLLs to ${TARGET_NAME} directory (Manual Mode)..."
                COMMAND_EXPAND_LISTS
            )
        endforeach()
    endif()
endif()

include(CPack)

# ---------------------------------------------------------
# 7. Тесты
# ---------------------------------------------------------
#enable_testing()
#add_subdirectory(tests)