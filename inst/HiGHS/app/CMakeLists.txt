# create highs binary using library without pic
add_executable(highs)

target_sources(highs PRIVATE RunHighs.cpp)

if (UNIX)
        target_compile_options(highs PUBLIC "-Wno-unused-variable")
        target_compile_options(highs PUBLIC "-Wno-unused-const-variable")
endif()

target_link_libraries(highs libhighs)

if(EMSCRIPTEN AND EMSCRIPTEN_HTML)
        set(CMAKE_EXECUTABLE_SUFFIX ".html")
        set_target_properties(highs PROPERTIES LINK_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/highs_webdemo_shell.html)
endif()

target_include_directories(highs PRIVATE
        $<BUILD_INTERFACE:${HIGHS_SOURCE_DIR}/app>  
        )

# install the binary
install(TARGETS highs EXPORT highs-targets
        RUNTIME)
