set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosOversetApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/external_includes )

## generate variables with the sources
set( KRATOS_OVERSET_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/overset_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/OversetCondition.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/overset_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_overset_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_overset_strategies_to_python.cpp
)

#find SplitTreeSearch library
find_library(SplitTreeSearch_LIBRARY SplitTreeSearch ${SplitTreeSearch_DIR}/lib)

if(SplitTreeSearch_LIBRARY STREQUAL "SplitTreeSearch_LIBRARY-NOTFOUND")
  	message (FATAL_ERROR "SplitTreeSearch library not found")
else(SplitTreeSearch_LIBRARY STREQUAL "SplitTreeSearch_LIBRARY-NOTFOUND")
	message ("-- Found SplitTreeSearch: " ${SplitTreeSearch_LIBRARY})
endif(SplitTreeSearch_LIBRARY STREQUAL "SplitTreeSearch_LIBRARY-NOTFOUND")


###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosOversetApplication SHARED ${KRATOS_OVERSET_APPLICATION_SOURCES})
target_link_libraries(KratosOversetApplication KratosCore KratosTrilinosApplication ${MPI_LIBRARIES} ${SplitTreeSearch_LIBRARY} ${TRILINOS_LIBRARIES} )
set_target_properties(KratosOversetApplication PROPERTIES PREFIX "")
install(TARGETS KratosOversetApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosOversetApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/OversetApplication.py" DESTINATION KratosMultiphysics )
