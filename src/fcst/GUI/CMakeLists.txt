#---------------------------------------------------------------------------
#
#    FCST: Fuel Cell Simulation Toolbox
#
#    Copyright (C) 2006-13 by Energy Systems Design Laboratory, University of Alberta
#
#    This software is distributed under the MIT License.
#    For more information, see the README file in /doc/LICENSE
#
#---------------------------------------------------------------------------

####################### 
# 
# OpenFCST GUI CMake File
#

# Developer: Philip Wardlaw, 2014
#
#######################

#Group source
FILE(GLOB srC *.cc contrib/*.cpp)
FILE(GLOB CsrC *.cc source/*.cc)

SET(SRCFiles 
      ${srC}
      ${CsrC}) 

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/contrib/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/source/)

#QT steps:
FIND_PACKAGE(Qt4 REQUIRED QtCore QtGui QtXml)
MARK_AS_ADVANCED(QT_QMAKE_EXECUTABLE)
INCLUDE(${QT_USE_FILE})

QT4_ADD_RESOURCES(SOURCE_RCC
  application.qrc
  )

QT4_WRAP_CPP(SOURCE_MOC
  include/main_window.h
  include/settings_window.h
  contrib/browse_lineedit.h
  contrib/parameter_delegate.h)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPENFCST_C11_FLAG}")

add_definitions(${CMAKE_CXX_FLAGS})


ADD_EXECUTABLE(${FCST_GUI_EXE} ${SRCFiles} 
  ${SOURCE_RCC} 
  ${SOURCE_MOC})


TARGET_LINK_LIBRARIES(${FCST_GUI_EXE} ${QT_LIBRARIES})

INSTALL(TARGETS ${FCST_GUI_EXE}
  RUNTIME DESTINATION ../bin
  COMPONENT parameter_gui
  )
file(GLOB ContribHFiles contrib/*.h)
file(GLOB includeHFiles include/*.h)
INSTALL(FILES ${ContribHFiles}
         DESTINATION include/GUI/contrib
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)  
INSTALL(FILES ${includeHFiles}
         DESTINATION include/GUI/include
         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)   