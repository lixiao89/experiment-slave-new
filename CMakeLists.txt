cmake_minimum_required(VERSION 2.8.3)
project(experiment_slave)

## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(catkin REQUIRED COMPONENTS
  cisst_msgs
  cisst_ros_bridge
  roscpp
  std_msgs
)

set (REQUIRED_CISST_LIBRARIES
  cisstCommon
  cisstCommonXML
  cisstVector
  cisstNumerical
  cisstOSAbstraction
  cisstMultiTask
  cisstRobot
  cisstParameterTypes)

find_package(cisst REQUIRED ${REQUIRED_CISST_LIBRARIES} )

catkin_package(
#  INCLUDE_DIRS include
#  LIBRARIES force_estimation
#  CATKIN_DEPENDS cisst_msgs cisst_ros_bridge roscpp std_msgs
#  DEPENDS system_lib
)

include_directories(
  ${catkin_INCLUDE_DIRS}
)

message(STATUS "Enabling profiling with gprof")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pg")
set(CMAKE_SHAREDBoost_USE_MULTITHREADED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
set(CMAKE_EXE_FLAGS "${CMAKE_EXE_FLAGS} -pg") 

if( cisst_FOUND )
  
  include( ${CISST_USE_FILE} )
  
  cisst_find_saw_component( sawKeyboard )
  cisst_find_saw_component( sawControllers )
  cisst_find_saw_component( sawJR3ForceSensor )
  cisst_find_saw_component( sawBarrett )
  cisst_find_saw_component( sawCANBus )
  
  # if( sawKeyboard_FOUND       AND 
  #    sawJR3ForceSensor_FOUND AND 
  #    sawControllers_FOUND    AND
  #    sawBarrett_FOUND        AND
  #    sawCANBus_FOUND ) 

    include_directories( 
     ${sawKeyboard_INCLUDE_DIR} 
     ${sawControllers_INCLUDE_DIR} 
     ${sawJR3ForceSensor_INCLUDE_DIR} 
     ${sawBarrett_INCLUDE_DIR} 
     ${sawCANBus_INCLUDE_DIR}) 

    link_directories( 
     ${sawKeyboard_LIBRARY_DIR} 
     ${sawControllers_LIBRARY_DIR} 
     ${sawJR3ForceSensor_LIBRARY_DIR} 
     ${sawBarrett_LIBRARY_DIR} 
     ${sawCANBus_LIBRARY_DIR}) 
    
     set(
        REQUIRED_SAW_COMPONENTS
        sawCANBus
        sawBarrett
        sawControllers
        sawJR3ForceSensor
        sawKeyboard )


    add_executable( 
      hybridctrl 
      src/main.cpp 
      src/osaHybridForcePosition.cpp 
      src/mtsUserStudySlave.cpp )

    target_link_libraries( 
      hybridctrl
      ${REQUIRED_SAW_COMPONENTS}
      ${catkin_LIBRARIES}
      pthread)

  cisst_target_link_libraries(hybridctrl ${REQUIRED_CISST_LIBRARIES} )
    
  # endif( sawKeyboard_FOUND       AND 
  #       sawJR3ForceSensor_FOUND AND 
  #       sawControllers_FOUND    AND
  #       sawBarrett_FOUND        AND
  #       sawCANBus_FOUND ) 

    add_executable( parsexml src/parsexml.cpp ${HEADER_FILES} )
    cisst_target_link_libraries( parsexml ${REQUIRED_CISST_LIBRARIES} )

endif ( cisst_FOUND )


## Generate messages in the 'msg' folder
# add_message_files(
#   FILES
#   Message1.msg
#   Message2.msg
# )

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
# generate_messages(
#   DEPENDENCIES
#   cisst_msgs#   std_msgs
# )

