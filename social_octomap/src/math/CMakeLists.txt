SET (octomath_SRCS
  Vector3.cpp
  Quaternion.cpp
  Pose6D.cpp
 )


ADD_LIBRARY( octomath SHARED ${octomath_SRCS})

SET_TARGET_PROPERTIES( octomath PROPERTIES
  VERSION ${SOCIAL_OCTOMAP_VERSION}
  SOVERSION ${SOCIAL_OCTOMAP_SOVERSION}
)

ADD_LIBRARY( octomath-static STATIC ${octomath_SRCS})
SET_TARGET_PROPERTIES(octomath-static PROPERTIES OUTPUT_NAME "octomath")

if(NOT EXISTS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap")
  file(MAKE_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap")
endif()

export(TARGETS octomath octomath-static
  APPEND FILE "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap/social_octomap-targets.cmake")

install(TARGETS octomath octomath-static
  EXPORT social_octomap-targets
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  ${INSTALL_TARGETS_DEFAULT_ARGS}
)
