SET (social_octomap_SRCS
  AbstractOcTree.cpp
  AbstractOccupancyOcTree.cpp
  Pointcloud.cpp
  ScanGraph.cpp
  CountingOcTree.cpp
  OcTree.cpp
  OcTreeNode.cpp
  OcTreeStamped.cpp
  ColorOcTree.cpp
  )

# dynamic and static libs, see CMake FAQ:
ADD_LIBRARY( social_octomap SHARED ${social_octomap_SRCS})
set_target_properties( social_octomap PROPERTIES
  VERSION ${SOCIAL_OCTOMAP_VERSION}
  SOVERSION ${SOCIAL_OCTOMAP_SOVERSION}
)
ADD_LIBRARY( social_octomap-static STATIC ${social_octomap_SRCS})
SET_TARGET_PROPERTIES(social_octomap-static PROPERTIES OUTPUT_NAME "social_octomap") 
add_dependencies(social_octomap-static octomath-static)

TARGET_LINK_LIBRARIES(social_octomap octomath)

if(NOT EXISTS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap")
  file(MAKE_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap")
endif()

export(TARGETS social_octomap social_octomap-static
  APPEND FILE "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap/social_octomap-targets.cmake")

ADD_SUBDIRECTORY( testing )

ADD_EXECUTABLE(graph2tree graph2tree.cpp)
TARGET_LINK_LIBRARIES(graph2tree social_octomap)

ADD_EXECUTABLE(log2graph log2graph.cpp)
TARGET_LINK_LIBRARIES(log2graph social_octomap)

ADD_EXECUTABLE(binvox2bt binvox2bt.cpp)
TARGET_LINK_LIBRARIES(binvox2bt social_octomap)

ADD_EXECUTABLE(bt2vrml bt2vrml.cpp)
TARGET_LINK_LIBRARIES(bt2vrml social_octomap)

ADD_EXECUTABLE(edit_octree edit_octree.cpp)
TARGET_LINK_LIBRARIES(edit_octree social_octomap)

ADD_EXECUTABLE(convert_octree convert_octree.cpp)
TARGET_LINK_LIBRARIES(convert_octree social_octomap)

ADD_EXECUTABLE(eval_octree_accuracy eval_octree_accuracy.cpp)
TARGET_LINK_LIBRARIES(eval_octree_accuracy social_octomap)

ADD_EXECUTABLE(compare_octrees compare_octrees.cpp)
TARGET_LINK_LIBRARIES(compare_octrees social_octomap)

ADD_EXECUTABLE(simple_example simple_example.cpp)
TARGET_LINK_LIBRARIES(simple_example social_octomap)

ADD_EXECUTABLE(normals_example normals_example.cpp)
TARGET_LINK_LIBRARIES(normals_example social_octomap)

ADD_EXECUTABLE(intersection_example intersection_example.cpp)
TARGET_LINK_LIBRARIES(intersection_example social_octomap)

ADD_EXECUTABLE(octree2pointcloud octree2pointcloud.cpp)
TARGET_LINK_LIBRARIES(octree2pointcloud social_octomap)

install(TARGETS social_octomap social_octomap-static
  EXPORT social_octomap-targets
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  ${INSTALL_TARGETS_DEFAULT_ARGS}
)
install(EXPORT social_octomap-targets DESTINATION "${CMAKE_INSTALL_DATADIR}/social_octomap")

install(TARGETS
	graph2tree
	log2graph
	binvox2bt
	bt2vrml
	edit_octree
	convert_octree
	eval_octree_accuracy
	compare_octrees
	${INSTALL_TARGETS_DEFAULT_ARGS}
)

