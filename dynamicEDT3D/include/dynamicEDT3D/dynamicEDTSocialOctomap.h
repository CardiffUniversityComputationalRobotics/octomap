/**
 * dynamicEDT3D:
 * A library for incrementally updatable Euclidean distance transforms in 3D.
 * @author C. Sprunk, B. Lau, W. Burgard, University of Freiburg, Copyright (C) 2011.
 * @see http://social_octomap.sourceforge.net/
 * License: New BSD License
 */

/*
 * Copyright (c) 2011-2012, C. Sprunk, B. Lau, W. Burgard, University of Freiburg
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University of Freiburg nor the names of its
 *       contributors may be used to endorse or promote products derived from
 *       this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef DYNAMICEDTSOCIAL_OCTOMAP_H_
#define DYNAMICEDTSOCIAL_OCTOMAP_H_

#include "dynamicEDT3D.h"
#include <social_octomap/OcTree.h>
#include <social_octomap/OcTreeStamped.h>

/// A DynamicEDTSocialOctomapBase object connects a DynamicEDT3D object to an social_octomap.
template <class TREE>
class DynamicEDTSocialOctomapBase : private DynamicEDT3D
{
public:
	/** Create a DynamicEDTSocialOctomapBase object that maintains a distance transform in the bounding box given by bbxMin, bbxMax and clamps distances at maxdist.
	 *  treatUnknownAsOccupied configures the treatment of unknown cells in the distance computation.
	 *
	 *  The constructor copies occupancy data but does not yet compute the distance map. You need to call udpate to do this.
	 *
	 *  The distance map is maintained in a full three-dimensional array, i.e., there exists a float field in memory for every voxel inside the bounding box given by bbxMin and bbxMax. Consider this when computing distance maps for large social_octomaps, they will use much more memory than the social_octomap itself!
	 */
	DynamicEDTSocialOctomapBase(float maxdist, TREE *_octree, social_octomap::point3d bbxMin, social_octomap::point3d bbxMax, bool treatUnknownAsOccupied);

	virtual ~DynamicEDTSocialOctomapBase();

	/// trigger updating of the distance map. This will query the social_octomap for the set of changes since the last update.
	/// If you set updateRealDist to false, computations will be faster (square root will be omitted), but you can only retrieve squared distances
	virtual void update(bool updateRealDist = true);

	/// retrieves distance and closestObstacle (closestObstacle is to be discarded if distance is maximum distance, the method does not write closestObstacle in this case).
	/// Returns DynamicEDTSocialOctomapBase::distanceValue_Error if point is outside the map.
	void getDistanceAndClosestObstacle(const social_octomap::point3d &p, float &distance, social_octomap::point3d &closestObstacle) const;

	/// retrieves distance at point. Returns DynamicEDTSocialOctomapBase::distanceValue_Error if point is outside the map.
	float getDistance(const social_octomap::point3d &p) const;

	/// retrieves distance at key. Returns DynamicEDTSocialOctomapBase::distanceValue_Error if key is outside the map.
	float getDistance(const social_octomap::OcTreeKey &k) const;

	/// retrieves squared distance in cells at point. Returns DynamicEDTSocialOctomapBase::distanceInCellsValue_Error if point is outside the map.
	int getSquaredDistanceInCells(const social_octomap::point3d &p) const;

	// variant of getDistanceAndClosestObstacle that ommits the check whether p is inside the area of the distance map. Use only if you are certain that p is covered by the distance map and if you need to save the time of the check.
	void getDistanceAndClosestObstacle_unsafe(const social_octomap::point3d &p, float &distance, social_octomap::point3d &closestObstacle) const;

	// variant of getDistance that ommits the check whether p is inside the area of the distance map. Use only if you are certain that p is covered by the distance map and if you need to save the time of the check.
	float getDistance_unsafe(const social_octomap::point3d &p) const;

	// variant of getDistance that ommits the check whether p is inside the area of the distance map. Use only if you are certain that p is covered by the distance map and if you need to save the time of the check.
	float getDistance_unsafe(const social_octomap::OcTreeKey &k) const;

	// variant of getSquaredDistanceInCells that ommits the check whether p is inside the area of the distance map. Use only if you are certain that p is covered by the distance map and if you need to save the time of the check.
	int getSquaredDistanceInCells_unsafe(const social_octomap::point3d &p) const;

	/// retrieve maximum distance value
	float getMaxDist() const
	{
		return maxDist * octree->getResolution();
	}

	/// retrieve squared maximum distance value in grid cells
	int getSquaredMaxDistCells() const
	{
		return maxDist_squared;
	}

	/// Brute force method used for debug purposes. Checks occupancy state consistency between social_octomap and internal representation.
	bool checkConsistency() const;

	/// distance value returned when requesting distance for a cell outside the map
	static float distanceValue_Error;
	/// distance value returned when requesting distance in cell units for a cell outside the map
	static int distanceInCellsValue_Error;

private:
	void initializeOcTree(social_octomap::point3d bbxMin, social_octomap::point3d bbxMax);
	void insertMaxDepthLeafAtInitialize(social_octomap::OcTreeKey key);
	void updateMaxDepthLeaf(social_octomap::OcTreeKey &key, bool occupied);

	void worldToMap(const social_octomap::point3d &p, int &x, int &y, int &z) const;
	void mapToWorld(int x, int y, int z, social_octomap::point3d &p) const;
	void mapToWorld(int x, int y, int z, social_octomap::OcTreeKey &key) const;

	TREE *octree;
	bool unknownOccupied;
	int treeDepth;
	double treeResolution;
	social_octomap::OcTreeKey boundingBoxMinKey;
	social_octomap::OcTreeKey boundingBoxMaxKey;
	int offsetX, offsetY, offsetZ;
};

typedef DynamicEDTSocialOctomapBase<social_octomap::OcTree> DynamicEDTSocialOctomap;
typedef DynamicEDTSocialOctomapBase<social_octomap::OcTreeStamped> DynamicEDTSocialOctomapStamped;

#include "dynamicEDTSocialOctomap.hxx"

#endif /* DYNAMICEDTSOCIAL_OCTOMAP_H_ */
