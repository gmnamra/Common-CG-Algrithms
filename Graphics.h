// Copyright (C) Common Graphics Algorithms e.U, Fang Hao
//
// This file is implementation of Common Graphics Algorithms.
//
// Please contact the author if any conditions of this file are
// not clear to you.
//
// Author: Fang Hao .Nanjing University ,VISG

#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <assert.h>
#include <glm.hpp>
#include <Fade_2D.h>

using namespace std;
using namespace GEOM_FADE25D;

/********************** Declare All Functions Here ***********************/

double Cos(glm::vec3 vec_1, glm::vec3 vec_2);
double Sin(glm::vec3 vec_1, glm::vec3 vec_2);
double Tan(glm::vec3 vec_1, glm::vec3 vec_2);

double pointToLineDistance(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);
double pointToPlaneDistance(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);

glm::vec3 pointToPlaneProjection(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
glm::vec3 pointToLineProjection(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);

glm::vec3 lineToLineIntersection(glm::vec3 l1_1, glm::vec3 l1_2, glm::vec3 l2_1, glm::vec3 l2_2);
glm::vec3 lineToPlaneIntersection(glm::vec3 l1, glm::vec3 l2, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);

bool isPointOnLine(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);
bool isPointOnSegment(glm::vec3 p, glm::vec3 s1, glm::vec3 s2);
bool isPointInPlane(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
bool isPointInTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);


//*************************** Cout Test **********************************
/*
* Cout glm::vector in standard form: vec2,vec3,vec4
*/
void Out(glm::vec2 vec_2)
{
	cout << "Vec2 -- X : " << vec_2.x << "  Y : " << vec_2.y << endl;
}

void Out(glm::vec3 vec_3)
{
	cout << "Vec3 -- X : " << vec_3.x << "  Y : " << vec_3.y << "  Z : " << vec_3.z  << endl;
}

void Out(glm::vec4 vec_4)
{
	cout << "Vec4 -- X : " << vec_4.x << "  Y : " << vec_4.y << "  Z : " << vec_4.z << " A : " << vec_4.a << endl;
}

//*************************** Delaunay ***********************************

/* \brief 2D Delaunay Triangulation
*
* Delaunay triangulation using exist contour points
*
* @param points are the exist contour points to delaunay
* @param segments are the constraint condition in delaunay
*  most of them are line segments
* @param is_loop specifies the line segments whether to loop
* @param max_length and min_length are the max and min limit
*  length of the segments that delaunay generated,default 0
* @param path is the path to save obj
* @param name is the name of obj
*
*	This method is to Delaunay Triangulation
*/

bool Delaunay(vector<glm::vec2> points, vector<glm::vec2> segments, bool is_loop,
	string name,double max_length = 0, double min_length = 0,string path = "defualt")
{
	// read the contour points
	Fade_2D dt;
	vector<Point2> vPoints;
	for (int i = 0; i < points.size(); ++i)
	{
		vPoints.push_back(Point2(points[i].x, points[i].y, 0.0f));
	}
	dt.insert(vPoints);
	// read the segment points
	vector<Point2> sPoints;
	for (int i = 0; i < segments.size(); ++i)
	{
		sPoints.push_back(Point2(segments[i].x, segments[i].y, 0.0f));
	}

	// read the constraint segments
	vector<Segment2> vSegments;
	if (is_loop == true)
	{
		for (int i = 0; i < sPoints.size(); ++i)
		{
			Point2& p0(sPoints[i]);
			Point2& p1(sPoints[(i + 1) % sPoints.size()]);
			vSegments.push_back(Segment2(p0, p1));
		}
	}
	else
	{
		for (int i = 1; i < sPoints.size(); ++i)
		{
			Point2& p0(sPoints[i-1]);
			Point2& p1(sPoints[i]);
			vSegments.push_back(Segment2(p0, p1));
		}
	}
	ConstraintGraph2* pCG = dt.createConstraint(vSegments, CIS_CONSTRAINED_DELAUNAY);

	// apply the segment constraint
	dt.applyConstraintsAndZones();
	if (path == "default")
	{
		string _name = name + ".ps";
		dt.show(_name);
	}
	else
	{
		string _name = path + name + ".ps";
		dt.show(_name);
	}

	// generate the seed point
	// we use mean value point here
	double x(0.0), y(0.0);
	for (int i = 0; i < vPoints.size(); ++i)
	{
		x += vPoints[i].x();
		y += vPoints[i].y();
	}
	x /= vPoints.size();
	y /= vPoints.size();

	Point2 seedPoint(x, y, 0.0f);
	vector<ConstraintGraph2*> vCG;
	vCG.push_back(pCG);
	Zone2* pGrowZone = dt.createZone(vCG, ZL_GROW, seedPoint);
	Zone2* pBoundedZone(pGrowZone->convertToBoundedZone());

	if (max_length == 0 && min_length == 0)
	{
		// Automatic calculation of the longest and shortest side lengths
		double mean_length = 0;
		for (int i = 0; i < vPoints.size(); ++i)
		{
			mean_length += sqrt(sqDistance2D(vPoints[i], vPoints[(i + 1) % vPoints.size()]));
		}
		mean_length /= vPoints.size();

		// 最后一个参数表示约束边是否可以被分割，通常设定为true
		// 但是在本程序中应该设定为false
		// 因为约束边在三角化之前是严格定义的，不可进行修改，
		// 因为一旦修改则会导致边界上出现多余的点则造成仿真出现错误
		dt.refine(pBoundedZone, 27, mean_length*3.0f, mean_length*1.5f, false);
	}
	else
	{
		dt.refine(pBoundedZone, 27, max_length, min_length, false);
	}

	// output mesh as .obj
	if (path == "default")
	{
		string _name = name + ".obj";
		dt.writeObj(_name, pBoundedZone);
	}
	else
	{
		string _name = path + name + ".obj";
		dt.writeObj(_name,pBoundedZone);
	} 
	return true;
}

//*************************** Trigonometric ***********************************

/* \brief Calculate cosα
*
* Calculate cosine value of two vectors
*
* @param vec_1 is first vector
* @param vec_2 is second vector
*
*	This method is to Calculate cosine value of two vectors
*/
double Cos(glm::vec3 vec_1, glm::vec3 vec_2)
{
	return (glm::dot(glm::normalize(vec_1), glm::normalize(vec_2)));
}

/* \brief Calculate sinα
*
* Calculate sine value of two vectors
*
* @param vec_1 is first vector
* @param vec_2 is second vector
*
*	This method is to Calculate sine value of two vectors
*/
double Sin(glm::vec3 vec_1, glm::vec3 vec_2)
{
	return (sqrt(1.0f - pow(Cos(vec_1,vec_2),2)));
}

/* \brief Calculate tanα
*
* Calculate tangent value of two vectors
*
* @param vec_1 is first vector
* @param vec_2 is second vector
*
*	This method is to Calculate tangent value of two vectors
*   Cosα could not be zero
*/
double Tan(glm::vec3 vec_1, glm::vec3 vec_2)
{
	double cos = Cos(vec_1, vec_2);
	double sin = Sin(vec_1, vec_2);

	assert(cos != 0);

	return (sin / cos);
}

//*************************** Geometry ***********************************

/* \brief Point to Line distance
*
* Calculate the distance of point and line in 3D space
*
* @param p is the point
* @param l1,l2 are points on the line
*
* This method is to Calculate the distance of point and line in 3D space
*
*/
double pointToLineDistance(glm::vec3 p, glm::vec3 l1, glm::vec3 l2)
{
	assert(p != l1 && p != l2);

	glm::vec3 line_vec = l2 - l1;
	glm::vec3 point_vec = p - l1;

	double d = glm::dot(line_vec, point_vec) / glm::length(line_vec);

	return (sqrt(pow(glm::length(point_vec), 2) - pow(d, 2)));
}

/* \brief Point to Plane distance
*
* Calculate the distance of point and plane in 3D space
*
* @param p is the point
* @param p1,p2,p3 are points on the plane
*
* This method is to Calculate the distance of point and plane in 3D space
*
*/
double pointToPlaneDistance(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 project_point = pointToPlaneProjection(p, p1, p2, p3);

	return (glm::distance(p, project_point));
}

/* \brief Point to Plane Project
*
* Calculate the project point of point to plane in 3D space
*
* @param p is the point
* @param p1,p2,p3 are points on the plane
*
* This method is to Calculate the project point of point to plane in 3D space
*
*/
glm::vec3 pointToPlaneProjection(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 plane_vec1 = p2 - p1;
	glm::vec3 plane_vec2 = p3 - p1;
	glm::vec3 point_vec  = p - p1;
	glm::vec3 unit_normal_vec = glm::normalize(glm::cross(plane_vec1, plane_vec2));

	float dot = glm::dot(point_vec, unit_normal_vec);
	assert(dot != 0);

	return (p - dot * unit_normal_vec);
}

/* \brief Point to Line Project
*
* Calculate the project point of point to line in 3D space
*
* @param p is the point
* @param l1,l2 are points on the plane
*
* This method is to Calculate the project point of point to line in 3D space
*
*/
glm::vec3 pointToLineProjection(glm::vec3 p, glm::vec3 l1, glm::vec3 l2)
{
	glm::vec3 line_vec = l2 - l1;
	glm::vec3 point_vec = p - l1;
	glm::vec3 unit_line_vec = glm::normalize(line_vec);

	float dot = glm::dot(point_vec, unit_line_vec);

	return (l1 + dot * unit_line_vec);
}

/* \brief Line & Line interseciton point
*
* Calculate the intersection point of line and line in 3D space
*
* @param l1_1,l1_2 are points of line_1
* @param l2_1,l2_2 are points of line_2
*
* This method is to Calculate the intersection point of line and line in 3D space
*
*/
glm::vec3 lineToLineIntersection(glm::vec3 l1_1, glm::vec3 l1_2, glm::vec3 l2_1, glm::vec3 l2_2)
{
	glm::vec3 line1 = l1_2 - l1_1;
	glm::vec3 line2 = l2_2 - l2_1;
	glm::vec3 norm1 = glm::normalize(line1);
	glm::vec3 norm2 = glm::normalize(line2);

	// whether two lines are collinear
	assert(glm::cross(line1.line2) != 0);

	// whether two lines are coplanar
	assert(glm::dot(glm::cross(l2_1 - l1_1, line1), line2) == 0);

	float t1 = glm::length(glm::cross(l1_1 - l2_1, norm2)) / glm::length(glm::cross(norm2, norm1));
	float t2 = glm::length(glm::cross(l2_1 - l1_1, norm1)) / glm::length(glm::cross(norm1, norm2));

	assert((l1_1 + t1 * norm1) == (l2_1 + t2 * norm2));

	return l1_1 + t1 * norm1;
}

/* \brief Plane & Line interseciton point
*
* Calculate the intersection point of line and plane in 3D space
*
* @param p1,p2,p3 are points in the plane
* @param l1,l2 are points on the line
*
* This method is to Calculate the intersection point of line and plane in 3D space
*  
*/
glm::vec3 lineToPlaneIntersection(glm::vec3 l1, glm::vec3 l2, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{
	glm::vec3 line_vec = l2 - l1;
	glm::vec3 plane_norm_vec = glm::normalize(glm::cross(p2 - p1, p3 - p1));

	assert(glm::dot(line_vec, plane_norm_vec) != 0);

	glm::vec3 point_vec = p1 - l1;

	float t = (glm::dot(point_vec, plane_norm_vec)) / (glm::dot(plane_norm_vec, line_vec));

	return l1 + t * line_vec;
}

/* \brief Point on Line
*
* Judge whether point on line
*
* @param p is point
* @param l1,l2 are points on the line
*
* This method is to Judge whether point on line
*
*/
bool isPointOnLine(glm::vec3 p, glm::vec3 l1, glm::vec3 l2)
{
	glm::vec3 line_vec = l2 - l1;
	glm::vec3 point_vec = p - l1;

//	return (glm::cross(line_vec, point_vec) == 0.0f);
}

/* \brief Point on Segment
*
* Judge whether point on segment
*
* @param p is point
* @param l1,l2 are points on the segment
*
* This method is to Judge whether point on segment
*
*/
bool isPointOnSegment(glm::vec3 p, glm::vec3 s1, glm::vec3 s2)
{

}

/* \brief Point in Plane
*
* Judge whether point in plane
*
* @param p is point
* @param p1,p2,p3 are points on the plane
*
* This method is to  Judge whether point in plane
*
*/
bool isPointInPlane(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{

}

/* \brief Point in Triangle
*
* Judge whether point in Triangle
*
* @param p is point
* @param p1,p2,p3 are points on the Triangle
*
* This method is to  Judge whether point in Triangle
*
*/
bool isPointInTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3)
{

}