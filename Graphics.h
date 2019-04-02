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
#include <map>
#include <stack>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <assert.h>
#include <glm\glm.hpp>
#include <Fade_2D.h>

using namespace std;
using namespace GEOM_FADE25D;

/************************ Struct Data *********************************/
struct Triangle
{
	int v0;
	int v1;
	int v2;
};

/********************** Declare All Functions Here ***********************/
void readOBJ(vector<glm::vec2>& point_vec, vector<Triangle>& triangle_vec, string path);
void writeOBJ(const vector<glm::vec2>& point_vec, const vector<Triangle>& triangle_vec, string path);

void Out(glm::vec2 vec_2);
void Out(glm::vec3 vec_3);
void Out(glm::vec4 vec_4);

double Cos(glm::vec3 vec_1, glm::vec3 vec_2);
double Cos2d(glm::vec2 vec_1, glm::vec2 vec_2);
double Sin(glm::vec3 vec_1, glm::vec3 vec_2);
double Sin2d(glm::vec2 vec_1, glm::vec2 vec_2);
double Tan(glm::vec3 vec_1, glm::vec3 vec_2);
double Tan2d(glm::vec2 vec_1, glm::vec2 vec_2);

vector<glm::vec2> oneBezierInterpolation(glm::vec2 p1, glm::vec2 p2, int num);
vector<glm::vec2> twoBezierInterpolation(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, int num);
vector<glm::vec2> threeBezierInterpolation(glm::vec2 p1, glm::vec2 p2, glm::vec2 p3, glm::vec2 p4, int num);

glm::vec2 getVerticalUnitVec(glm::vec2 vec, bool is_left);

double pointToLineDistance(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);
double pointToPlaneDistance(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);

glm::vec3 pointToPlaneProjection(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
glm::vec3 pointToLineProjection(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);

glm::vec3 lineToLineIntersection(glm::vec3 l1_1, glm::vec3 l1_2, glm::vec3 l2_1, glm::vec3 l2_2);
glm::vec3 lineToPlaneIntersection(glm::vec3 l1, glm::vec3 l2, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
glm::vec2 segToSegIntersection2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2);

bool isPointOnLine(glm::vec3 p, glm::vec3 l1, glm::vec3 l2);
bool isPointOnSegment(glm::vec3 p, glm::vec3 s1, glm::vec3 s2);
bool isPointInPlane(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
bool isPointInTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3);
bool isPointByTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3, float& u, float& v);
bool isPointInPolygon2D(glm::vec2 p, vector<glm::vec2> polygon);
bool isSegmentInPolygon2D(glm::vec2 s1, glm::vec2 s2, vector<glm::vec2> polygon);
bool isSegmentIntersect2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2);

vector<glm::vec2> getConvexHull(vector<glm::vec2> points);
vector<glm::vec2> getCoutourOfNonConvex2dMesh(string path);
vector<int> getCoutourIndexOfNonConvex2dMesh(string path);
vector<glm::vec2> getCoutourOfNonConvex2dMesh(const vector<glm::vec2>& point_vec, const vector<Triangle>& triangle_vec);
vector<int> getCoutourIndexOfNonConvex2dMesh(const vector<glm::vec2>& point_vec, const vector<Triangle>& triangle_vec);

/********************** Declare All Functions Here ***********************/

