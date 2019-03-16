# Common-CG-Algrithms
  This is a very simple c++ head file which contains some common used 
  algirthms in Compute Graphics.
## Development Environment
  * [Fade2d](https://www.geom.at/fade2d/html/)
  * [GLM](https://github.com/g-truc/glm)
## Algrithms
  * double Cos(glm::vec3 vec_1, glm::vec3 vec_2) : 
      * Calculate cosine value of two vectors
  * double S(glm::vec3 vec_1, glm::vec3 vec_2) : 
      * Calculate sine value of two vectors
  * double Tan(glm::vec3 vec_1, glm::vec3 vec_2) : 
      * Calculate tangent value of two vectors
  * glm::vec2 getVerticalUnitVec(glm::vec2 vec, bool is_left) :
      * Get the Vertical Unit vector of a given vector
  * double pointToLineDistance(glm::vec3 p, glm::vec3 l1, glm::vec3 l2) : 
      * Calculate the distance of point and line in 3D space
  * double pointToPlaneDistance(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) : 
      * Calculate the distance of point and plane in 3D space
  * glm::vec3 pointToPlaneProjection(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) : 
      * Calculate the project point of point to plane in 3D space
  * glm::vec3 pointToLineProjection(glm::vec3 p, glm::vec3 l1, glm::vec3 l2) : 
      * Calculate the project point of point to line in 3D space
  * glm::vec3 lineToLineIntersection(glm::vec3 l1_1, glm::vec3 l1_2, glm::vec3 l2_1, glm::vec3 l2_2) : 
      * Calculate the intersection point of line and line in 3D space
  * glm::vec3 lineToPlaneIntersection(glm::vec3 l1, glm::vec3 l2, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) : 
      * Calculate the intersection point of line and plane in 3D space
  * glm::vec2 segToSegIntersection2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2);
      * Calculate the intersection point of Segment and Segment in 2D space
  * bool isPointOnLine(glm::vec3 p, glm::vec3 l1, glm::vec3 l2) : 
      * Judge whether point on line
  * bool isPointOnSegment(glm::vec3 p, glm::vec3 s1, glm::vec3 s2) : 
      * Judge whether point on segment
  * bool isPointInPlane(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) : 
      * Judge whether point in plane
  * bool isPointInTriangle(glm::vec3 p, glm::vec3 p1, glm::vec3 p2, glm::vec3 p3) : 
      * Judge whether point in Triangle
  * bool isPointInPolygon2D(glm::vec2 p, vector<glm::vec2> polygon) :
      * Judge whether point in or out polygon
  * bool isSegmentInPolygon2D(glm::vec2 s1, glm::vec2 s2, vector<glm::vec2> polygon) :
      * Judge whether segment in or out polygon
  * bool isSegmentIntersect2D(glm::vec2 s1_1, glm::vec2 s1_2, glm::vec2 s2_1, glm::vec2 s2_2) :
      * Judge whether two segments intersect
      
  * To be continue .....
