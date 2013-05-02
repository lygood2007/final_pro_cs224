/** utils.cpp
 ** Brief: This is the source file of those useful functions.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "utils.h"
#include "CS123Common.h"
#include "debug_marco.h"
/**
  * @brief bilinearInterp 2D dimension bilinear interpolation
  * @param vec The 2D array
  * @param x The x coordinate
  * @param z The z coordiante
  * @return The result
  */
 float bilinearInterp( QVector<QVector<float > > &vec, float x, float z )
 {
     /**
      * Clamp firstly
      */

     /**
      *  Vec must be square!
      */
     assert( vec.size() > 0 && vec[0].size() > 0 );
     const int width = vec[0].size();
     const int height = vec.size();

     if( x < 0 )
         x = 0.f;
     if( z < 0 )
         z = 0.f;
     if( x > width - 1 )
         x = width - 1;
     if( z > height - 1 )
         z = height -1;

     const int X = (int)x;
     const int Y = (int)z;
     const float s1 = x - X;
     const float s0 = 1.f - s1;
     const float t1 = z - Y;
     const float t0 = 1.f-t1;
     float e1, e2, e3,e4;
     e1 = e2 = e3 = e4 = 0;
     e1 = vec[Y][X];
     if( Y+1 <= height- 1 )
         e2 = vec[Y+1][X];
     if( X +1 <= width -1 )
         e3 = vec[Y][X+1];
     if( Y+1 <= height - 1 && X + 1 <= width - 1)
         e4 = vec[Y+1][X+1];

     float result = s0*(t0*e1 + t1*e2 )+
             s1*(t0*e3  + t1*e4 );

     return  result;
 }

 /**
  * @brief bilinearInterpReal 2D dimension bilinear interpolation based on real vector, not the index.
  * @brief This function will do bilinear intepolation based on x,z value of the vector
  * @param vec The 2D array
  * @param pos The position
  * @param dx The delta x
  * @param width The width of the grid
  * @param height The height of the grid
  * @return The result y value
  */
 float bilinearInterpReal( const float* vec, Vector3 pos, float origX, float origZ, float dx, const int width, const int height )
 {
     float x = pos.x - origX;
     float z = pos.z - origZ;
     if( x < 0 )
         x = 0.f;
     if( z < 0 )
         z = 0.f;
     if( x > (width - 1)*dx )
         x = (width - 1)*dx;
     if( z > (height - 1)*dx )
         z = (height -1)*dx;

     // Get the index
     const int X = (int)(x/dx);
     const int Y = (int)(z/dx);
     const float s1 = x - X*dx;
     const float s0 = 1.f - s1;
     const float t1 = z - Y*dx;
     const float t0 = 1.f-t1;
     float e1, e2, e3,e4;
     e1 = e2 = e3 = e4 = 0;
     e1 = vec[Y*width+X];
     if( Y+1 <= height- 1 )
     {
         e2 = vec[(Y+1)*width + X];
     }
     if( X +1 <= width -1 )
     {
         e3 = vec[Y*width + X+1];
     }
     if( Y+1 <= height - 1 && X + 1 <= width - 1)
     {
         e4 = vec[(Y+1)*width + X+1];
     }

     float result = s0*(t0*e1 + t1*e2 )+
             s1*(t0*e3  + t1*e4 );

     return  result;
 }

 /**
  * @brief randomFloatGenerator Generate a float number between min and max
  * @param min The lower bound
  * @param max The upper bound
  * @return The random float number we get
  */
 float randomFloatGenerator( float min, float max )
 {
     float result = min + (rand() / (RAND_MAX + 0.01f)*(max - min));
     return result;
 }

 /**
  * @brief doIntersectTriangles Check if the ray intersects the triangle
  * @param eyePos The start point of the ray
  * @param d The direction of ray
  * @param v0 The vertex 0 of triangle
  * @param v1 The vertex 1 of triangle
  * @param v2 The vertex 2 of triangle
  * @return True if the ray intersects
  */
 bool doIntersectTriangles( const Vector3& eyePos, const Vector3& d, const Vector3& v0, const Vector3& v1, const Vector3& v2 )
 {
         float a,f,u,v;
         Vector3 e1 = v1 - v0;
         Vector3 e2 = v2 - v0;
         Vector3 h = d.cross(e2);

         a = e1.dot(h);

         if (a > -EPSILON && a < EPSILON)
             return false;
         f = 1/a;
         Vector3 s = eyePos - v0;
         u = f * s.dot(h) ;

         if (u < 0.0 || u > 1.0)
             return false;

         Vector3 q = s.cross(e1);
         v = f * d.dot(q);

         if (v < 0.0 || u + v > 1.0)
             return false;

         float t = f * e2.dot(q);

         if (t > EPSILON) // ray intersection
             return true;
 }

