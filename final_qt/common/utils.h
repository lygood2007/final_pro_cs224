/** utils.h
 ** Brief: This is the header file of those useful functions.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#ifndef UTILS_H
#define UTILS_H

#include <QVector>
#include "vector.h"
/**
  * @brief bilinearInterp 2D dimension bilinear interpolation
  * @param vec The 2D array
  * @param x The x coordinate
  * @param z The z coordiante
  * @return The result
  */
 float bilinearInterp( QVector<QVector<float > > &vec, float x, float z );

 /**
  * @brief randomFloatGenerator Generate a float number between min and max
  * @param min The lower bound
  * @param max The upper bound
  * @return The random float number we get
  */
 float randomFloatGenerator( float min = 0.f, float max = 1.f);

 /**
  * @brief doIntersectTriangles Check if the ray intersects the triangle
  * @param eyePos The start point of the ray
  * @param d The direction of ray
  * @param v0 The vertex 0 of triangle
  * @param v1 The vertex 1 of triangle
  * @param v2 The vertex 2 of triangle
  * @return True if the ray intersects
  */
 bool doIntersectTriangles( const Vector3& eyePos, const Vector3& d, const Vector3& v0, const Vector3& v1, const Vector3& v2 );

#endif // UTILS_H
