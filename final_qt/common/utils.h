/** utils.h
 ** Brief: This is the header file of those useful functions.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/
#ifndef UTILS_H
#define UTILS_H

#include <QVector>

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

#endif // UTILS_H
