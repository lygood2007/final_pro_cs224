/** utils.cpp
 ** Brief: This is the source file of those useful functions.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "utils.h"
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
