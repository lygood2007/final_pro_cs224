/** random_terrain.h
 ** Brief: The header file of RandomTerrain Class
 ** Project: large-scale fluids
 ** Date: 04/15/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef RANDOMTERRAIN_H
#define RANDOMTERRAIN_H

#include "terrain.h"

#define DEFAULT_DECAY 2
#define DEFAULT_ROUGHNESS 2


class RandomTerrain : public Terrain
{
public:

    RandomTerrain(const int decay = DEFAULT_DECAY, const float roughness = DEFAULT_ROUGHNESS);
    virtual ~RandomTerrain();

private:


protected:

       /**
        * Computes the amount to perturb the height of the vertex currently being processed.
        * Feel free to modify this.
        *
        * @param depth The current recursion depth
        */
       double getPerturb(const int cur_depth) const;

       /**
       * Subdivides a square by finding the vertices at its corners, the midpoints of each side, and
       * the center (as the algorithm describes). Then recurs on each of the four sub-squares created.
       *
       * @param topLeft The grid coordinate of the top-left corner of the square to subdivide
       * @param bottomRight The grid coordinate of the bottom-right corner of the square to subdivide
       * @param depth The current recursion depth, decreasing as this function recurses deeper. The
       *              function should stop recurring when this value reaches zero.
       */
      void subdivideSquare(const GridIndex tlg, const GridIndex brg, GLint curDepth);

      /**
       * Sets default values for the four corners of the terrain grid and calls subdivideSquare()
       * to begin the terrain generation process. You do not need to modify this function.
       */
      virtual void populateTerrain();

protected:
      int m_decay; // Controls how much heights can vary per recursion depth level. Higher values generate smoother terrain.
      float m_roughness; // Controls the height of the terrain. Large values generate higher peaks and lower valleys.
};

#endif // RANDOMTERRAIN_H
