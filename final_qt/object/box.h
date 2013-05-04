/** box.h
 ** Brief: This is the header file of the class Box
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef BOX_H
#define BOX_H

#include "object.h"

#define BOX_DENSITY 800

// A box with specified length, width, height

class FluidGPU;

class Box : public Object
{
public:
    Box();
    Box(  FluidGPU* fluid,Vector3 position, float length,
          float height, float width, float dx, GLuint texID, const Colorf color = Colorf(1.f,1.f,1.f,1.f));
    virtual ~Box();

protected:
    /**
     * @brief compute the tessellation
     */
    virtual void computeTessel(); // Compute the best tesselation parameter (slices)
    /**
     * @brief buildTriangleList build the triangle list
     */
    virtual void buildTriangleList();

    /**
     * @brief computeMass Compute the mass of the object
     */
    virtual void computeMass();

private:

    float m_length;
    float m_height;
    float m_width;
};

#endif
