/** sphere.h
 ** Brief: This is the header file of the class Sphere
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef SPHERE_H
#define SPHERE_H

#include "object.h"

// A sphere with specified length, width, height

#define SPHERE_DENSITY 600

class FluidGPU;

class Sphere : public Object
{
public:
    Sphere();
    Sphere(  FluidGPU* fluid,Vector3 position, float radius, float dx, GLuint texID,
             const Colorf color = Colorf(1.f,1.f,1.f,1.f));
    virtual ~Sphere();

protected:
    /**
     * @brief Compute the tessellation
     */
    virtual void computeTessel(); // Compute the best tesselation parameter (slices)
    /**
     * @brief buildTriangleList Build the triangle list
     */
    virtual void buildTriangleList();

    /**
     * @brief computeMass Compute the mass of the object
     */
    virtual void computeMass();

private:

    float m_radius;
};
#endif // SPHERE_H
