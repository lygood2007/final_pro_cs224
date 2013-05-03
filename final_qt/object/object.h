#ifndef SHAPE_H
#define SHAPE_H

#include <qgl.h>
#include "GL/glut.h"
#include "types.h"
#include "vector.h"
#include "fluid_global.h"
#include <QList>
#include "CS123Algebra.h"
#include "object_defs.h"


#define JITTER_ORIGIN

#define MAX_INIT_U 20
#define MAX_INIT_W 20
#define MAG_U 2.2
#define MAG_W 2.2
#define MIN_DENSITY 300

class FluidGPU;

class Object
{
public:

    Object();
    Object( FluidGPU* fluid, Vector3 position, float dx, GLuint texID, float density, const Colorf color = Colorf(1.f,1.f,1.f,1.f) );
    virtual ~Object();

    void setFluid( FluidGPU* fluid );
    /**
     * @brief draw Draw the box
     */
    void draw();

    /**
     * @brief setColor set the object's color
     * @param color
     */
    inline void setColor( const Colorf color ) {
        m_color =  color;
    }

    /**
     * Enable drawing normals
     */
    inline void enableNormal(){m_renderNormals = true;}
    /**
     * Disable drawing normals
     */
    inline void disableNormal(){m_renderNormals = false;}
    /**
     * Get m_renderNormal
     */
    inline bool isRenderingNormal() const {return m_renderNormals;}

    /**
     * @brief update update the box's position
     */
    void update( float dt, FluidGPU* heightField );

    inline float getDensity() const { return m_density; }

    inline void setDensity( float density ){ if( density <= MIN_DENSITY || density >= 2000)return; m_density = density; computeMass(); }

    /**
     * @brief Compute mass, tesselation,etc
     */
    void initPhysics();

    /**
     * @brief setW Set the w parameter
     * @param w
     */

    inline void setW( float w ){ m_w = w; }

    /**
     * @brief setCoeffDrag Set the drag coefficient
     * @param dragCoef
     */
    inline void setCoeffDrag( float dragCoef ){ m_dragCoeff = dragCoef;c1 = -WATER_DENSITY*m_dragCoeff*0.5; }

    /**
     * @brief setCoeffLift Set the lift coefficient
     * @param liftCoef
     */
    inline void setCoeffLift( float liftCoef ){ m_liftCoeff = liftCoef; c2 = -WATER_DENSITY*m_liftCoeff*0.5;}

private:

    /**
     * @brief initFluidInfo Initialize the data from fluid
     */
    void initFluidInfo( FluidGPU* fluid );
    /**
     * @brief drawNormal Draw the normals of this box
     */
    void drawNormal() const;
    /**
         *  Update the triangles' postion and norms (average position and normal)
         *  The are used for buoyance computation
         */
    void updatePosAndNorm();
    /**
     * @brief updatePosWithBoundaryCheck update the next position, avoiding moving out of boundary
     */
    void updatePosWithBoundaryCheck( float dt );

    /**
     * @brief bilinearInterpReal Get
     * @brief This function will do bilinear intepolation based on x,z value of the vector
     * @param fluid The fluid
     * @param pos The position
     * @param fv the fluid velocity
     * @param the interpolated height
     */
     void getInterpVelocityAndHeight( Vector3 pos,
                                       Vector3& fv, float&h );
     /**
      * @brief computeOrigin compute the origin position of x and z
      */
     void computeOrigin();

private:

    GLuint m_texID;
    bool m_renderNormals;
    float m_domainSize;
    float m_origX;
    float m_origZ;

    float m_dragCoeff;
    float m_liftCoeff;

    float m_angle[6];
    Matrix4x4 m_rotMat[3];

    Colorf m_color;

    bool m_upwards;
    bool m_decay;
    float m_w;

    Vector3 m_lastBuoAcc;

    // Only for speeding up
    float c1;
    float c2;
    float c3;

    float* m_curH;
    float* m_curU;
    float* m_curW;

    int m_gridSize;
    int m_uwidth;
    int m_wheight;

protected:

    /**
     * @brief compute the tessellation
     */
    virtual void computeTessel() = 0; // Compute the best tesselation parameter (slices)
    /**
     * @brief buildTriangleList build the triangle list
     */
    virtual void buildTriangleList() = 0;

    /**
     * @brief computeMass Compute the mass of the object
     */
    virtual void computeMass() = 0;

    /**
     * @brief generate the rotation randomly
     **/
     void generateRotation();

     /**
      * Get the rotated position using m_angle which defines the angle of ration along x,y,z     *
      */
     Vector3 getRotVec( Vector3 vec );

protected:

      float m_dx;
      float m_dxInv;
      //float m_fluiddx; // backup the dx of fluid's grid

     Vector3 m_velocity;
     Vector3 m_position;

    QList<Tri> m_tris; // Triangles;

    int m_tessell[3];

    float m_mass;
    float m_massInv;

    float m_density;

};

#endif // SHAPE_H
