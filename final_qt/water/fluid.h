/** fluid.h
 ** Brief: This is the header file of the Fluid class.
 ** Project: large-scale fluids
 ** Date: 04/16/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef FLUID_H
#define FLUID_H

#include <QVector>
#include "types.h"
#include "terrain.h"
//#define FLUID_DEBUG

#define GRID_SIZE TERRAIN_BOUND*2
#define DOMAIN_SIZE TERRAIN_BOUND
#define GRAVITY -10

extern float bilinearInterp( QVector<QVector<float > > &vec, const float x, const float z );
extern float randomFloatGenerator( float min = 0.f, float max = 1.f);

class Fluid
{

public:

    enum FieldType
    {
        HEIGHT = 0,
        VELOCITY_U,
        VELOCITY_W,
        VELOCITY
    };

    enum DrawMethod
    {
        DRAW_POINTS = 0,
        DRAW_MESH
    };

    Fluid();
    Fluid(const int gridSize, const float domainSize );
    ~Fluid();
    void draw() const; //the name says it all - draw some fluid

    /**
     * Update the simulation at each time step
    **/
    void update(const float dt);

    /**
     * @brief Add the drop to the water
     * @param posX The x position
     * @param posZ The y position
     */
    void addDrop(const int posX, const int posZ);

    /**
     * @brief Add drop to random positions
     *
     */
    void addRandomDrop( const float freq = 0.025 );

    /**
     * @brief setColor Set the fluid's color
     * @param color The parameter
     */
    void setColor(const Colorf color);

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

 private:

    /**
     * @brief init Initialize the variables
     * @param gridSize The length of the grid
     */
    void init(const int gridSize = GRID_SIZE, const float domainSize = DOMAIN_SIZE);

    /**
     * @brief Advect the array
     */
    void advect(  FieldType type, QVector<QVector<float> >& vec );

    /**
     * @brief updateVelocities Update the velocities
     */
    void updateVelocities();

    /**
     * @brief updateHeight Update the height field
     */
    void updateHeight();

    /**
     * @brief Apply boundary condition
     */
    void applyBoundary();

    /**
     * @brief Check boundary
     */
    void checkBoundary();

    /**
     * @brief Write the height field or velocity to image
     */
    void writeToImage( FieldType type );

    /**
     * @brief drawFluid Draw the fluid with different method
     * @param method The method for drawing, could be DRAW_POINTS or DRAW_MESH
     */
    void drawFluid( DrawMethod method ) const;
    /**
     * @brief Computet the normal of each points
     */
    void computeNormal();

    /**
     * @brief Draw the normals of the fluid points
     */
    void drawNormal() const;


    int m_gridSize;

    float m_domainSize;

    float m_dx;
    float m_dt;
    float m_dxInv;

    QVector<QVector<float> > m_heightField; // Stores the height
    QVector<QVector<Vector3> > m_normalField; // Stroes the normal, a better way is to declare this as a 2D array,
                                                  // To make it compatible with terrain, I declare as pointer
    QVector<QVector<float> > m_velocityU;
    QVector<QVector<float> > m_velocityW;

    float m_timeElapsed; // For debug
    int m_updateCount; // For debug

    bool m_renderNormals; // For debug
    Colorf m_color; // The color of the water

    float** m_tempBuffer;
};

#endif // FLUID_H
