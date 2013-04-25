/** fluid.h
 ** Brief: This is the header file of the Fluid class.
 ** Project: large-scale fluids
 ** Date: 04/16/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef FLUIDGPU_H
#define FLUIDGPU_H

#include <QVector>
#include "types.h"
#include "terrain.h"
#include "glwidget.h"

/**
 *  The default parameter for fluids are now define here
 **/
#include "fluid_global.h"
//#define FLUID_DEBUG

class FluidGPU
{

public:

    FluidGPU();
    FluidGPU(const int gridSize, const float domainSize );
    // Initialize from terrain, we should use this
    FluidGPU( Terrain* t );
    ~FluidGPU();
    void draw() const; //the name says it all - draw some fluid

    /**
     * Update the simulation at each time step
    **/
    void update(const float dt);

    /**
     * @brief Increment the height of the rectangular region around (posX, posZ) by incHeight
     * @param posX The x position
     * @param posZ The z position
     * @param radius The radius
     * @param incHeight The amount of height added
     */
    void incrementH(const int posX, const int posZ, const int radius, const float incHeight );

    /**
     * @brief addDrop Add the drop to specified region
     * @param posX The x position
     * @param posZ The z position
     */
    void addDrop( const int posX, const int posZ );
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

    /**
     * @brief backupHeight backUp the height information into fluid
     * @param t The terrain
     */
    void backupHeight( Terrain* t );

 private:

     friend void GLWidget::intersectFluid(const int x, const int y);
    /**
     * @brief init Initialize the variables
     * @param gridSize The length of the grid
     */
    void init(const int gridSize = GRID_SIZE, const float domainSize = DOMAIN_SIZE);

    /**
     * @brief Advect the array
     */
    void advect(  FieldType type, float* vec );

    /**
     * @brief updateVelocities Update the velocities
     */
    void updateVelocities();

    /**
     * @brief updateDepth Update the depth field
     */
    void updateDepth();

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
     * @brief Write the depth field or velocity to image
     */
    void saveToImage( FieldType type );

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

    /**
     * @brief initDepthField Initialize the depth field
     */
    void initDepthField( );

    /**
     * @brief build the triangle List
     */
    void buildTriangleList();

    /**
     * @brief apply when dampening the waves to simulate open water scenes
     * implements section 2.1.4 and appendix 3.1
     * to use, define DAMPEN_WAVES above
     */
    void dampenWaves();

    /**
     * @brief computes the resting height of the fluid
     * the resting height is computed as the average across the depth field
     * the result of this function is used in dampenWaves()
     * @return h_rest
     */
    float computeHRest();

    /**
     * @brief applies stability enhancements on h, u, w, as described in 2.1.5
     */
    void clampFields();

    /**
     * @brief getIndex1D Return the corresponding 1D index based on which type
     * @param i The row number
     * @param j The col number
     * @param type The type of the field
     * @return The 1D index
     */
    inline int getIndex1D( int i, int j, FieldType type) const;

private:

// Variables
    int m_gridSize;
    int m_uWidth; // width for velocityU (m_gridSize+1)

    float m_domainSize;
    float m_dx;
    float m_dt;
    float m_dxInv;
    float m_timeElapsed; // For debug
    int m_updateCount; // For debug

    bool m_renderNormals; // For debug
    Colorf m_color; // The color of the water

    /**
     * Design for gpu fluid. 2D vector is not suitable for cuda
     **/
    Vector3* m_normalField;
    float* m_tempBuffer;
    float* m_depthField;
    float* m_velocityU;
    float* m_velocityW;
    float* m_terrainHeightField;
    float* m_heightField;
    float* m_sigmaField;
    float* m_gammaField;
    float* m_phiField;
    float* m_psiField;

     QVector<Tri> m_triangles;
};

#endif // FLUIDGPU_H
