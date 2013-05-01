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
#include "particle.h"
#include "particlesource.h"

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
     * @brief addDrop Add the drop to specified region
     * @param posX The x position
     * @param posZ The z position
     */
    void addDrop( const int posX, const int posZ );

    /**
     * @brief drop particles on a small area
     * @param posX The x position
     * @param posZ The z position
     */
    void addDroppingParticles( const int posX, const int posZ );

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

    /**
     * @brief getGridSize Return the gird size
     * @return The grid size
     */
    inline int getGridSize() const {return m_gridSize;}

    /**
     * @brief getFieldArray Get the pointer to the array with specified type
     * @param buffLength The buffer length returned
     * @return The buffer
     */
    float* getFieldArray( FieldType type, int& buffLength ) const;

 private:

     friend void GLWidget::intersectFluid(const int x, const int y, QMouseEvent *event);
    /**
     * @brief init Initialize the variables
     * @param gridSize The length of the grid
     */
    void init(const int gridSize = GRID_SIZE, const float domainSize = DOMAIN_SIZE);

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
     * @brief Draw the normals of the fluid points
     */
    void drawNormal() const;

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

    /**
     * @brief updates the positions and velocities of the current particles
     */
    void updateParticles();

    /**
     * @brief checks to see which particles should be removed from the list
     */
    void removeParticles();

    /**
     * @brief checks if this particle has hit the fluid, if so, returns height to fluid
     * DEPRACATED, USE FUNCTIONS BELOW
     */
    bool fluidParticleInteraction(Particle *particle);

    /**
     * @brief better particle/fluid interaction
     */
    Vector3 fluidParticleInteractionCheck(Particle *particle);
    void fluidParticleUpdate(QVector<Vector3> values);

    /**
     * @brief renders the particles
     */
    void drawParticles() const;

    /**
     * @brief create particle sources
     */
    void initParticleSources();

    /**
     * @brief generate new particles for this timestep from the sources
     */
    void updateParticleSources();

    void drawParticles2() const;

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
    Vector3* m_normalField; // Essential
    float* m_depthField; // Essential
    float* m_velocityU; // I preserve this for future use
    float* m_velocityW; // I preserve this for future use
    float* m_terrainHeightField; // Essential
    float* m_heightField; // Essential

    /**
     * dampening waves structures
     */
    float* m_sigmaField;
    float* m_gammaField;
    float* m_phiField;
    float* m_psiField;

    QVector<Tri> m_triangles;

    /**
     * first attempt at particles
     */
    QVector<Particle*> m_particles; // stores the particles
    QVector<ParticleSource*> m_particle_sources; // stores sources of particles (faucets)

    /**
     * second attempt at particles
     */
    Vector3 *m_particle_positions; // positions of particles, y = -1 is inactive
    Vector3 *m_particle_velocities; // velocities of particles
    Vector3 m_particle_acceleration; // acceleration of particles (same for all)
};

#endif // FLUIDGPU_H
