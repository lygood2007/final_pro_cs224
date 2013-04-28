/** fluid.h
 ** Brief: This is the header file of the Fluid class.
 ** Project: large-scale fluids
 ** Date: 04/16/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef FLUIDCPU_H
#define FLUIDCPU_H

#include <QVector>
#include "types.h"
#include "terrain.h"
#include "glwidget.h"
#include "particle.h"

/**
 *  The default parameter for fluids are now define here
 **/
#include "fluid_global.h"
//#define FLUID_DEBUG

/*
#define GRID_SIZE TERRAIN_BOUND*2
#define DOMAIN_SIZE TERRAIN_BOUND
#define GRAVITY -10
#define SAVE_NAME_HEIGHT "height"
#define SAVE_NAME_VELOCITY "velocity"
//#define SAVE_IMAGE

#define DAMPEN_WAVES
#define LAMBDA_DECAY 0.9
#define LAMBDA_UPDATE 0.1
#define DAMPENING_REGION 5
#define QUADRATIC_A 10.0f
#define QUADRATIC_B 0.0f
#define QUADRATIC_C 0.0f
#define INIT_SIGMA_GAMMA 0.0f
#define INIT_PHI_PSI 0.0f

#define CLAMP_ALPHA 0.5f*/

#define C_DEPOSIT 1
#define SPLASH_PARTICLE_RADIUS 0.1

#define ALPHA_MIN_SPLASH 0.45
#define V_MIN_SPLASH 4
#define L_MIN_SPLASH -4
#define BREAKING_WAVE_NUM_SPLASH_PARTICLES 50

//extern float bilinearInterp( QVector<QVector<float > > &vec, const float x, const float z );
//extern float randomFloatGenerator( float min = 0.f, float max = 1.f);

class FluidCPU
{

public:

    FluidCPU();
    FluidCPU(const int gridSize, const float domainSize );
    // Initialize from terrain, we should use this
    FluidCPU( Terrain* t );
    ~FluidCPU();
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
     * @brief updates the positions and velocities of the current particles
     */
    void updateParticles();

    /**
     * @brief checks to see which particles should be removed from the list
     */
    void removeParticles();

    /**
     * @brief checks if this particle has hit the fluid, if so, returns height to fluid
     */
    bool fluidParticleInteraction(Particle *particle);

    /**
     * @brief renders the particles
     */
    void drawParticles() const;

    /**
     * @brief generate multiple splash particles here
     */
    void generateSplashParticles(int i, int j, int numParticles);

    /**
     * @brief generate a splash particle at (i, j) on the grid
     */
    void generateSplashParticle(int i, int j, Vector3 position);

    /**
     * @brief check grid for breaking waves, generate splash particles if so
     */
    void checkForBreakingWaves();

    /**
     * @brief computes the values of the three conditions for breaking waves
     */
    double computeBreakingWaveCondition1(int i, int j);
    double computeBreakingWaveCondition2(int i, int j);
    double computeBreakingWaveCondition3(int i, int j);

private:
// Variables
    int m_gridSize;

    float m_domainSize;

    float m_dx;
    float m_dt;
    float m_dxInv;
    float m_timeElapsed; // For debug
    int m_updateCount; // For debug

    bool m_renderNormals; // For debug
    Colorf m_color; // The color of the water

    float** m_tempBuffer;

    /**
     * Design for cpu fluid. 2D vector is not suitable for cuda
     **/

    QVector<QVector<float> > m_depthField; // Stores the height
    QVector<QVector<Vector3> > m_normalField; // Stroes the normal, a better way is to declare this as a 2D array,

    QVector<QVector<float> > m_velocityU;
    QVector<QVector<float> > m_velocityW;
    QVector<QVector<float> > m_terrainHeightField;
    QVector<Tri> m_triangles;
    QVector<QVector<float> > m_sigmaField; // stores sigma values for wave dampening
    QVector<QVector<float> > m_gammaField; // stores gamma values for wave dampening
    QVector<QVector<float> > m_phiField; // stores phi values for wave dampening
    QVector<QVector<float> > m_psiField; // stores psi values for wave dampening

    QVector<Particle*> m_particles; // stores the particles
    QVector<QVector<float> > m_depthFieldPrev; // stores the previous values of the heights
};

#endif // FLUIDCPU_H
