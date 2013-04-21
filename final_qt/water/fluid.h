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
#include "glwidget.h"
//#define FLUID_DEBUG

#define GRID_SIZE TERRAIN_BOUND*2
#define DOMAIN_SIZE TERRAIN_BOUND
#define GRAVITY -10
#define SAVE_NAME_HEIGHT "height"
#define SAVE_NAME_VELOCITY "velocity"
//#define SAVE_IMAGE

//#define DAMPEN_WAVES
#define LAMBDA_DECAY 0.9
#define LAMBDA_UPDATE 0.1
#define DAMPENING_REGION 10
#define QUADRATIC_A 10.0f
#define QUADRATIC_B 0.0f
#define QUADRATIC_C 0.0f
#define INIT_SIGMA_GAMMA 0.0f
#define INIT_PHI_PSI 0.0f

#define CLAMP_ALPHA 0.5f

//extern float bilinearInterp( QVector<QVector<float > > &vec, const float x, const float z );
//extern float randomFloatGenerator( float min = 0.f, float max = 1.f);


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
    // Initialize from terrain, we should use this
    Fluid( Terrain* t );
    ~Fluid();
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

private:
// Variables
    int m_gridSize;

    float m_domainSize;

    float m_dx;
    float m_dt;
    float m_dxInv;

    QVector<QVector<float> > m_depthField; // Stores the height
    QVector<QVector<Vector3> > m_normalField; // Stroes the normal, a better way is to declare this as a 2D array,
                                                  // To make it compatible with terrain, I declare as pointer
    QVector<QVector<float> > m_velocityU;
    QVector<QVector<float> > m_velocityW;
    QVector<QVector<float> > m_terrainHeightField;
    QVector<Tri> m_triangles;
    float m_timeElapsed; // For debug
    int m_updateCount; // For debug

    bool m_renderNormals; // For debug
    Colorf m_color; // The color of the water

    float** m_tempBuffer;

    QVector<QVector<float> > m_sigmaField; // stores sigma values for wave dampening
    QVector<QVector<float> > m_gammaField; // stores gamma values for wave dampening
    QVector<QVector<float> > m_phiField; // stores phi values for wave dampening
    QVector<QVector<float> > m_psiField; // stores psi values for wave dampening
};

#endif // FLUID_H
