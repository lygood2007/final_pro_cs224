/** glwidget.h
 ** Brief: This is the header file of the class containing the framework of the project.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef GL_WIDGET_H
#define GL_WIDGET_H

#include <qgl.h>
#include <QTime>
#include <QTimer>
#include <QHash>
#include <QString>
#include "utils.h"
#include "camera.h"
#include "fluid_global.h"
#include "resourceloader.h"
#include "object_defs.h"
// We fix the size
#define WIN_W 1000.0
#define WIN_H 700.0

// Flag for testing
#define DRAW_TERRAIN

#define USE_HEIGHTMAP

#define RENDER_FLUID
//#define USE_FBO
//#define USE_SKYBOX

// Colors to use when rendering
#define SEA_WATER 0.0f,0.42f,0.58f,0.9f
#define TIME_STEP 0.03 //0.03 //is max for gridsize 80
/**
    Uncomment this if you don't want to use CUDA to compute
**/
#define USE_GPU_FLUID


class QGLShaderProgram;
class QGLFramebufferObject;
class Terrain;
class FluidCPU;
class FluidGPU;
class Object;

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent);
    ~GLWidget();

    //Not the best idea but I'm placing these here in case I need them elsewhere
    bool m_useShaders, m_useFBO, m_useSimpleCube, m_useAxis,
        m_useSkybox, m_useParticles, m_useDampening,
        m_useParticleSources, m_useRectangularParticleSources;

private:
    // Private variables
    QTime m_time; // The time variable
    QTimer m_timer; // The timer variable
    OrbitCamera m_camera; // Camera
    Terrain* m_terrain;
    QList<Object*> m_objects;
    GLuint m_boxTexID;
    GLuint m_sphereTexID;

#ifdef USE_GPU_FLUID
    FluidGPU* m_fluid;
#else
    FluidCPU* m_fluid;
#endif

    QHash<QString, QGLShaderProgram *> m_shaderPrograms; // hash map of all shader programs
    QHash<QString, QGLFramebufferObject *> m_framebufferObjects; // hash map of all framebuffer objects
    GLuint m_skybox; // skybox call list ID
    GLuint m_cubeMap; // cubeMap texture ID
    QFont m_font; // font for rendering tex


    bool m_mouseLeftDown; // True if mouse left is down
    bool m_mouseRightDown; // True if mouse right is down

    bool m_drawFrame; // True if draw in wireframe mode
    bool m_animate;

    int m_prevTime;
    float m_prevFps, m_fps;
    float m_delta;

    float m_timeStep;

    GLuint m_waterNormalMap;


public:
    // Private functions
    /** Initliaze the variables in the class*/
    void init();
    /** Initialize the GL states, put all your GL initialization code here */
    void initializeGL();
    /** The render function, will be called every frame */
    void paintGL();

    /** Putting all our resource init's in one place*/
    void initializeResources();

    /** Copied from CS123 Bloom lab - these really aren't the best and could be improved*/
    void applyOrthogonalCamera(float width, float height);
    void applyPerspectiveCamera(float width, float height);

    /** Placing all visible geometry rendering in one method*/
    void renderScene();

    /** The splitting these out to improve rendering and shader passes*/
    void renderTexturedQuad(int width, int height);
    void renderSkybox();

    /**
     * @brief GLWidget::renderObjects render the objects
     */
    void renderObjects();

    /** Placeing all visible geometry rendering in one method*/
    void renderFluid();
    void renderParticles();
    /** So we can render the particles separately for shading purposes*/
    void renderSpray();
    void renderSplash();
    void renderFoam();

    /** Updates the current OpenGL state to avoid object distortion when the window is resized. */
    void resizeGL(int w, int h);

    /** Handles all the mouse press event */
    void mousePressEvent(QMouseEvent *event);
    /** Handles all the mouse move event */
    void mouseMoveEvent(QMouseEvent *event);
    /** Handles all the mouse release event */
    void mouseReleaseEvent(QMouseEvent *event);
    /** Handles all the mouse wheel evet */
    void wheelEvent(QWheelEvent *event);

    /** Handles all the key press event */
    void keyPressEvent(QKeyEvent *event);
    /** Handles all the key release event */
    void keyReleaseEvent(QKeyEvent *event);

    /**
     * Called whenever m_camera is modified or the canvas is resized. Sets the current OpenGL projection
     * and modelview matrices to match the values in m_camera.
     */
    void updateCamera();

    /**
     * Render the text. If you want to add more texts, put it here
     **/
    void paintText();
    /**
     * Update the time variables
     **/
    void timeUpdate();

    /**
     * @brief intersectFluid Check if the ray shooting from position (x, y)  intersects the fluid
     * @param x, The x position in screen space
     * @param y, The y position in screen space
     * @return Return if it is intersected
     */
    bool intersectFluid(  const int x, const int y, int& indexRow, int& indexCol, Vector3& pos );

    //copied from CS123 lab 09 - SH
    void loadCubeMap();
    void loadObjectTexMap();
    void createShaderPrograms();
    void createFramebufferObjects(int width, int height);
    void createBlurKernel(int radius, int width, int height, GLfloat* kernel, GLfloat* offsets);

    void renderBlur(int width, int height);

    /**
     * @brief addObject Drop a object from the air with object's type specified by type
     * @param x The x position
     * @param z The z position
     * @param type The object's type
     * @param Height The height
     */
    void addObject( const float x, const float z, const ObjectType type, const float y = OBJECT_ORIGIN_HEIGHT );

    /**
     * @brief updateObjects Update the objects' positions
     * @param dt the time step
     */
    void updateObjects( float dt );

    /**
     * @brief resetObjects Delete the objects
     */
    void resetObjects();

private slots:
    /** Callback function, will be called whenever the timer ticks*/
    void tick();

};

#endif // VIEW_H

