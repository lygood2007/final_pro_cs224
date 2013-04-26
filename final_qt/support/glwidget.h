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
#include "resourceloader.h"

// We fix the size
#define WIN_W 800.0
#define WIN_H 800.0

// Flag for testing
#define DRAW_TERRAIN

//added by hcreynol
//#define USE_HEIGHTMAP

#define RENDER_FLUID

class QGLShaderProgram;
class QGLFramebufferObject;
class Terrain;
class Fluid;

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent);
    ~GLWidget();

private:
    // Private variables
    QTime m_time; // The time variable
    QTimer m_timer; // The timer variable
    OrbitCamera m_camera; // Camera
    Terrain* m_terrain;
    Fluid* m_fluid;

    QHash<QString, QGLShaderProgram *> m_shaderPrograms; // hash map of all shader programs
    QHash<QString, QGLFramebufferObject *> m_framebufferObjects; // hash map of all framebuffer objects
    GLuint m_skybox; // skybox call list ID
    GLuint m_cubeMap; // cubeMap texture ID
    QFont m_font; // font for rendering tex


    bool m_mouseLeftDown; // True if mouse left is down
    bool m_mouseRightDown; // True if mouse right is down

    bool m_drawFrame; // True if draw in wireframe mode

    int m_prevTime;
    float m_prevFps, m_fps;
    float m_delta;
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
    void renderScene();
    void renderTexturedQuad(int width, int height);
    void renderSkybox();

    /** Placeing all visible geometry rendering in one method*/
    void renderFluid();

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
    void intersectFluid( const int x, const int y);

    //copied from CS123 lab 09 - SH
    void loadCubeMap();
    void createShaderPrograms();
    void createFramebufferObjects(int width, int height);
    void createBlurKernel(int radius, int width, int height, GLfloat* kernel, GLfloat* offsets);

    void renderBlur(int width, int height);


private slots:
    /** Callback function, will be called whenever the timer ticks*/
    void tick();
};

#endif // VIEW_H

