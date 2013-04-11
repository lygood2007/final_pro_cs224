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
#include "camera.h"
#include "terrain.h"

// We fix the size
#define WIN_W 800.0
#define WIN_H 800.0

// Flag for testing
#define DRAW_TERRAIN

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

    bool m_mouseLeftDown; // True if mouse left is down
    bool m_mouseRightDown; // True if mouse right is down

    bool m_drawFrame; // True if draw in wireframe mode
private:
    // Private functions
    /** Initliaze the variables in the class*/
    void init();
    /** Initialize the GL states, put all your GL initialization code here */
    void initializeGL();
    /** The render function, will be called every frame */
    void paintGL();
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

private slots:
    /** Callback function, will be called whenever the timer ticks*/
    void tick();
};

#endif // VIEW_H

