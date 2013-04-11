/** glwidget.cpp
 ** Brief: This is the source file of glwidget.h.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "GL/glut.h"
#include "glwidget.h"
#include <QApplication>
#include <QKeyEvent>

// Declaration of Cuda functions
extern "C"
{
    void testVector();
}

GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent)
{
    init();
}

GLWidget::~GLWidget()
{
    if( m_terrain )
        delete m_terrain;
}

void GLWidget::init()
{
    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));

    m_terrain = new Terrain();

    // Start a timer that will try to get 60 frames per second (the actual
    // frame rate depends on the operating system and other running programs)
    m_time.start();
    m_timer.start(1000 / 60);

    m_mouseLeftDown = false;
    m_mouseRightDown = false;
    m_drawFrame = false;

    // Center the mouse, which is explained more in mouseMoveEvent() below.
    // This needs to be done here because the mouse may be initially outside
    // the fullscreen window and will not automatically receive mouse move
    // events. This occurs if there are two monitors and the mouse is on the
    // secondary monitor.
    QCursor::setPos(mapToGlobal(QPoint(width() / 2, height() / 2)));

    setFixedSize( WIN_W, WIN_H );
    m_camera.setRatio(WIN_W/WIN_H);
    updateCamera();

    /** ONLY A TEST FOR USING CUDA*/
    testVector();
}

void GLWidget::initializeGL()
{
    // All OpenGL initialization *MUST* be done during or after this
    // method. Before this method is called, there is no active OpenGL
    // context and all OpenGL calls have no effect.

    glClearColor(0, 0, 0, 0);   // Always reset the screen to black before drawing anything
    glEnable(GL_DEPTH_TEST);    // When drawing a triangle, only keep pixels closer to the camera than what's already been drawn
     glDisable(GL_DITHER);
     glShadeModel(GL_SMOOTH);
    // Make things pretty
    glEnable(GL_MULTISAMPLE);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);

    // Bind the ambient and diffuse color of each vertex to the current glColor() value
   glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

    // Cull triangles that are facing away from the camera
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Set up a single light
    glEnable(GL_LIGHTING);
    // Light's color
    GLfloat ambientColor[] = { 0.3f, 0.3f, 0.3f, 1.0f };
    GLfloat diffuseColor[] = { 1.0f, 1.0f, 1.0, 1.0f };
    GLfloat specularColor[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat lightPosition[] = { 0.f, 0.f, 10.f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientColor);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseColor);
    glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);
    m_terrain->generate();
    m_camera.applyPerspectiveCamera(WIN_W,WIN_H);
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // TODO: Implement the demo rendering here
#ifdef DRAW_TERRAIN
    m_terrain->draw();
#endif

}

void GLWidget::resizeGL(int w, int h)
{
    if (w < 1) w = 1;
    if (h < 1) h = 1;

    m_camera.setRatio((float)w/(float)h);
        updateCamera();
    glViewport(0, 0, w, h);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    if( event->button() == Qt::RightButton )
    {
        m_camera.mouseDown(event->x(),event->y());
        updateCamera();
        m_mouseRightDown = true;
    }
    else
    {
        m_mouseLeftDown = true;
    }
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    // This starter code implements mouse capture, which gives the change in
    // mouse position since the last mouse movement. The mouse needs to be
    // recentered after every movement because it might otherwise run into
    // the edge of the screen, which would stop the user from moving further
    // in that direction. Note that it is important to check that deltaX and
    // deltaY are not zero before recentering the mouse, otherwise there will
    // be an infinite loop of mouse move events.
    int deltaX = event->x() - width() / 2;
    int deltaY = event->y() - height() / 2;
    if (!deltaX && !deltaY) return;

    if( m_mouseRightDown )
       {
           m_camera.mouseMove(event->x(),event->y());
           updateCamera();

       }
       else if( m_mouseLeftDown )
       {
       }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_mouseLeftDown = false;
    m_mouseRightDown = false;
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    m_camera.mouseWheel(event->delta());
    updateCamera();
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Escape) QApplication::quit();

    // TODO: Handle keyboard presses here
    switch(event->key())
    {
        case Qt::Key_W:
     {
        if( !m_drawFrame )
        {
            glPolygonMode(GL_FRONT, GL_LINE);
            m_drawFrame = true;
        }
        else
        {
            glPolygonMode(GL_FRONT,GL_FILL);
            m_drawFrame = false;
        }
        break;
        }
    case Qt::Key_N:
    {
        if( m_terrain-> isRenderingNormal())
        {
            m_terrain->disableNormal();
        }
        else
        {
            m_terrain->enableNormal();
        }
        break;
    }
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
}

/**
 * Called whenever m_camera is modified or the canvas is resized. Sets the current OpenGL projection
 * and modelview matrices to match the values in m_camera.
 */
void GLWidget::updateCamera()
{
    m_camera.applyPerspectiveCamera(width(),height());
}

void GLWidget::tick()
{
    // Get the number of seconds since the last tick (variable update rate)
    float seconds = m_time.restart() * 0.001f;

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
