/** glwidget.cpp
 ** Brief: This is the source file of glwidget.h.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "GL/glut.h"
#include "glwidget.h"
#include "types.h"
#include "fluid.h"
#include "terrain.h"
#include "random_terrain.h"
#include "heightmap_terrain.h"
#include <QApplication>
#include <QKeyEvent>
//added by SH
#include <QGLFramebufferObject>
#include <QGLShaderProgram>

// Declaration of Cuda functions
extern "C"
{
    void testVector();
//    extern void APIENTRYP glActiveTexture(GLenum);
}

/**
 * Local variables in this cpp scope
 */
static Colorf clearColor = Colorf(0.f,0.f,0.f,0.f);

GLWidget::GLWidget(QWidget *parent) : QGLWidget(parent)
  ,m_timer(this),m_prevTime(0), m_prevFps(0.f), m_fps(0.f),
    m_font("Deja Vu Sans Mono", 8,  4)
{
    init();
}

GLWidget::~GLWidget()
{
    if( m_terrain )
        delete m_terrain;
    if( m_fluid )
        delete m_fluid;

    foreach (QGLShaderProgram *sp, m_shaderPrograms)
        delete sp;
    foreach (QGLFramebufferObject *fbo, m_framebufferObjects)
        delete fbo;
}

void GLWidget::init()
{
    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));

    m_fluid =  new Fluid();
#ifdef USE_HEIGHTMAP
    m_terrain = new HeightmapTerrain(); //added by hcreynol
#else
    m_terrain = new Terrain();
#endif

    // Start a timer that will try to get 60 frames per second (the actual
    // frame rate depends on the operating system and other running programs)
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

    glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);   // Always reset the screen to black before drawing anything
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
    timeUpdate();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

#ifdef DRAW_TERRAIN
    m_terrain->draw();
#endif

#ifdef RENDER_FLUID
    // Fluid part
    m_fluid->update( 0.02 );
    m_fluid->draw();
#endif
    //The lighting stuff - SH
    int width = this->width();
    int height = this->height();


    // Render the scene to a framebuffer
//    m_framebufferObjects["fbo_0"]->bind();
//    applyPerspectiveCamera(width, height);
//    renderScene();
//    m_framebufferObjects["fbo_0"]->release();

//    // Copy the rendered scene into framebuffer 1
//    m_framebufferObjects["fbo_0"]->blitFramebuffer(m_framebufferObjects["fbo_1"],
//                                                   QRect(0, 0, width, height), m_framebufferObjects["fbo_0"],
//                                                   QRect(0, 0, width, height), GL_COLOR_BUFFER_BIT, GL_NEAREST);

//    // TODO: Step 0 - draw the scene to the screen as a textured quad
//    applyOrthogonalCamera(width, height);
//    glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_1"]->texture());
//    renderTexturedQuad(width, height);
//    glBindTexture(GL_TEXTURE_2D, 0);

//    // TODO: Step 1 - use the brightpass shader to render bright areas
//    // only to fbo_2
//    m_framebufferObjects["fbo_2"]->bind();
//    m_shaderPrograms["brightpass"]->bind();
//    glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_1"]->texture());
//    renderTexturedQuad(width, height);
//    m_shaderPrograms["brightpass"]->release();
//    glBindTexture(GL_TEXTURE_2D, 0);
//    m_framebufferObjects["fbo_2"]->release();

//    // TODO: Uncomment this section in step 2 of the lab
//    float scales[] = {4.f,8.f};
//    for (int i = 0; i < 2; ++i)
//    {
//        // Render the blurred brightpass filter result to fbo 1
//       renderBlur(width / scales[i], height / scales[i]);

//       // Bind the image from fbo to a texture
//        glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_1"]->texture());
//        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

//        // Enable alpha blending and render the texture to the screen
//        glEnable(GL_BLEND);
//        glBlendFunc(GL_ONE, GL_ONE);
//        renderTexturedQuad(width * scales[i], height * scales[i]);
//        glDisable(GL_BLEND);
//        glBindTexture(GL_TEXTURE_2D, 0);
//    }

    paintText();
}

/**
  Renders the scene.  May be called multiple times by paintGL() if necessary.
**/
void GLWidget::renderScene()
{
    // Enable depth testing
//    glEnable(GL_DEPTH_TEST);
//    glClear(GL_DEPTH_BUFFER_BIT);

    // Enable cube maps and draw the skybox
//    glEnable(GL_TEXTURE_CUBE_MAP);
//    glBindTexture(GL_TEXTURE_CUBE_MAP, m_cubeMap);
//    glCallList(m_skybox);

    // Enable culling (back) faces for rendering the fluid
//    glEnable(GL_CULL_FACE);

    // Render the fluid with the refraction shader bound
//    glActiveTexture(GL_TEXTURE0);
//    m_shaderPrograms["refract"]->bind();
//    m_shaderPrograms["refract"]->setUniformValue("CubeMap", GL_TEXTURE0);
////    glPushMatrix();
////    glTranslatef(-1.25f, 0.f, 0.f);
////    glCallList(m_dragon.idx);
////    glPopMatrix();
//    m_shaderPrograms["refract"]->release();

//    // Render the fluid with the reflection shader bound
//    m_shaderPrograms["reflect"]->bind();
//    m_shaderPrograms["reflect"]->setUniformValue("CubeMap", GL_TEXTURE0);
////    glPushMatrix();
////    glTranslatef(1.25f,0.f,0.f);
////    glCallList(m_dragon.idx);
////    glPopMatrix();
//    m_shaderPrograms["reflect"]->release();
    // Disable culling, depth testing and cube maps
//    glDisable(GL_CULL_FACE);
//    glDisable(GL_DEPTH_TEST);
//    glBindTexture(GL_TEXTURE_CUBE_MAP,0);
//    glDisable(GL_TEXTURE_CUBE_MAP);
}


void GLWidget::resizeGL(int w, int h)
{
    if (w < 1) w = 1;
    if (h < 1) h = 1;

    m_camera.setRatio((float)w/(float)h);
        updateCamera();
    glViewport(0, 0, w, h);
}

/**
  Load a cube map for the skybox
 **/
void GLWidget::loadCubeMap()
{
//    QList<QFile *> fileList;
//    fileList.append(new QFile("/course/cs123/bin/textures/astra/posx.jpg"));
//    fileList.append(new QFile("/course/cs123/bin/textures/astra/negx.jpg"));
//    fileList.append(new QFile("/course/cs123/bin/textures/astra/posy.jpg"));
//    fileList.append(new QFile("/course/cs123/bin/textures/astra/negy.jpg"));
//    fileList.append(new QFile("/course/cs123/bin/textures/astra/posz.jpg"));
//    fileList.append(new QFile("/course/cs123/bin/textures/astra/negz.jpg"));
//    m_cubeMap = ResourceLoader::loadCubeMap(fileList);
}

/**
  Create shader programs.
 **/
void GLWidget::createShaderPrograms()
{
    const QGLContext *ctx = context();
    m_shaderPrograms["reflect"] = ResourceLoader::newShaderProgram(ctx, "shaders/reflect.vert", "shaders/reflect.frag");
    m_shaderPrograms["refract"] = ResourceLoader::newShaderProgram(ctx, "shaders/refract.vert", "shaders/refract.frag");
    m_shaderPrograms["brightpass"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/brightpass.frag");
    m_shaderPrograms["blur"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/blur.frag");
}

/**
  Allocate framebuffer objects.

  @param width: the viewport width
  @param height: the viewport height
 **/
void GLWidget::createFramebufferObjects(int width, int height)
{
    // Allocate the main framebuffer object for rendering the scene to
    // This needs a depth attachment
    m_framebufferObjects["fbo_0"] = new QGLFramebufferObject(width, height, QGLFramebufferObject::Depth,
                                                             GL_TEXTURE_2D, GL_RGB16F_ARB);
    m_framebufferObjects["fbo_0"]->format().setSamples(16);
    // Allocate the secondary framebuffer obejcts for rendering textures to (post process effects)
    // These do not require depth attachments
    m_framebufferObjects["fbo_1"] = new QGLFramebufferObject(width, height, QGLFramebufferObject::NoAttachment,
                                                             GL_TEXTURE_2D, GL_RGB16F_ARB);

    m_framebufferObjects["fbo_2"] = new QGLFramebufferObject(width, height, QGLFramebufferObject::NoAttachment,
                                                             GL_TEXTURE_2D, GL_RGB16F_ARB);
}

/**
  Run a gaussian blur on the texture stored in fbo 2 and
  put the result in fbo 1.  The blur should have a radius of 2.

  @param width: the viewport width
  @param height: the viewport height
**/
void GLWidget::renderBlur(int width, int height)
{
    int radius = 2;
    int dim = radius * 2 + 1;
    GLfloat kernel[dim * dim];
    GLfloat offsets[dim * dim * 2];
    createBlurKernel(radius, width, height, &kernel[0], &offsets[0]);

    // TODO: Step 2 - Finish filling this in
    m_framebufferObjects["fbo_1"]->bind();

    m_shaderPrograms["blur"]->bind();
    m_shaderPrograms["blur"]->setUniformValue("arraySize", dim * dim);
    m_shaderPrograms["blur"]->setUniformValueArray("offsets", offsets, dim * dim * 2, 2);
    m_shaderPrograms["blur"]->setUniformValueArray("kernel", kernel, dim * dim, 1);

    glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_2"]->texture());

    renderTexturedQuad(width, height);

    m_shaderPrograms["blur"]->release();;
    glBindTexture(GL_TEXTURE_2D, 0);
    m_framebufferObjects["fbo_1"]->release();
}

/**
Creates a gaussian blur kernel with the specified radius.  The kernel values
and offsets are stored.

@param radius: The radius of the kernel to create.
@param width: The width of the image.
@param height: The height of the image.
@param kernel: The array to write the kernel values to.
@param offsets: The array to write the offset values to.
**/
void GLWidget::createBlurKernel(int radius, int width, int height,
                                                  GLfloat* kernel, GLfloat* offsets)
{
  int size = radius * 2 + 1;
  float sigma = radius / 3.0f;
  float twoSigmaSigma = 2.0f * sigma * sigma;
  float rootSigma = sqrt(twoSigmaSigma * M_PI);
  float total = 0.0f;
  float xOff = 1.0f / width, yOff = 1.0f / height;
  int offsetIndex = 0;
  for (int y = -radius, idx = 0; y <= radius; ++y)
  {
      for (int x = -radius; x <= radius; ++x,++idx)
      {
          float d = x * x + y * y;
          kernel[idx] = exp(-d / twoSigmaSigma) / rootSigma;
          total += kernel[idx];
          offsets[offsetIndex++] = x * xOff;
          offsets[offsetIndex++] = y * yOff;
      }
  }
  for (int i = 0; i < size * size; ++i)
  {
      kernel[i] /= total;
  }
}

/**
  Draws a textured quad. The texture must be bound and unbound
  before and after calling this method - this method assumes that the texture
  has been bound beforehand.

  @param w: the width of the quad to draw
  @param h: the height of the quad to draw
**/
void GLWidget::renderTexturedQuad(int width, int height) {
    // Clamp value to edge of texture when texture index is out of bounds
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // Draw the  quad
    glBegin(GL_QUADS);
    glTexCoord2f(0.0f, 0.0f);
    glVertex2f(0.0f, 0.0f);
    glTexCoord2f(1.0f, 0.0f);
    glVertex2f(width, 0.0f);
    glTexCoord2f(1.0f, 1.0f);
    glVertex2f(width, height);
    glTexCoord2f(0.0f, 1.0f);
    glVertex2f(0.0f, height);
    glEnd();
}

/**
  Called to switch to an orthogonal OpenGL camera.
  Useful for rending a textured quad across the whole screen.

  @param width: the viewport width
  @param height: the viewport height
**/
void GLWidget::applyOrthogonalCamera(float width, float height)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.f, width, 0.f, height);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

/**
  Called to switch to a perspective OpenGL camera.

  @param width: the viewport width
  @param height: the viewport height
**/
void GLWidget::applyPerspectiveCamera(float width, float height)
{
    float ratio = ((float) width) / height;
    Vector3 dir(-Vector3::fromAngles(m_camera.m_theta, m_camera.m_phi));
    Vector3 eye(m_camera.m_center - dir * m_camera.m_zoom);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(m_camera.m_fovy, ratio, 0.1f, 1000.f);
    gluLookAt(eye.x, eye.y, eye.z, eye.x + dir.x, eye.y + dir.y, eye.z + dir.z,
              m_camera.m_up.x, m_camera.m_up.y, m_camera.m_up.z);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
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
#ifdef RENDER_FLUID
    m_fluid->addRandomDrop();
#endif
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
        if( m_fluid->isRenderingNormal() )
        {
            m_fluid->disableNormal();
        }
        else
        {
            m_fluid->enableNormal();
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

/**
 * Render the text. If you want to add more texts, put it here
 **/
void GLWidget::paintText()
{
    // Combine the previous and current framerate
    if (m_fps >= 0 && m_fps < 1000)
    {
       m_prevFps *= 0.95f;
       m_prevFps += m_fps * 0.05f;
    }

    glColor3f(1.f,1.f,1.f);
    // QGLWidget's renderText takes xy coordinates, a string, and a font
    renderText(10, 20, "FPS: " + QString::number((int) (m_prevFps)), m_font);
    renderText(10, 35, "Delta: " + QString::number((float) (m_delta)), m_font);
    renderText(10, 50, "S: Save screenshot", m_font);
}

/**
 * Update the time variables
 **/
void GLWidget::timeUpdate()
{
    long time = m_time.elapsed();
    m_delta = (time-m_prevTime)/1000.f;
    m_fps = 1000.f /  (time-m_prevTime);
    m_prevTime = time;
}

void GLWidget::tick()
{
    // Get the number of seconds since the last tick (variable update rate)
 //   float seconds = m_time.restart() * 0.001f;

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
