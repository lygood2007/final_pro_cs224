/** glwidget.cpp
 ** Brief: This is the source file of glwidget.h.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "GL/glut.h"
#include "glwidget.h"
#include "CS123Algebra.h"
#include "types.h"
#include "fluid.h"
#include "utils.h"
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
    extern void glActiveTexture(GLenum);
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
    glDeleteLists(m_skybox, 1);
    const_cast<QGLContext *>(context())->deleteTexture(m_cubeMap);
}

void GLWidget::init()
{
    // View needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // View needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // The game loop is implemented using a timer
    connect(&m_timer, SIGNAL(timeout()), this, SLOT(tick()));


#ifdef USE_HEIGHTMAP
    m_terrain = new HeightmapTerrain(); //added by hcreynol
#else
    m_terrain = new RandomTerrain();
#endif
    m_fluid =  new Fluid(m_terrain);
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
    m_camera.applyPerspectiveCamera(WIN_W,WIN_H);
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

    initializeResources();
}

/** Putting all of our other inits in one place*/
void GLWidget::initializeResources()
{
    cout << "--- Loading Resources ---" << endl;

    m_terrain->generate();
    cout << "  Generated Terrain ->" << endl;

    m_fluid->backupHeight(m_terrain);
    cout << "  Calculated Fluid Height ->" << endl;

    m_skybox = ResourceLoader::loadSkybox();
    loadCubeMap();
    cout << "  Loaded Skymap ->" << endl;

    createShaderPrograms();
    cout << "  Loaded Shaders ->" << endl;

    createFramebufferObjects(width(), height());
    cout << "  Loaded FBO's->" << endl;

    cout << " --- Finish Loading Resources ---" << endl;
}


void GLWidget::paintGL()
{
    int width = this->width();
    int height = this->height();
    timeUpdate();

    updateCamera();
    m_camera.applyPerspectiveCamera(width,height);
    m_framebufferObjects["fbo_0"]->bind();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    renderScene();
    m_framebufferObjects["fbo_0"]->release();

    // Copy the rendered scene into framebuffer 1
    m_framebufferObjects["fbo_0"]->blitFramebuffer(m_framebufferObjects["fbo_1"],
                                                   QRect(0, 0, width, height), m_framebufferObjects["fbo_0"],
                                                   QRect(0, 0, width, height), GL_COLOR_BUFFER_BIT, GL_NEAREST);

    //draw the scene to the screen as a textured quad
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDisable(GL_LIGHTING);
    glEnable(GL_TEXTURE_2D);
    glViewport(0,0,width,height);
    applyOrthogonalCamera(width, height);
    glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_0"]->texture());
    renderTexturedQuad(width, height);
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);


    // use the brightpass shader to render bright area only to fbo_2 for bloom effects
    m_framebufferObjects["fbo_2"]->bind();
    m_shaderPrograms["brightpass"]->bind();
    glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_1"]->texture());
    renderTexturedQuad(width, height);
    m_shaderPrograms["brightpass"]->release();
    glBindTexture(GL_TEXTURE_2D, 0);
    m_framebufferObjects["fbo_2"]->release();

//      running a blurring effect over the bright areas
    float scales[] = {4.f,8.f};
    for (int i = 0; i < 2; ++i)
    {
        // Render the blurred brightpass filter result to fbo 1
       renderBlur(width / scales[i], height / scales[i]);

       // Bind the image from fbo to a texture
        glBindTexture(GL_TEXTURE_2D, m_framebufferObjects["fbo_1"]->texture());
        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//        glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

        // Enable alpha blending and render the texture to the screen
        glEnable(GL_BLEND);
        glDisable(GL_LIGHTING);
        glEnable(GL_TEXTURE_2D);
        glViewport(0,0,width,height);
        glBlendFunc(GL_ONE, GL_ONE);
        renderTexturedQuad(width * scales[i], height * scales[i]);
        glDisable(GL_BLEND);
        glBindTexture(GL_TEXTURE_2D, 0);
        glDisable(GL_TEXTURE_2D);
        glEnable(GL_LIGHTING);
    }

    paintText();
}

/**
 * @brief GLWidget::renderSkybox - renders a skybox that is always centered around the camera
 */
void GLWidget::renderSkybox()
{
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);;
    glDisable(GL_LIGHTING); //so the map will be uniformly bright
    //I should just be able to ask the camera it's position but that doesn't seem to work, so this
//    Vector4 temp = m_camera.getEyePos(); Vector3 eye = Vector3(temp.x, temp.y, temp.z);
    Vector3 dir(-Vector3::fromAngles(m_camera.m_theta, m_camera.m_phi));
    Vector3 eye(m_camera.m_center - dir * m_camera.m_zoom);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(eye.x,eye.y, eye.z); //keeps the skybox centered around the camera

    // Enable cube maps and draw the skybox
    glEnable(GL_TEXTURE_CUBE_MAP);
    glBindTexture(GL_TEXTURE_CUBE_MAP, m_cubeMap);
    glCallList(m_skybox);
    glPopMatrix();

    //turn the lights back on and release textures
    glBindTexture(GL_TEXTURE_CUBE_MAP,0);
    glDisable(GL_TEXTURE_CUBE_MAP);
    glEnable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
}

/**
  Renders the visible geometry, terrain and fluid only at this point
**/
void GLWidget::renderFluid()
{

#ifdef RENDER_FLUID
    // Fluid part
    m_fluid->update( 0.02 );
    m_fluid->draw();
#endif


}


/**
  Renders the scene.  May be called multiple times by paintGL() if necessary.
**/
void GLWidget::renderScene()
{

    renderSkybox();//@NOTE - This must go first!!

#ifdef DRAW_TERRAIN
   m_terrain->draw();
#endif

   if(false) // make true to draw some axis'
   {
       static GLUquadric * quad = gluNewQuadric();
       glColor3f(0, 0, 1);
       gluCylinder(quad, 1, 1, 10, 10, 10); // Z
       glPushMatrix();
       glRotatef(90, 0, 1, 0);
       glColor3f(1, 0, 0);
       gluCylinder(quad, 1, 1, 10, 10, 10); // X
       glPopMatrix();
       glPushMatrix();
       glRotatef(-90, 1, 0, 0);
       glColor3f(0, 1, 0);
       gluCylinder(quad, 1, 1, 10, 10, 10); // Y
       glPopMatrix();
    }

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // Enable depth testing

    //glEnable(GL_DEPTH_TEST);
    // Enable culling (back) faces for rendering the fluid and terrain
 //   glEnable(GL_CULL_FACE);
    glBindTexture(GL_TEXTURE_CUBE_MAP, m_cubeMap);

//     Render the fluid with the refraction shader bound
//    glActiveTexture(GL_TEXTURE0);
//    m_shaderPrograms["refract"]->bind();
//    m_shaderPrograms["refract"]->setUniformValue("CubeMap", GL_TEXTURE0);
//    glPushMatrix();
//    glTranslatef(-1.25f, 0.f, 0.f);
//    renderFluid();
//    glPopMatrix();
//    m_shaderPrograms["refract"]->release();

    if(false) //true for perfect reflection, false for fresnel
    {
        // Render the fluid with the reflection shader bound
        m_shaderPrograms["reflect"]->bind();
        m_shaderPrograms["reflect"]->setUniformValue("CubeMap", GL_TEXTURE0);
        m_shaderPrograms["reflect"]->setUniformValue("CurrColor", 0.1f,0.4f,0.8f,0.5f);
        glPushMatrix();
        glTranslatef(1.25f,0.f,0.f);
        renderFluid();
        glPopMatrix();
        m_shaderPrograms["reflect"]->release();
    }
    else
    {
        // Render the fluid with the fresnel shader bound for reflection and refraction
        m_shaderPrograms["fresnel"]->bind();
        m_shaderPrograms["fresnel"]->setUniformValue("CubeMap", GL_TEXTURE0);
        m_shaderPrograms["fresnel"]->setUniformValue("CurrColor", 0.1f,0.4f,0.8f,1.0f);
        glPushMatrix();
        glTranslatef(1.25f,0.f,0.f);
        renderFluid();
        glPopMatrix();
        m_shaderPrograms["fresnel"]->release();
    }


//     Disable culling, depth testing and cube maps
    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
    glDisable(GL_TEXTURE_CUBE_MAP);
    //glDisable(GL_CULL_FACE);
   // glDisable(GL_DEPTH_TEST);
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
    QList<QFile *> fileList;

    if(true) //true for real, false for testing color reflections
    {
        fileList.append(new QFile("resource/posx.jpg"));
        fileList.append(new QFile("resource/negx.jpg"));
        fileList.append(new QFile("resource/posy.jpg"));
        fileList.append(new QFile("resource/sandy_sea_floor.jpg"));
        fileList.append(new QFile("resource/posz.jpg"));
        fileList.append(new QFile("resource/negz.jpg"));
    }
    else
    {
        fileList.append(new QFile("resource/red.jpg"));
        fileList.append(new QFile("resource/orange.jpg"));
        fileList.append(new QFile("resource/green.jpg"));
        fileList.append(new QFile("resource/sandy_sea_floor.jpg"));
        fileList.append(new QFile("resource/purple.jpg"));
        fileList.append(new QFile("resource/yellow.jpg"));

    }
    m_cubeMap = ResourceLoader::loadCubeMap(fileList);
}

/**
  Create shader programs.
 **/
void GLWidget::createShaderPrograms()
{
    const QGLContext *ctx = context();
    m_shaderPrograms["reflect"] = ResourceLoader::newShaderProgram(ctx, "shaders/reflect.vert", "shaders/reflect.frag");
    m_shaderPrograms["refract"] = ResourceLoader::newShaderProgram(ctx, "shaders/refract.vert", "shaders/refract.frag");
    m_shaderPrograms["fresnel"] = ResourceLoader::newShaderProgram(ctx, "shaders/fresnel.vert", "shaders/fresnel.frag");
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
    //m_framebufferObjects["fbo_0"]->format().setSamples(16);
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
    glColor3f(1.0,1.0,1.0);
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

//    Vector4 temp = m_camera.getEyePos();
//    Vector3 eye = Vector3(temp.x, temp.y, temp.z);

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
         intersectFluid( event->x(), event->y() );
    //m_fluid->addRandomDrop();
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

/**
 * @brief intersectFluid Check if the ray shooting from position (x, y)  intersects the fluid
 * @param x, The x position in screen space
 * @param y, The y position in screen space
 * @return Return if it is intersected
 */
void GLWidget::intersectFluid( const int x, const int y)
{
    Vector4 eyePos = m_camera.getEyePos();
    Vector4 pFilmCam;
    Matrix4x4 invViewTransMat = m_camera.getInvViewTransMatrix();
    pFilmCam.x = ((REAL)(2*x))/width() - 1; pFilmCam.y = 1- ((REAL)(2*y))/height(); pFilmCam.z = -1;
    pFilmCam.w = 1;

    Vector4 pFilmWorld = invViewTransMat*pFilmCam;
    Vector4 d = pFilmWorld - eyePos;
    d = d.getNormalized();
    Vector3 dir3 = Vector3(d.x,d.y,d.z);
    Vector3 eye3 = Vector3(eyePos.x,eyePos.y,eyePos.z);
    float halfDomain = m_fluid->m_domainSize/2;
    float dx = m_fluid->m_dx;

    const QVector<Tri>& temp = m_fluid->m_triangles;
    const QVector<QVector<float> >& tempHeight = m_fluid->m_terrainHeightField;
    const QVector<QVector<float> >& tempDepth = m_fluid->m_depthField;
    int indexRow = -1;
    int indexCol = -1;
    const int gridSize = m_fluid->m_gridSize;

    for( int i = 0; i < temp.size(); i++ )
    {
        Tri curTri = temp[i];
        // Firstly check if the triangle is visible
        int count = 0;
        int r[3]; int c[3];
        r[0] = curTri.a2D.indRow; c[0] = curTri.a2D.indCol;
        r[1] = curTri.b2D.indRow; c[1] = curTri.b2D.indCol;
        r[2] = curTri.c2D.indRow; c[2] = curTri.c2D.indCol;
        for( int m = 0; m < 3; m++ )
        {
            if( m_fluid->m_depthField[r[m]][c[m]] > EPSILON )
                count++;
        }
        if( count == 0 )
            continue;
        else
        {
            const Vector3 p0 = Vector3(-halfDomain + c[0]*dx, tempHeight[r[0]][c[0]] + tempDepth[r[0]][c[0]],
                                      - halfDomain + r[0]*dx );
            const Vector3 p1 = Vector3(-halfDomain + c[1]*dx, tempHeight[r[1]][c[1]] + tempDepth[r[1]][c[1]],
                                      - halfDomain + r[1]*dx );
            const Vector3 p2 = Vector3(-halfDomain + c[2]*dx, tempHeight[r[2]][c[2]] + tempDepth[r[2]][c[2]],
                                      - halfDomain + r[2]*dx );

             if( doIntersectTriangles( eye3, dir3, p0, p1, p2 ) )
             {
                 indexRow =  gridSize - c[0];
                 indexCol = r[0];
                 break;
             }
        }
    }

    if( indexRow != -1 && indexCol != -1 )
    {
        m_fluid->addDrop( indexCol, indexRow );
    }
}

void GLWidget::tick()
{
    // Get the number of seconds since the last tick (variable update rate)
 //   float seconds = m_time.restart() * 0.001f;

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
