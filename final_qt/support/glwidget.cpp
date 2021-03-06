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
#include "fluidCPU.h"
#include "fluidGPU.h"
#include "utils.h"
#include "terrain.h"
#include "sphere.h"
#include "box.h"
#include "random_terrain.h"
#include "heightmap_terrain.h"
#include <QApplication>
#include <QKeyEvent>
//added by SH
#include <QGLFramebufferObject>
#include <QGLShaderProgram>
#include <fstream>
#include <sstream>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>

// Declaration of Cuda functions
extern "C"
{
    void testVector();
    bool findSupportDevice();
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
#ifdef RENDER_FLUID
    if( m_fluid )
        delete m_fluid;
#endif

    // Delete Objects
    foreach( Object* o, m_objects )
    {
        if( o )
        {
            delete o;
            o = NULL;
        }
    }

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

    initConfig();
    //use everything
    m_useShaders = m_useFBO = m_useSimpleCube = m_useSkybox = m_useParticles = true;
    m_useParticleSources = m_useRectangularParticleSources = false;
    //except these
    m_useAxis = false; m_useDampening = false;

    if( m_useHeightMap )
    {
        m_terrain = new HeightmapTerrain(m_heightMapFileName.c_str());
    }
    else
        m_terrain = new RandomTerrain();


    // Start a timer that will try to get 100 frames per second (the actual
    // frame rate depends on the operating system and other running programs)
    m_timer.start(1000 / 60);

    m_mouseLeftDown = false;
    m_mouseRightDown = false;
    m_drawFrame = false;
    m_animate = false;
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

#ifdef RENDER_FLUID
#ifdef USE_GPU_FLUID
        m_fluid = new FluidGPU(m_terrain, this); //I also passed a pointer to glwidget so can control things easier via bool's;
#else
     m_fluid =  new FluidCPU(m_terrain);
#endif
#endif
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
//    GLfloat specularColor[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat lightPosition[] = { 0.f, 0.f, 10.f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambientColor);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuseColor);
//    glLightfv(GL_LIGHT0, GL_SPECULAR, specularColor);
    glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);
    glEnable(GL_LIGHT0);

    glEnable(GL_PROGRAM_POINT_SIZE_EXT); //so we can adjust particle sizes when drawing gl_points

    initializeResources();
}

/** Putting all of our other inits in one place*/
void GLWidget::initializeResources()
{
    cout << "--- Loading Resources ---" << endl;

    m_terrain->generate();
    cout << "  Generated Terrain ->" << endl;
#ifdef RENDER_FLUID
    m_fluid->backupHeight(m_terrain);
    cout << "  Calculated Fluid Height ->" << endl;
#endif
    loadObjectTexMap();
    m_skybox = ResourceLoader::loadSkybox();
    loadCubeMap();
    cout << "  Loaded Skymap ->" << endl;

    createShaderPrograms();
    m_foamTex = loadTexture("./resource/foam.jpg"); //my foam texture to try and smooth it out
    cout << "  Loaded Shaders ->" << endl;

    createFramebufferObjects(width(), height());
    cout << "  Loaded FBO's->" << endl;

    cout << " --- Finish Loading Resources ---" << endl;
}


void GLWidget::paintGL()
{
    int width = this->width();
    int height = this->height();


    if(m_useFBO)
    {
        updateCamera();
        m_camera.applyPerspectiveCamera(width,height);
        m_framebufferObjects["fbo_0"]->bind();
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    renderScene();

    if(m_useFBO)
    {
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
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
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
            glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

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
    }//end use FBO

    paintText(); //update text information in upper left corner

}

/**
 * @brief GLWidget::renderObjects render the objects
 */
void GLWidget::renderObjects()
{
    for( int i = 0; i < m_objects.size(); i++ )
    {
        if( m_objects[i])
            m_objects[i]->draw();
    }
}

/**
  Renders the scene.  May be called multiple times by paintGL() if necessary.
**/
void GLWidget::renderScene()
{
    if(m_useSkybox) renderSkybox();//@NOTE - This must go first!!
    renderObjects();
#ifdef DRAW_TERRAIN
    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST );
  m_terrain->draw();
  glDisable(GL_DEPTH_TEST);
#endif
   if(m_useAxis) // make true to draw some axis'
   {
       static GLUquadric * quad = gluNewQuadric();
       glColor3f(0, 0, 1);
       gluCylinder(quad, 1, 1, 20, 20, 20); // Z
       glPushMatrix();
       glRotatef(90, 0, 1, 0);
       glColor3f(1, 0, 0);
       gluCylinder(quad, 1, 1, 20, 20, 20); // X
       glPopMatrix();
       glPushMatrix();
       glRotatef(-90, 1, 0, 0);
       glColor3f(0, 1, 0);
       gluCylinder(quad, 1, 1, 20, 20, 20); // Y
       glPopMatrix();
    }

   if(m_useSkybox)
   {
       if(m_useShaders){

        //this is entirely awful but needed
        //this takes the correct rotation matrix from the camera and applies it to the
        //cube map so the fresnel shader works correctly
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_CUBE_MAP, m_cubeMap);
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        Vector3 dir(-Vector3::fromAngles(m_camera.m_theta, m_camera.m_phi));
        Vector3 eye( - dir * m_camera.m_zoom);
        gluPerspective(m_camera.m_fovy, (float)WIN_H/WIN_W, m_camera.m_near, m_camera.m_far);
        gluLookAt(eye.x, eye.y, eye.z, eye.x + dir.x, eye.y + dir.y, eye.z + dir.z,
                  m_camera.m_up.x, m_camera.m_up.y, m_camera.m_up.z);
        glTranslatef(eye.x,eye.y, eye.z);
        double matrix[16];
        glGetDoublev(GL_MODELVIEW_MATRIX, matrix);
        glPopMatrix();
        //some of this could maybe be simplified but this spells it out correctly at least
        Matrix4x4 temp = Matrix4x4(matrix);
        Matrix4x4 tempT = temp.getTranspose();
        Matrix4x4 tempI = tempT.getInverse();
        Matrix4x4 tempIT = tempI.getTranspose();
        glMatrixMode(GL_TEXTURE);
        glPushMatrix();
        glLoadMatrixd(tempIT.data);

        // Render the fluid with the fresnel shader bound for reflection and refraction
        m_shaderPrograms["fresnel"]->bind();
        m_shaderPrograms["fresnel"]->setUniformValue("CubeMap", 0);
        m_shaderPrograms["fresnel"]->setUniformValue("CurrColor", SEA_WATER);
        glPushMatrix();
        glTranslatef(0.f,1.25f,0.f);
        renderFluid();
        glPopMatrix();
        glActiveTexture(GL_TEXTURE0);
        m_shaderPrograms["fresnel"]->release();

        //release the cube map
        glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
        glDisable(GL_TEXTURE_CUBE_MAP);
        glPopMatrix();

//        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE); //this should allow us to change point size in the shader but it never quite worked
        // Render the spray with the point shader
        m_shaderPrograms["spray"]->bind();
        m_shaderPrograms["spray"]->setUniformValue("CurrColor", SEA_WATER);
        glPushMatrix();
        glTranslatef(0.f,1.25f,0.f);
        renderSpray();
        glPopMatrix();
        m_shaderPrograms["spray"]->release();

        // Render the splash with the splash shader
        m_shaderPrograms["splash"]->bind();
        m_shaderPrograms["splash"]->setUniformValue("CurrColor", SEA_WATER);
        glPushMatrix();
        glTranslatef(0.f,1.25f,0.f);
        renderSplash();
        glPopMatrix();
        m_shaderPrograms["splash"]->release();


        // Render the foam with the foam shader
        glBindTexture(GL_TEXTURE_2D, m_foamTex);
        m_shaderPrograms["foam"]->bind();
        m_shaderPrograms["foam"]->setUniformValue("texture", GL_TEXTURE0);
        m_shaderPrograms["foam"]->setUniformValue("CurrColor", SEA_WATER);
        glPushMatrix();
        glTranslatef(0.f,1.25f,0.f);
        renderFoam();
        glPopMatrix();
        m_shaderPrograms["foam"]->release();
        glBindTexture(GL_TEXTURE_2D, 0);

        }
        else //plain old fluid, nothing special
        {
            renderFluid();
            renderParticles();
        }

    }
   else //plain old fluid, nothing special
   {
       renderFluid();
       renderParticles();
   }
}

/**
 * @brief GLWidget::renderSkybox - renders a skybox that is always centered around the camera
 */
void GLWidget::renderSkybox()
{
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
    glDisable(GL_LIGHTING); //so the map will be uniformly bright
    //I should just be able to ask the camera it's position but that doesn't seem to work, so this
    Vector3 dir(-Vector3::fromAngles(m_camera.m_theta, m_camera.m_phi));
    Vector3 eye(m_camera.m_center - dir * m_camera.m_zoom);

    glMatrixMode(GL_MODELVIEW);
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
    glEnable(GL_CULL_FACE);
//    glDisable(GL_DEPTH_TEST);
}

/**
  Renders the fluid only
**/
void GLWidget::renderFluid()
{

    //the GL_BLEND here allows for slightly transparent water
    //ultimately I'd like to have the water opacity based on the depth
    glEnable (GL_BLEND);
    glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

#ifdef RENDER_FLUID
    glEnable(GL_DEPTH_TEST );
    m_fluid->draw();
    glDisable(GL_DEPTH_TEST );
#endif

}

/**
    Render all the particles
**/
void GLWidget::renderParticles()
{
#ifdef RENDER_FLUID
    if(m_useParticles)
    {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //        glBlendFunc(GL_ONE,GL_ONE);
        glEnable(GL_DEPTH_TEST );
//        m_fluid->drawParticles2();
        m_fluid->drawSpray();
        m_fluid->drawSplash();
        m_fluid->drawFoam();
        glDisable(GL_DEPTH_TEST );
    }
}

/**
    Render just spray particles
**/
void GLWidget::renderSpray()
{
    if(m_useParticles)
    {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //        glBlendFunc(GL_ONE,GL_ONE);
        glEnable(GL_DEPTH_TEST );
        m_fluid->drawSpray();
        glDisable(GL_DEPTH_TEST );
    }
#endif
}

/**
    Render just splash particles
**/
void GLWidget::renderSplash()
{
    if(m_useParticles)
    {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //        glBlendFunc(GL_ONE,GL_ONE);
        glEnable(GL_DEPTH_TEST );
        m_fluid->drawSplash();
        glDisable(GL_DEPTH_TEST );
    }
}

/**
    Render just foam particles
**/
void GLWidget::renderFoam()
{
    if(m_useParticles)
    {
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //        glBlendFunc(GL_ONE,GL_ONE);
        glEnable(GL_DEPTH_TEST );
        m_fluid->drawFoam();
        glDisable(GL_DEPTH_TEST );
    }
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

    if(m_useSimpleCube) //true for real, false for testing color reflections
    {
        fileList.append(new QFile("resource/posx.jpg"));
        fileList.append(new QFile("resource/negx.jpg"));
        fileList.append(new QFile("resource/posy.jpg"));
        fileList.append(new QFile("resource/rock_texture.jpg"));
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
 * @brief GLWidget::loadObjectTexMap Load the object's texture map into memory
 */
void GLWidget::loadObjectTexMap()
{
   m_boxTexID = loadTexture("resource/box_tex.jpg");
   m_sphereTexID = loadTexture("resource/rock_tex.jpg");
}

/**
  Create shader programs.
 **/
void GLWidget::createShaderPrograms()
{
    const QGLContext *ctx = context();
//    m_shaderPrograms["reflect"] = ResourceLoader::newShaderProgram(ctx, "shaders/reflect.vert", "shaders/reflect.frag");
//    m_shaderPrograms["refract"] = ResourceLoader::newShaderProgram(ctx, "shaders/refract.vert", "shaders/refract.frag");
    m_shaderPrograms["fresnel"] = ResourceLoader::newShaderProgram(ctx, "shaders/f2.vert", "shaders/f2.frag");
    m_shaderPrograms["brightpass"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/brightpass.frag");
    m_shaderPrograms["blur"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/blur.frag");
    m_shaderPrograms["spray"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/spray.frag");
    m_shaderPrograms["splash"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/splash.frag");
    m_shaderPrograms["foam"] = ResourceLoader::newFragShaderProgram(ctx, "shaders/foam.frag");
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
    else if ( event->button() == Qt::MiddleButton )
    {
        m_camera.mouseDown(event->x(),event->y());
        updateCamera();
        m_mouseMiddleDown = true;
    }
    else
    {
#ifdef RENDER_FLUID
        if( m_animate )
        {
            int indexRow;
            int indexCol;
            Vector3 pos;
            if(intersectFluid( event->x(), event->y(), indexRow, indexCol,pos ))
            {
                if(event->button() == Qt::LeftButton){
                    m_fluid->addDrop( indexCol, indexRow );
                } else {
#ifdef USE_GPU_FLUID
                    m_fluid->addDroppingParticles( indexCol, indexRow );
#endif
                }
            }
        }
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
    else if (m_mouseMiddleDown)
    {
        m_camera.mouseMovePan(event->x(),event->y());
        updateCamera();
    }
   else if( m_mouseLeftDown )
   {
   }
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_mouseLeftDown = false;
    m_mouseMiddleDown = false;
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

    switch(event->key())
    {
        case Qt::Key_W:
     {
        if( !m_drawFrame )
        {
            glPolygonMode(GL_FRONT, GL_LINE);
            m_drawFrame = true;
            m_useFBO = false; //cannot draw wireframes if FBO is on

        }
        else
        {
            glPolygonMode(GL_FRONT,GL_FILL);
            m_drawFrame = false;
            m_useFBO = true;
        }
        //update and repaint everything
        updateCamera();
        m_camera.applyPerspectiveCamera(WIN_W,WIN_H);
        paintGL();
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
#ifdef RENDER_FLUID
        if( m_fluid->isRenderingNormal() )
        {            m_fluid->disableNormal();
        }
        else
        {
            m_fluid->enableNormal();
        }
#endif
        for( int i = 0; i < m_objects.size(); i++ )
        {
            if( m_objects[i] )
            {
                if( m_objects[i]->isRenderingNormal() )
                {
                    m_objects[i]->disableNormal();
                }
                else
                {
                    m_objects[i]->enableNormal();
                }
            }
        }
        break;
    }
    case Qt::Key_R:
      {
            m_animate = !m_animate;
            break;
        }
    case Qt::Key_C:
    {
        m_useSkybox = !m_useSkybox;
        break;
    }
    case Qt::Key_S:
    {
        m_useShaders = !m_useShaders;
        break;
    }
    case Qt::Key_A:
    {
        m_useAxis = !m_useAxis;
        break;
    }
    case Qt::Key_F:
    {
        m_useFBO = !m_useFBO;
        //update and repaint everything
        updateCamera();
        m_camera.applyPerspectiveCamera(WIN_W,WIN_H);
        paintGL();
        break;
    }
    case Qt::Key_X:
    {
        m_useSimpleCube = !m_useSimpleCube;
        loadCubeMap(); //need to reload the textures
        break;
    }
    case Qt::Key_P:
    {
        m_useParticles = !m_useParticles;
        break;
    }
    case Qt::Key_O: //add an box
    {
        addObject( 0, 0, BOX );
        break;
    }
    case Qt::Key_Y:
    {
        addObject( 0, 0, SPHERE );
        break;
    }
    case Qt::Key_U: //make the boxes heavier
    {
        foreach( Object* o, m_objects)
        {
            if( dynamic_cast<Box*>( o ) )
            {
                o->setDensity(o->getDensity()+50);
            }
        }
        break;
    }
    case Qt::Key_I: //reduce all object density
    {
        foreach( Object* o, m_objects)
        {
            if( dynamic_cast<Box*>( o ) )
            {
                o->setDensity(o->getDensity()-50);
            }
        }
        break;
    }
    case Qt::Key_D:
    {
        m_useDampening = !m_useDampening;
        break;
    }
    case Qt::Key_1:
    {
        m_fluid->createWave(1);
        break;
    }
    case Qt::Key_2:
    {
        m_fluid->createWave(2);
        break;
    }
    case Qt::Key_3:
    {
        m_fluid->createWave(3);
        break;
    }
    case Qt::Key_4:
    {
        m_fluid->createWave(4);
        break;
    }
    case Qt::Key_Q:
    {
        m_fluid->resetFluid();
        // Reset the objects
        resetObjects();
        break;
    }
    case Qt::Key_K:
    {
        m_useParticleSources = !m_useParticleSources;
        break;
    }
    case Qt::Key_L:
    {
        m_useRectangularParticleSources = !m_useRectangularParticleSources;
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
    if(!m_animate) renderText(10, 65, "Rendering Off!", m_font);
    if(m_useParticles) renderText(10, 80, "Particles On", m_font);
    if(m_useFBO) renderText(10, 95, "FrameBufers On", m_font);
    if(m_useDampening) renderText(10, 110, "Dampening On", m_font);


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
 * @param x, The x position in screen spacem_fluid->addDrop( event->x(), event->y() );
 * @param y, The y position in screen space
 */
bool GLWidget::intersectFluid( const int x, const int y, int& indexRow, int& indexCol, Vector3& pos )
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

    const QVector<TriIndex>& temp = m_fluid->m_triangles;

#ifdef USE_GPU_FLUID
    const float* tempHeight = m_fluid->m_heightField;
    const float* tempDepth = m_fluid->m_depthField;
#else
    const QVector<QVector<float> >& tempHeight = m_fluid->m_terrainHeightField;
    const QVector<QVector<float> >& tempDepth = m_fluid->m_depthField;
#endif
    indexRow = -1;
    indexCol = -1;
    const int gridSize = m_fluid->m_gridSize;

    for( int i = 0; i < temp.size(); i++ )
    {
        TriIndex curTri = temp[i];
        // Firstly check if the triangle is visible
        //int count = 0;
        int r[3]; int c[3];
        r[0] = curTri.a2D.indRow; c[0] = curTri.a2D.indCol;
        r[1] = curTri.b2D.indRow; c[1] = curTri.b2D.indCol;
        r[2] = curTri.c2D.indRow; c[2] = curTri.c2D.indCol;
      /*  for( int m = 0; m < 3; m++ )
        {
#ifdef USE_GPU_FLUID
            if( tempDepth[m_fluid->getIndex1D(r[m],c[m],DEPTH)] > EPSILON )
                count++;
#else
            if( m_fluid->m_depthField[r[m]][c[m]] > EPSILON )
                count++;
#endif
        }*/
        /*if( count == 0 )
            continue;
        else
        {*/
#ifdef USE_GPU_FLUID
            const Vector3 p0 = Vector3(-halfDomain + c[0]*dx, tempHeight[m_fluid->getIndex1D(r[0],c[0],HEIGHT)]
                                      ,- halfDomain + r[0]*dx );
            const Vector3 p1 = Vector3(-halfDomain + c[1]*dx, tempHeight[m_fluid->getIndex1D(r[1],c[1],HEIGHT)]
                                       ,- halfDomain + r[1]*dx );
            const Vector3 p2 = Vector3(-halfDomain + c[2]*dx, tempHeight[m_fluid->getIndex1D(r[2],c[2],HEIGHT)]
                                      ,- halfDomain + r[2]*dx );
#else
            const Vector3 p0 = Vector3(-halfDomain + c[0]*dx, tempHeight[r[0]][c[0]] + tempDepth[r[0]][c[0]],
                                      - halfDomain + r[0]*dx );
            const Vector3 p1 = Vector3(-halfDomain + c[1]*dx, tempHeight[r[1]][c[1]] + tempDepth[r[1]][c[1]],
                                      - halfDomain + r[1]*dx );
            const Vector3 p2 = Vector3(-halfDomain + c[2]*dx, tempHeight[r[2]][c[2]] + tempDepth[r[2]][c[2]],
                                      - halfDomain + r[2]*dx );
#endif
             if( doIntersectTriangles( eye3, dir3, p0, p1, p2 ) )
             {
                 // ******************************************************************************************************
                 // Here, it's a hack that it shouldn't be like this !
                 // I don't know where is wrong, the indexRow should equal to r[0] and indexCol should
                 // equal to c[0]
                 // I guess the intersect check is wrong
                 // But the result is not correct
                 // ******************************************************************************************************
                 indexRow =  gridSize - c[0];
                 indexCol = r[0];
                 pos = (p0 + p1 + p2)/3.f;
                 break;
             }
        //}
    }

    if( indexRow != -1 && indexCol != -1 )
    {
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * @brief addObject Drop a object from the air with object's type specified by type
 * @param x The x position
 * @param z The z position
 * @param type The object's type
 * @param Height The height
 */
void GLWidget::addObject( const float x, const float z, const ObjectType type, const float y )
{
    switch( type )
    {
    case BOX:
    {
        const float length = 3.f;
        const float height = 3.f;
        const float width = 3.f;

        Object* newBox = new Box( m_fluid, Vector3(x,y,z), length, height, width, m_terrain->getdx(), m_boxTexID );
        newBox->initPhysics();
        m_objects.push_back( newBox );
        break;
    }
    case SPHERE:
    {
        const float radius = 2.f;
        Object* newSphere= new Sphere( m_fluid, Vector3(x,y,z), radius, m_terrain->getdx(), m_sphereTexID );
        newSphere->initPhysics();
        m_objects.push_back( newSphere );
        break;
    }
    default:
        assert(0);
        break;
    }
}

/**
 * @brief updateObjects Update the objects' positions
 * @param dt the time step
 */
void GLWidget::updateObjects( float dt )
{
    foreach( Object* o, m_objects )
    {
        if( o )
        {
            o->update( dt, m_fluid );
        }
    }
}

/**
 * @brief resetObjects Delete the objects
 */
void GLWidget::resetObjects()
{

    foreach( Object*o, m_objects )
    {
        if( o )
        {
            delete o;
            o = NULL;
        }
    }
    m_objects.clear();
}

/**
 * @brief loadConfig load the configuration from configure file
 * @return True if it load successfully
 */
 bool GLWidget::loadConfig( std::string fileName )
 {
    ifstream inFile;
    inFile.open( fileName.c_str() );
    if( !inFile.is_open() )
    {
        printf("Cannot find the configuration file:" CONFIG );
        return false;
    }
    string line;
    std::getline( inFile, line );
    bool correctHeader =!line.compare("#FLUID#");
    // Firstly check the header
    if( !correctHeader )
    {
        printf("Wrong type\n");
        inFile.close();
        return false;
    }
    while( !inFile.eof() )
    {
        std::string subLine;
        std::getline( inFile, subLine );
        boost::char_separator<char> sep("\t :");
        boost::tokenizer< boost::char_separator<char> > tokens(subLine, sep);
        vector<string> texts;
        BOOST_FOREACH (const std::string& t, tokens) {
            texts.push_back(t);
        }
        if( texts.size() <= 0 )
            continue;
        else
        {
            if( texts[0] == "#USE_HEIGHTMAP#" )
            {
                if( texts.size() != 3 )
                {
                    printf("Bad file\n");
                    return false;
                }
                else
                {
                    if( texts[1] == "1")
                    {
                        m_useHeightMap = true;
                        m_heightMapFileName = texts[2];
                    }
                    else if( texts[1] == "0" )
                    {
                        m_useHeightMap = false;
                    }
                    else
                    {
                        printf("Bad file\n");
                        return false;
                    }
                }
            }
            else if( texts[0] == "#TIME_STEP#")
            {
                if( texts.size() != 2 )
                {
                    printf("Bad file\n");
                    return false;
                }
                else
                {
                    std::stringstream ss;
                    ss<<texts[1];
                    ss>>m_timeStep;
                }
            }
            else
            {
                // Add more tokens
            }
        }
    }
    inFile.close();
    return true;
 }

 /**
  * @brief initConfig Init the configuration
  */
 void GLWidget::initConfig()
 {
     if( !loadConfig(CONFIG) )
     {
         m_timeStep = TIME_STEP;
         m_useHeightMap = false;
     }
 }

void GLWidget::tick()
{
    // Get the number of seconds since the last tick (variable update rate)
 //   float seconds = m_time.restart() * 0.001f;
    if( m_animate )
    {
        updateObjects( m_timeStep);
        timeUpdate();
#ifdef RENDER_FLUID
        m_fluid->update( m_timeStep );
#endif
    }
    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
