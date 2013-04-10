#include "glwidget.h"
#include <qgl.h>
#include <GL/glu.h>
#include <iostream>
#include <stdio.h>
#include <QLayout>

// CHANGE THESE ACCORDING TO YOUR DESIRED RESOLUTION
#define VIDEO_WIDTH 640
#define VIDEO_HEIGHT 480

using namespace std;
extern "C"
{
void invertImage(unsigned char *bits, int width, int height);
}
GLWidget::GLWidget(QWidget *parent) :
    QGLWidget(parent)
{
    texflag = false;
    initVLC();

    setMinimumSize(VIDEO_WIDTH, VIDEO_HEIGHT);
}

void GLWidget::initVLC()
{
    // List of parameters used to initialize libvlc.
    // These arguments are same as those you can pass
    // the the VLC command line.
    char const* vlc_argv[] =
    {
        //"--verbose", "3",
        // Edit this line if libvlc can't locate your plugins directory
        //"--plugin-path", "/path/to/vlc",
    };
    int vlc_argc = sizeof(vlc_argv) / sizeof(*vlc_argv);

    // Create a libvlc instance
    m_vlcInstance = libvlc_new( vlc_argc, vlc_argv);
    // Create the mediaplayer used to play a media
    m_vlcMediaplayer = libvlc_media_player_new( m_vlcInstance);

    // We're done with the initialization!
}

GLWidget::~GLWidget()
{
    delete image_;
}

void GLWidget::refreshtex()
{
    glGenTextures(1, &tex_);
    glBindTexture(GL_TEXTURE_2D, tex_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    QImage img = QGLWidget::convertToGLFormat(*image_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, img.width(), img.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, img.bits());
}

void GLWidget::initializeGL()
{
    glClearColor(0,0,0,0);
    image_ = new QImage(VIDEO_WIDTH, VIDEO_HEIGHT, QImage::Format_RGB32);

    refreshtex();

    glEnable(GL_TEXTURE_2D);
    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, width(), height());

    glOrtho(-1, 1, 1, -1, 0, 1);

// CHANGE THIS PATH TO YOUR VIDEO FILE
    QString path = \
            "plague.avi";
    playvideo(path);


}

void GLWidget::playvideo(QString path)
{
    // Create a new media from the path
    m_vlcMedia = libvlc_media_new_path( m_vlcInstance, path.toAscii());

    // We now need a struct for storing the video buffer
    // and a mutex to protect it.
    // The structure will be given as an arguments for the
    // lock/unlock callbacks.
    struct ctx* context;
    // Allocating the space for the structure
    context = ( struct ctx* )malloc( sizeof( *context ) );
    // Allocating the video buffer
    context->pixels = ( uchar* )malloc( ( sizeof( *( context->pixels ) ) * VIDEO_WIDTH * VIDEO_HEIGHT ) * 4 );
    // Allocating the mutex
    context->mutex = new QMutex();
    context->parent = this;

    libvlc_video_set_callbacks(m_vlcMediaplayer, lock, unlock, newframe, context);
    libvlc_video_set_format(m_vlcMediaplayer, "RV32", VIDEO_WIDTH, VIDEO_HEIGHT, 4*VIDEO_WIDTH);

    // Put the media into the mediaplayer
    libvlc_media_player_set_media( m_vlcMediaplayer, m_vlcMedia);

    // Finally, start the playback.
    libvlc_media_player_play( m_vlcMediaplayer);
}

void *GLWidget::lock(void *data, void **buffer)
{
    struct ctx *context = (ctx*) data;
    context->mutex->lock();
    *buffer = context->pixels;
}

void GLWidget::unlock(void *data, void *id, void *const *buffer)
{
    struct ctx *context = (ctx *) data;

    context->parent->updatetex(context->pixels);

    context->mutex->unlock();
}

void GLWidget::newframe(void *data, void *id)
{

}

void GLWidget::updatetex(uchar *pixels)
{
// CALL YOUR IMAGE PROCESSING FUNCTIONS HERE
invertImage(pixels, image_->width(), image_->height());
    memcpy(image_->bits(), pixels, 4*sizeof(unsigned char)*image_->width()*image_->height());
    texflag = true;
}

void GLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, width(), height());
    glOrtho(-1, 1, 1, -1, 0, 1);
}

void GLWidget::paintGL()
{
    // get the current frame
    if(texflag)
    {
        QImage img = QGLWidget::convertToGLFormat(*image_);
        glBindTexture(GL_TEXTURE_2D, tex_);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, img.width(), img.height(), GL_RGBA, GL_UNSIGNED_BYTE, img.bits());
        texflag = false;
    }
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBindTexture(GL_TEXTURE_2D, tex_);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(-1, -1);
    glTexCoord2f(0, 1);
    glVertex2f(-1, 1);
    glTexCoord2f(1, 1);
    glVertex2f(1, 1);
    glTexCoord2f(1, 0);
    glVertex2f(1, -1);
    glEnd();
    glBindTexture(GL_TEXTURE_2D, 0);
    //*/
    update();
}

void GLWidget::mousePressEvent(QMouseEvent *)
{

}
