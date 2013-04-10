#ifndef GLWIDGET_H
#define GLWIDGET_H

#include <QGLWidget>
#include <QMutex>
#include <vlc/vlc.h>
#include <vlc/libvlc.h>

class GLWidget : public QGLWidget
{
    Q_OBJECT
public:
    GLWidget(QWidget *parent);
    virtual ~GLWidget();

    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();

    void refreshtex();
    void updatetex(uchar *pixels);

    void mousePressEvent(QMouseEvent *);

private:
    unsigned int tex_;
    QImage *image_;

    void initVLC();

    static void *lock(void*, void**);
    static void unlock(void*, void*, void *const*);
    static void newframe(void*, void*);
    libvlc_instance_t*      m_vlcInstance;
    libvlc_media_t*         m_vlcMedia;
    libvlc_media_player_t*  m_vlcMediaplayer;

    bool texflag;

    void playvideo(QString path);
    
};

struct ctx
{
    uchar*                  pixels;
    QMutex*                 mutex;
    GLWidget*               parent;
};

#endif // GLWIDGET_H
