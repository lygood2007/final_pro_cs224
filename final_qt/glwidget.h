#ifndef VIEW_H
#define VIEW_H

#include <qgl.h>
#include <QTime>
#include <QTimer>

class GLWidget : public QGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent);
    ~GLWidget();

private:
    QTime time;
    QTimer timer;

    void initializeGL();
    void paintGL();
    void resizeGL(int w, int h);

    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);

    void keyPressEvent(QKeyEvent *event);
    void keyReleaseEvent(QKeyEvent *event);

private slots:
    void tick();
};

#endif // VIEW_H

