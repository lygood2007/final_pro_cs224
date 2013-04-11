/** camera.cpp
 ** Brief: This is the header file of the camera structure.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#include "GL/glut.h"
#include "CS123Common.h"
#include "camera.h"
#include <qgl.h>

OrbitCamera::OrbitCamera()
{
    m_up = Vector3(0.f, 1.f, 0.f);
    m_zoom = 30.f;
    m_theta = M_PI * 1.5f;
    m_phi = 1.2f;
    m_fovy = 60.f;
    m_near = 0.1f;
    m_far = 1000.0f;
    m_ratio = 1.0f;

    updateMatrices();
}

OrbitCamera::~OrbitCamera()
{
    //nothing to release
}

void OrbitCamera::mouseDown(const int x, const int y)
{
    m_oldX = x;
    m_oldY = y;
}

void OrbitCamera::mouseMove( int x, int y )
{
    Vector2 delta;
    delta.y = y - m_oldY;
    delta.x = x - m_oldX;
    m_oldX = x;
    m_oldY = y;

   m_theta += delta.x * 0.01f;
    m_phi += delta.y * 0.01f;

    m_theta -= floorf(m_theta / (2*M_PI)) * (2*M_PI);
    m_phi = max((float)(0.01f - M_PI / 2), min((float)(M_PI / 2 - 0.01f), m_phi));

    updateModelviewMatrix();
}

void OrbitCamera::mouseWheel(float delta)
{
    m_zoom *= powf(0.999f, delta);
    updateModelviewMatrix();
}

void OrbitCamera::updateMatrices()
{
    updateProjectionMatrix();
    updateModelviewMatrix();
}

void OrbitCamera::updateProjectionMatrix()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double matrix[16];
    glPushMatrix();
    glLoadIdentity();
    gluPerspective(m_fovy, m_ratio, m_near, m_far);
    glGetDoublev(GL_MODELVIEW_MATRIX, matrix);
    glPopMatrix();
    //because it's column order
    m_projectionMatrix = Matrix4x4(matrix).getTranspose();
}

void OrbitCamera::updateModelviewMatrix()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    double matrix[16];
    glPushMatrix();
    glLoadIdentity();

    // Move the object forward by m_zoomZ units before we rotate, so it will rotate about a point in front of us
    glTranslatef(0, 0, -m_zoom );
    // Now rotate the object, pivoting it about the new origin in front of us

    glRotatef(m_phi*180/M_PI, 1, 0, 0);
    glRotatef(m_theta*180/M_PI, 0, 1, 0);

    glGetDoublev(GL_MODELVIEW_MATRIX, matrix);
    glPopMatrix();
    // becaues it's column order
   m_modelviewMatrix = Matrix4x4(matrix).getTranspose();
}

void OrbitCamera::applyPerspectiveCamera( const float width, const float height)
{
    setRatio(width/height);
    Vector3 dir(-Vector3::fromAngles(m_theta, m_phi));
    Vector3 eye( - dir * m_zoom);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(m_fovy, m_ratio, m_near, m_far);
    gluLookAt(eye.x, eye.y, eye.z, eye.x + dir.x, eye.y + dir.y, eye.z + dir.z,
              m_up.x, m_up.y, m_up.z);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

Matrix4x4 OrbitCamera::getInvViewTransMatrix()
{

    Matrix4x4 invModelView = m_modelviewMatrix.getInverse();

    Matrix4x4 invScale = Matrix4x4::identity();

    REAL m_heightAngle = m_fovy;
    REAL heightRadians = m_heightAngle/360*M_PI;
    REAL h = m_far*tan(heightRadians);
    REAL w = h*m_ratio;

    invScale.data[0] = w;
    invScale.data[5] = h;
    invScale.data[10] = m_far;

    Matrix4x4 invViewTransMat = invModelView*invScale;
    return invViewTransMat;
}

Vector4 OrbitCamera::getEyePos()
{
    Vector4 eyePos;
    Vector3 dir2(-Vector3::fromAngles(m_theta, m_phi));

    updateModelviewMatrix();
    eyePos = m_modelviewMatrix.getInverse()*Vector4(0,0,0,1);
    return eyePos;
}

void OrbitCamera::setRatio(float ratio)
{
    m_ratio = ratio;
    updateProjectionMatrix();
}
