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
    m_center = Vector3(0.f, 0.f, 0.f); //added this and it might be wrong - SH
    m_up = Vector3(0.f, 1.f, 0.f);
    m_zoom = 30.f;
    m_theta = M_PI * 1.5f;
    m_phi = 1.2f;
    m_fovy = 60.f;
    m_near = 0.1f;
    m_far = 100000000.0f;
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


void OrbitCamera::mouseMovePan(const int x, const int y)
{

    float dx, dy;
    dy = y - m_oldY;
    dx = x - m_oldX;
    m_oldX = x;
    m_oldY = y;

    glMatrixMode(GL_MODELVIEW);
    double matrix[16];
    glPushMatrix();
    glLoadIdentity();

        //translate the camera by x and y
    glTranslatef(dx, dy, 0.f);

    glGetDoublev(GL_MODELVIEW_MATRIX, matrix);
    glPopMatrix();
    // because it's column order
   m_modelviewMatrix = Matrix4x4(matrix).getTranspose();
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
//    glLoadIdentity();
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
//    glLoadIdentity();
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

void OrbitCamera::translateCamera()
{
    setRatio(width/height);
    Vector3 dir(-Vector3::fromAngles(m_theta, m_phi));
    Vector3 eye( - dir * m_zoom);

    float dx, dy;
    dy = y - m_oldY;
    dx = x - m_oldX;
    m_oldX = x;
    m_oldY = y;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(m_fovy, m_ratio, m_near, m_far);
    //from here
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //moved this to here
    gluLookAt(dx, dy, eye.z, dx + dir.x, dy + dir.y, eye.z + dir.z,
              m_up.x, m_up.y, m_up.z); //dave said it's okay

}


void OrbitCamera::applyPerspectiveCamera( const float width, const float height)
{
    setRatio(width/height);
    Vector3 dir(-Vector3::fromAngles(m_theta, m_phi));
    Vector3 eye( - dir * m_zoom);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(m_fovy, m_ratio, m_near, m_far);
    //from here
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    //moved this to here
    gluLookAt(eye.x, eye.y, eye.z, eye.x + dir.x, eye.y + dir.y, eye.z + dir.z,
              m_up.x, m_up.y, m_up.z); //dave said it's okay
}

Matrix4x4 OrbitCamera::getInvViewTransMatrix()
{

    Matrix4x4 invModelView = m_modelviewMatrix.getInverse();

    Matrix4x4 invScale = Matrix4x4::identity();

    REAL m_heightAngle = m_fovy;
    REAL heightRadians = m_heightAngle/360*M_PI; //shouldn't this be 180??
    REAL h = m_far*tan(heightRadians);
    REAL w = h*m_ratio;

    invScale.data[0] = w;
    invScale.data[5] = h;
    invScale.data[10] = m_far;

    Matrix4x4 invViewTransMat = invModelView*invScale;
    return invViewTransMat;

    //These last two lines are what Dave thinks this entire method should be
    //however changing these breaks the fluid intersection

//    Matrix4x4 temp = m_modelviewMatrix * m_projectionMatrix;

//    return temp.getInverse();

}

Vector4 OrbitCamera::getEyePos()
{
    Vector4 eyePos;
//    Vector3 dir2(-Vector3::fromAngles(m_theta, m_phi)); //commented out, unused

//    updateModelviewMatrix(); //why are you updating in a get method?
    eyePos = m_modelviewMatrix.getInverse()*Vector4(0,0,0,1);
    return eyePos;
}

void OrbitCamera::setRatio(float ratio)
{
    m_ratio = ratio;
    updateProjectionMatrix();
}
