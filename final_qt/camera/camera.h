/** camera.h
 ** Brief: This is the header file of the camera structure.
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef CAMERA_H
#define CAMERA_H
#include "vector.h"
#include <QMouseEvent>
#include "CS123Algebra.h"

class OrbitCamera
{

public:
    OrbitCamera();
    ~OrbitCamera();
    void mouseDown(const int x, const int y);
    void mouseMove(const int x, const int y);
    void mouseMovePan(const int x, const int y);
    void mouseWheel(const float delta);
    void updateMatrices();
    void updateProjectionMatrix();
    void updateModelviewMatrix();
    void translateCamera();
    void applyPerspectiveCamera(const float width, const float height);
    Matrix4x4 getInvViewTransMatrix();
    Vector4 getEyePos();
    void setRatio( float ratio );

//private: //I need all of these for rendering the lighting - SH
    Vector3 m_center, m_up;
    float m_theta, m_phi;
    float m_fovy;
    float m_zoom;

    float m_ratio;
    float m_near;
    float m_far;
    int m_oldX, m_oldY;

    Matrix4x4 m_projectionMatrix;
    Matrix4x4 m_modelviewMatrix;
};

#endif // CAMERA_H
