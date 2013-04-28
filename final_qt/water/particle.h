#ifndef PARTICLE_H
#define PARTICLE_H

#include <GL/glut.h>
#include <GL/glu.h>

#include <QVector>
#include "types.h"
#include "terrain.h"
#include "glwidget.h"

#define DEFAULT_PARTICLE_RADIUS 0.1
#define DEFAULT_PARTICLE_VOLUME 1.0
#define DEFAULT_PARTICLE_SLICES 10

class Particle
{
public:
    Particle();
    Particle(double radius, double volume);
    Particle(double radius, double volume, Vector3 position, Vector3 velocity, Vector3 acceleration);

    ~Particle();

    void init();

    inline double getRadius(){ return m_radius; }
    inline double getVolume(){ return m_volume; }
    inline Vector3 getPosition(){ return m_position; }
    inline Vector3 getVelocity(){ return m_velocity; }
    inline Vector3 getAcceleration(){ return m_acceleration; }
    inline Colorf getColor(){ return m_color; }

    inline void setRadius(double radius){ m_radius = radius; }
    inline void setVolume(double volume){ m_volume = volume; }
    inline void setPosition(Vector3 position){ m_position = position; }
    inline void setVelocity(Vector3 velocity){ m_velocity = velocity; }
    inline void setAcceleration(Vector3 acceleration){ m_acceleration = acceleration; }
    inline void setColor(Colorf color){ m_color = color; }

    Vector3 computeNextPosition(double dt);
    Vector3 computeNextVelocity(double dt);
    void updateParticle(double dt);

    void drawParticle();

private:
    double m_radius; // particle radius
    double m_volume; // particle volume, V_eff in paper

    Vector3 m_position; // position vector
    Vector3 m_velocity; // velocity vector
    Vector3 m_acceleration; // acceleration vector

    GLUquadric *m_quadric; // quadric for rendering
    Colorf m_color; // color for rendering
};

#endif // PARTICLE_H
