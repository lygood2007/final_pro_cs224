#ifndef PARTICLESOURCE_H
#define PARTICLESOURCE_H

#include <QVector>
#include "types.h"
#include "terrain.h"
#include "glwidget.h"

#define DEFAULT_GENERATION_RADIUS 1.0f
#define DEFAULT_VERTICAL_RANGE 1.0f
#define DEFAULT_NUM_PARTICLES_PER_TIMESTEP 10

class ParticleSource
{
public:
    ParticleSource();
    ParticleSource(Vector3 centerPosition, float radius, float verticalRange, int numParticlesPerTimestep);

    inline Vector3 getCenterPosition(){ return m_center_position; }
    inline float getRadius(){ return m_radius; }
    inline float getVerticalRange(){ return m_vertical_range; }
    inline int getNumParticlesPerTimestep(){ return m_num_particles_per_timestep; }
    inline void setCenterPosition(Vector3 centerPosition){ m_center_position = centerPosition; }
    inline void setRadius(float radius){ m_radius = radius; }
    inline void setVerticalRange(float verticalRange){ m_vertical_range = verticalRange; }
    inline void setNumParticlesPerTimestep(int numParticlesPerTimestep){ m_num_particles_per_timestep = numParticlesPerTimestep; }

    void generateParticles(Vector3 **particlePositions, Vector3 **particleVelocities, int totalNumParticles);

private:
    Vector3 m_center_position;
    float m_radius;
    float m_vertical_range;
    int m_num_particles_per_timestep;
};

#endif // PARTICLESOURCE_H
