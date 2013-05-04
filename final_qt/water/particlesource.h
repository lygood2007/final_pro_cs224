#ifndef PARTICLESOURCE_H
#define PARTICLESOURCE_H

#include <QVector>
#include "types.h"
#include "terrain.h"
#include "glwidget.h"

#define DEFAULT_GENERATION_RADIUS 1.0f
#define DEFAULT_VERTICAL_RANGE 1.0f
#define DEFAULT_NUM_PARTICLES_PER_TIMESTEP 10
#define DEFAULT_RECTANGLE_DISTANCE 1.0f

class ParticleSource
{
public:
    ParticleSource();
    ParticleSource(Vector3 centerPosition, float radius, float verticalRange, int numParticlesPerTimestep);
    ParticleSource(float minX, float maxX, float minY, float maxY, float minZ, float maxZ,
                   int numParticlesPerTimestep);
    ParticleSource(Vector3 centerPosition, float radius, float verticalRange,
                   float minX, float maxX, float minY, float maxY, float minZ, float maxZ,
                   int numParticlesPerTimestep);

    inline Vector3 getCenterPosition(){ return m_center_position; }
    inline float getRadius(){ return m_radius; }
    inline float getMinX(){ return m_min_x; }
    inline float getMaxX(){ return m_max_x; }
    inline float getMinY(){ return m_min_y; }
    inline float getMaxY(){ return m_max_y; }
    inline float getMinZ(){ return m_min_z; }
    inline float getMaxZ(){ return m_max_z; }
    inline float getVerticalRange(){ return m_vertical_range; }
    inline int getNumParticlesPerTimestep(){ return m_num_particles_per_timestep; }
    inline void setCenterPosition(Vector3 centerPosition){ m_center_position = centerPosition; }
    inline void setRadius(float radius){ m_radius = radius; }
    inline void setMinX(float minX){ m_min_x = minX; }
    inline void setMaxX(float maxX){ m_max_x = maxX; }
    inline void setMinY(float minY){ m_min_y = minY; }
    inline void setMaxY(float maxY){ m_max_y = maxY; }
    inline void setMinZ(float minZ){ m_min_z = minZ; }
    inline void setMaxZ(float maxZ){ m_max_z = maxZ; }
    inline void setVerticalRange(float verticalRange){ m_vertical_range = verticalRange; }
    inline void setNumParticlesPerTimestep(int numParticlesPerTimestep){ m_num_particles_per_timestep = numParticlesPerTimestep; }

    void generateParticles(Vector3 **particlePositions, Vector3 **particleVelocities, int totalNumParticles, bool circular);

private:
    Vector3 m_center_position;
    float m_radius;
    float m_vertical_range;

    float m_min_x;
    float m_max_x;
    float m_min_y;
    float m_max_y;
    float m_min_z;
    float m_max_z;

    int m_num_particles_per_timestep;
};

#endif // PARTICLESOURCE_H
