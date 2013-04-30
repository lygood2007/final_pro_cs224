#ifndef PARTICLESOURCE_H
#define PARTICLESOURCE_H

#include <QVector>
#include "types.h"
#include "terrain.h"
#include "glwidget.h"
#include "particle.h"

#define DEFAULT_NUM_PARTICLES_PER_TIMESTEP 10

class ParticleSource
{
public:
    ParticleSource();
    ParticleSource(Vector3 startingCorner, Vector3 endingCorner);
    ParticleSource(Vector3 startingCorner, Vector3 endingCorner, int numParticlesPerTimestep);

    inline Vector3 getStartingCorner(){ return m_starting_corner; }
    inline Vector3 getEndingCorner(){ return m_ending_corner; }
    inline int getNumParticlesPerTimestep(){ return m_num_particles_per_timestep; }
    inline void setStartingCorner(Vector3 startingCorner){ m_starting_corner = startingCorner; }
    inline void setEndingCorner(Vector3 endingCorner){ m_ending_corner = endingCorner; }
    inline void setNumParticlesPerTimestep(int numParticlesPerTimestep){ m_num_particles_per_timestep = numParticlesPerTimestep; }

    QVector<Particle*> generateParticles();

private:
    Vector3 m_starting_corner;
    Vector3 m_ending_corner;
    int m_num_particles_per_timestep;
};

#endif // PARTICLESOURCE_H
