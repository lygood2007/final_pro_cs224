#include "particlesource.h"

ParticleSource::ParticleSource(){
    m_starting_corner = Vector3::zero();
    m_ending_corner = Vector3(1, 1, 1);
    m_num_particles_per_timestep = DEFAULT_NUM_PARTICLES_PER_TIMESTEP;
}

ParticleSource::ParticleSource(Vector3 startingCorner, Vector3 endingCorner){
    m_starting_corner = startingCorner;
    m_ending_corner = endingCorner;
    m_num_particles_per_timestep = DEFAULT_NUM_PARTICLES_PER_TIMESTEP;
}

ParticleSource::ParticleSource(Vector3 startingCorner, Vector3 endingCorner, int numParticlesPerTimestep){
    m_starting_corner = startingCorner;
    m_ending_corner = endingCorner;
    m_num_particles_per_timestep = numParticlesPerTimestep;
}

QVector<Particle*> ParticleSource::generateParticles(){
    //volume of the particles
    float Veff = C_DEPOSIT * (4 / 3) * M_PI *
            SPLASH_PARTICLE_RADIUS * SPLASH_PARTICLE_RADIUS * SPLASH_PARTICLE_RADIUS;

    QVector<Particle*> particles;
    for(int i = 0; i < m_num_particles_per_timestep; i++){
        float randX = randomFloatGenerator(m_starting_corner.x, m_ending_corner.x);
        float randY = randomFloatGenerator(m_starting_corner.y, m_ending_corner.y);
        float randZ = randomFloatGenerator(m_starting_corner.z, m_ending_corner.z);

        Vector3 position = Vector3(randX, randY, randZ);
        Vector3 velocity = Vector3::zero();
        Vector3 acceleration = Vector3(0, GRAVITY, 0);

        //make new particle
        Particle *newParticle = new Particle(SPLASH_PARTICLE_RADIUS, Veff, position, velocity, acceleration);
        particles.append(newParticle);
    }

    return particles;
}
