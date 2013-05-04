#include "particlesource.h"

ParticleSource::ParticleSource(){
    m_center_position = Vector3::zero();
    m_radius = DEFAULT_GENERATION_RADIUS;
    m_vertical_range = DEFAULT_VERTICAL_RANGE;
    m_num_particles_per_timestep = DEFAULT_NUM_PARTICLES_PER_TIMESTEP;
}

ParticleSource::ParticleSource(Vector3 centerPosition, float radius, float verticalRange, int numParticlesPerTimestep){
    m_center_position = centerPosition;
    m_radius = radius;
    m_vertical_range = verticalRange;
    m_num_particles_per_timestep = numParticlesPerTimestep;
}

void ParticleSource::generateParticles(Vector3 **particlePositions, Vector3 **particleVelocities, int totalNumParticles){
    Vector3 *currParticlePositions = (*particlePositions);
    Vector3 *currParticleVelocities = (*particleVelocities);

    int currIndex = 0;
    for(int i = 0; i < m_num_particles_per_timestep; i++){
        if(currIndex >= totalNumParticles){
            break;
        }

        //jitter particle positions
        //angle and radius of the cylinder's circle
        float angle = randomFloatGenerator(0.0f, 2.0f * M_PI);
        float radius = randomFloatGenerator(0.0f, m_radius);

        float randX = radius * cos(angle);
        float randY = randomFloatGenerator(0.0f, m_vertical_range);
        float randZ = radius * sin(angle);

        Vector3 position = Vector3(m_center_position.x + randX, m_center_position.y + randY, m_center_position.z + randZ);
        Vector3 velocity = Vector3::zero();

        //find new inactive particle
        while(currIndex < totalNumParticles){
            if(currParticlePositions[currIndex].y < TERRAIN_MIN_HEIGHT){
                //set new particle as active
                currParticlePositions[currIndex] = position;
                currParticleVelocities[currIndex] = velocity;

                currIndex++;
                break;
            }

            currIndex++;
        }
    }
}
