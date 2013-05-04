#include "particlesource.h"

ParticleSource::ParticleSource(){
    m_center_position = Vector3::zero();
    m_radius = DEFAULT_GENERATION_RADIUS;
    m_vertical_range = DEFAULT_VERTICAL_RANGE;

    m_min_x = -DEFAULT_RECTANGLE_DISTANCE;
    m_max_x = DEFAULT_RECTANGLE_DISTANCE;
    m_min_y = -DEFAULT_RECTANGLE_DISTANCE;
    m_max_y = DEFAULT_RECTANGLE_DISTANCE;
    m_min_z = -DEFAULT_RECTANGLE_DISTANCE;
    m_max_z = DEFAULT_RECTANGLE_DISTANCE;

    m_num_particles_per_timestep = DEFAULT_NUM_PARTICLES_PER_TIMESTEP;
}

ParticleSource::ParticleSource(Vector3 centerPosition, float radius, float verticalRange, int numParticlesPerTimestep){
    m_center_position = centerPosition;
    m_radius = radius;
    m_vertical_range = verticalRange;

    m_min_x = -DEFAULT_RECTANGLE_DISTANCE;
    m_max_x = DEFAULT_RECTANGLE_DISTANCE;
    m_min_y = -DEFAULT_RECTANGLE_DISTANCE;
    m_max_y = DEFAULT_RECTANGLE_DISTANCE;
    m_min_z = -DEFAULT_RECTANGLE_DISTANCE;
    m_max_z = DEFAULT_RECTANGLE_DISTANCE;

    m_num_particles_per_timestep = numParticlesPerTimestep;
}

ParticleSource::ParticleSource(float minX, float maxX, float minY, float maxY, float minZ, float maxZ, int numParticlesPerTimestep){
    m_center_position = Vector3::zero();
    m_radius = DEFAULT_GENERATION_RADIUS;
    m_vertical_range = DEFAULT_VERTICAL_RANGE;

    m_min_x = minX;
    m_max_x = maxX;
    m_min_y = minY;
    m_max_y = maxY;
    m_min_z = minZ;
    m_max_z = maxZ;

    m_num_particles_per_timestep = numParticlesPerTimestep;
}

ParticleSource::ParticleSource(Vector3 centerPosition, float radius, float verticalRange,
               float minX, float maxX, float minY, float maxY, float minZ, float maxZ,
               int numParticlesPerTimestep){
    m_center_position = centerPosition;
    m_radius = radius;
    m_vertical_range = verticalRange;

    m_min_x = minX;
    m_max_x = maxX;
    m_min_y = minY;
    m_max_y = maxY;
    m_min_z = minZ;
    m_max_z = maxZ;

    m_num_particles_per_timestep = numParticlesPerTimestep;
}

void ParticleSource::generateParticles(Vector3 **particlePositions, Vector3 **particleVelocities, int totalNumParticles, bool circular){
    Vector3 *currParticlePositions = (*particlePositions);
    Vector3 *currParticleVelocities = (*particleVelocities);

    int currIndex = 0;
    for(int i = 0; i < m_num_particles_per_timestep; i++){
        if(currIndex >= totalNumParticles){
            break;
        }

        //position and velocity
        Vector3 position = Vector3::zero();
        Vector3 velocity = Vector3::zero();

        if(circular){
            //jitter particle positions
            //angle and radius of the cylinder's circle
            float angle = randomFloatGenerator(0.0f, 2.0f * M_PI);
            float radius = randomFloatGenerator(0.0f, m_radius);

            float randX = radius * cos(angle);
            float randY = randomFloatGenerator(0.0f, m_vertical_range);
            float randZ = radius * sin(angle);

            position = Vector3(m_center_position.x + randX, m_center_position.y + randY, m_center_position.z + randZ);
        } else {
            //jitter particle positions
            //in the rectangular area
            float randX = randomFloatGenerator(m_min_x, m_max_x);
            float randY = randomFloatGenerator(m_min_y, m_max_y);
            float randZ = randomFloatGenerator(m_min_z, m_max_z);

            position = Vector3(randX, randY, randZ);
        }

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
