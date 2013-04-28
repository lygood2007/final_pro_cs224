#include "particle.h"

Particle::Particle(){
    m_radius = DEFAULT_PARTICLE_RADIUS;
    m_volume = DEFAULT_PARTICLE_VOLUME;

    m_position = Vector3::zero();
    m_velocity = Vector3::zero();
    m_acceleration = Vector3::zero();

    init();
}

Particle::Particle(double radius, double volume){
    m_radius = radius;
    m_volume = volume;

    m_position = Vector3::zero();
    m_velocity = Vector3::zero();
    m_acceleration = Vector3::zero();

    init();
}

Particle::Particle(double radius, double volume, Vector3 position, Vector3 velocity, Vector3 acceleration){
    m_radius = radius;
    m_volume = volume;

    m_position = position;
    m_velocity = velocity;
    m_acceleration = acceleration;

    init();
}

Particle::~Particle(){
    gluDeleteQuadric(m_quadric);
}

void Particle::init(){
    m_quadric = gluNewQuadric();
    m_color = Colorf(0.1f,0.4f,0.8f,1.0f);
}

Vector3 Particle::computeNextPosition(double dt){
    Vector3 newPosition = m_position;
    newPosition += dt * m_velocity;
    newPosition += dt * dt * m_acceleration;
    return newPosition;
}

Vector3 Particle::computeNextVelocity(double dt){
    Vector3 newVelocity = m_velocity;
    newVelocity += dt * m_acceleration;
    return newVelocity;
}

void Particle::updateParticle(double dt){
    m_position = computeNextPosition(dt);
    m_velocity = computeNextVelocity(dt);
}

void Particle::drawParticle(){
    //start
    glPushMatrix();
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE,GL_ONE);
    glColor4f(m_color.r, m_color.g, m_color.b, m_color.a);

    //translate to particle position
    glTranslatef(m_position.x, m_position.y, m_position.z);

    //draw sphere
    gluQuadricDrawStyle(m_quadric, GLU_FILL);
    gluSphere(m_quadric, m_radius, DEFAULT_PARTICLE_SLICES, DEFAULT_PARTICLE_SLICES);

    //end
    glDisable(GL_BLEND);
    glPopMatrix();
}
