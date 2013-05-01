/** fluid_global.h
 ** Brief: This is the header file of all the global macros for fluid simulation.
 **           Originally it's in fluid.h. But now we have two version of fluid ( CPU and GPU ) so I set them here
 ** Project: large-scale fluids
 ** Date: 04/16/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef FLUID_GLOBAL_H
#define FLUID_GLOBAL_H

#define TERRAIN_BOUND 50
#define TERRAIN_MIN_HEIGHT -30
#define TERRAIN_MAX_HEIGHT 20
#define DEFAULT_TERRAIN_DEPTH 7

#define GRID_SIZE TERRAIN_BOUND*2
#define DOMAIN_SIZE TERRAIN_BOUND
#define GRAVITY -10
#define SAVE_NAME_DEPTH "depth"
#define SAVE_NAME_VELOCITY "velocity"
//#define SAVE_IMAGE

//#define DAMPEN_WAVES
#define LAMBDA_DECAY 0.9
#define LAMBDA_UPDATE 0.1
#define DAMPENING_REGION 5
#define QUADRATIC_A 10.0f
#define QUADRATIC_B 0.0f
#define QUADRATIC_C 0.0f
#define INIT_SIGMA_GAMMA 0.0f
#define INIT_PHI_PSI 0.0f
#define CLAMP_ALPHA 0.5f

#define USE_PARTICLES
#define C_DEPOSIT 1
#define SPLASH_PARTICLE_RADIUS 0.1
#define ALPHA_MIN_SPLASH 0.45
#define V_MIN_SPLASH 4
#define L_MIN_SPLASH -4
#define BREAKING_WAVE_NUM_SPLASH_PARTICLES 50
#define LAMBDA_Y 0.1

#define NUM_DROPPING_PARTICLES 4000
#define PARTICLE_DROPPING_RADIUS 2.0
#define PARTICLE_DROP_HEIGHT 60
#define PARTICLE_DROP_RANGE 20

//#define USE_PARTICLES_2
#define TOTAL_NUM_PARTICLES 50000

const float defaultHeight = TERRAIN_MAX_HEIGHT-20;
const float defaultU = 0.f;
const float defaultW = 0.f;
const float maxHeight = TERRAIN_MAX_HEIGHT+10;

enum FieldType
{
    DEPTH = 0,
    VEL_U,
    VEL_W,
    VEL,
    TMP,
    TERRAINH,
    SIGMA,
    GAMMA,
    PHI,
    PSI,
    NORMAL,
    HEIGHT,
    PARTICLE_POSITIONS,
    PARTICLE_VELOCITIES
};

enum DrawMethod
{
    DRAW_POINTS = 0,
    DRAW_MESH
};

#endif // FLUID_GLOBAL_H
