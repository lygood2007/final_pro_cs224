/** fluid_global.h
 ** Brief: This is the header file of all the global macros for fluid simulation.
 **           Originally it's in fluid.h. But now we have two version of fluid ( CPU and GPU ) so I set them here
 ** Project: large-scale fluids
 ** Date: 04/16/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef FLUID_GLOBAL_H
#define FLUID_GLOBAL_H

#define TERRAIN_BOUND 80
#define TERRAIN_BOTTOM_BOUND 0
#define TERRAIN_MIN_HEIGHT 0
#define TERRAIN_MAX_HEIGHT 40
#define DEFAULT_TERRAIN_DEPTH 8

#define GRID_SIZE TERRAIN_BOUND*2
#define DOMAIN_SIZE TERRAIN_BOUND/2
#define GRAVITY -10
#define SAVE_NAME_DEPTH "depth"
#define SAVE_NAME_VELOCITY "velocity"
//#define SAVE_IMAGE

//#define DAMPEN_WAVES
#define LAMBDA_DECAY 0.9
#define LAMBDA_UPDATE 0.1
#define DAMPENING_REGION 16
#define QUADRATIC_A 10.0f
#define QUADRATIC_B 0.0f
#define QUADRATIC_C 0.0f
#define INIT_SIGMA_GAMMA 0.0f
#define INIT_PHI_PSI 0.0f
#define CLAMP_ALPHA 0.5f

#define C_DEPOSIT 1
#define SPRAY_PARTICLE_RADIUS 0.05
#define SPLASH_PARTICLE_RADIUS 0.15
#define FOAM_PARTICLE_RADIUS 0.25
#define ALPHA_MIN_SPLASH 0.45
#define V_MIN_SPLASH 4
#define L_MIN_SPLASH -4
#define BREAKING_WAVE_NUM_SPLASH_PARTICLES 4
#define LAMBDA_Y 0.1

#define BREAKING_WAVE_VEL_MULTIPLIER 1.0f
#define SPRAY_VEL_MULTIPLIER 2.0f
#define FOAM_TTL 2.5f
#define FOAM_TTL_VARIANCE_MULTIPLIER 4.0f
#define FOAM_GENERATION_PROBABILITY 0.2f

#define NUM_DROPPING_PARTICLES 10000
#define PARTICLE_DROPPING_RADIUS 2.0
#define PARTICLE_DROP_HEIGHT 60
#define PARTICLE_DROP_RANGE 5

#define TOTAL_NUM_SPRAY_PARTICLES 240000
#define TOTAL_NUM_SPLASH_PARTICLES 120000
#define TOTAL_NUM_FOAM_PARTICLES 40000

#define USE_PARTICLE_SOURCES //uncomment to use particle sources //TODO: make this a hot key

//particle colors - @NOTE these awful colors are only temporary
#define SPRAY_COLOR 1.0f,0.f,0.f,0.9f
#define SPLASH_COLOR 1.0f,1.0f,1.0f,1.0f //1.0f,0.f,1.0f,0.9f
#define FOAM_COLOR 0.0f,1.0f,0.0f,0.9f

#define OBJECT_ORIGIN_HEIGHT TERRAIN_MAX_HEIGHT
#define WATER_DENSITY 1000


// If you don't want to render the volumn, comment this
#define RENDER_VOLUME

const float defaultHeight = TERRAIN_MAX_HEIGHT-10;
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
    PARTICLE_VELOCITIES,
    SPRAY_POSITIONS,
    SPRAY_VELOCITIES,
    FOAM_POSITIONS,
    FOAM_TTLS,
    SPLASH_TO_FOAM,
    PAINT,
    BREAKING_WAVES
};

enum DrawMethod
{
    DRAW_POINTS = 0,
    DRAW_MESH_STRIP,
    DRAW_MESH,
    DRAW_MESH_VBO
};

#endif // FLUID_GLOBAL_H
