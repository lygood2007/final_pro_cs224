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
#define TERRAIN_MIN_HEIGHT 0
#define TERRAIN_MAX_HEIGHT 50
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

const float defaultHeight = TERRAIN_MAX_HEIGHT- 25;
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
    HEIGHT
};

enum DrawMethod
{
    DRAW_POINTS = 0,
    DRAW_MESH
};

#endif // FLUID_GLOBAL_H
